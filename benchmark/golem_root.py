#!/usr/bin/env python3
"""gCastle GOLEM with a hard exogenous-root mask.

Matrix convention:
    B[i, j] != 0 means node i -> node j.
The allowed mask fixes B[:, root] to zero in every forward pass, so the root
constraint is part of optimization rather than post-processing.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("CASTLE_BACKEND", "pytorch")

import numpy as np
import pandas as pd
import torch

from castle.algorithms.gradient.notears.torch.golem_utils.golem_model import GolemModel
from castle.algorithms.gradient.notears.torch.golem_utils.train import postprocess
from castle.algorithms.gradient.notears.torch.golem_utils.utils import set_seed


class RootMaskedGolemModel(GolemModel):
    def __init__(
        self,
        n: int,
        d: int,
        lambda_1: float,
        lambda_2: float,
        equal_variances: bool,
        allowed_mask: np.ndarray,
        B_init: np.ndarray | None = None,
        device: torch.device | None = None,
    ) -> None:
        if B_init is not None:
            B_init = np.asarray(
                B_init,
                dtype=np.float32
            )

        super().__init__(
            n=n,
            d=d,
            lambda_1=lambda_1,
            lambda_2=lambda_2,
            equal_variances=equal_variances,
            B_init=B_init,
            device=device,
        )

        # Explicitly force the learnable adjacency matrix to float32.
        # This avoids dependence on PyTorch's global default dtype.
        self._B = torch.nn.Parameter(
            self._B.detach().to(
                dtype=torch.float32,
                device=self._B.device
            )
        )

        mask = torch.as_tensor(
            allowed_mask,
            dtype=torch.float32,
            device=self._B.device,
        )

        self.register_buffer(
            "allowed_mask",
            mask
        )

    def _preprocess(
        self,
        B: torch.Tensor
    ) -> torch.Tensor:
        # allowed_mask already contains:
        # 1. zero diagonal;
        # 2. zero root column.
        #
        # Do not call super()._preprocess(), because the gCastle
        # implementation creates torch.ones()/torch.eye() without
        # specifying dtype, which can promote B to float64.
        return B * self.allowed_mask



def load_numeric_csv(path: str | Path) -> pd.DataFrame:
    data = pd.read_csv(path)
    unnamed = [column for column in data.columns if str(column).startswith("Unnamed:")]
    if unnamed:
        data = data.drop(columns=unnamed)

    numeric = data.apply(pd.to_numeric, errors="coerce")
    if numeric.shape[1] < 2:
        raise ValueError("The input matrix must contain at least two node columns.")
    if numeric.isna().any().any():
        bad = data.columns[numeric.isna().any(axis=0)].tolist()
        raise ValueError(f"Non-numeric or missing values were found in: {bad}")
    if not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError("The input matrix contains NaN or infinite values.")
    return numeric


def standardize_columns(data: pd.DataFrame) -> tuple[np.ndarray, dict[str, list[float]]]:
    X = data.to_numpy(dtype=np.float32)
    means = X.mean(axis=0)
    scales = X.std(axis=0, ddof=0)
    zero_variance = np.flatnonzero(scales == 0)
    if zero_variance.size:
        names = [str(data.columns[i]) for i in zero_variance]
        raise ValueError(f"Zero-variance nodes: {', '.join(names)}")
    X = (X - means) / scales
    return X.astype(np.float32), {"mean": means.tolist(), "scale": scales.tolist()}


def build_allowed_mask(n_nodes: int, root_index: int | None) -> np.ndarray:
    mask = np.ones((n_nodes, n_nodes), dtype=np.float32)
    np.fill_diagonal(mask, 0.0)
    if root_index is not None:
        mask[:, root_index] = 0.0
    return mask


def train_golem_stage(
    X: np.ndarray,
    allowed_mask: np.ndarray,
    lambda_1: float,
    lambda_2: float,
    equal_variances: bool,
    learning_rate: float,
    num_iter: int,
    checkpoint_iter: int | None,
    seed: int,
    device: torch.device,
    B_init: np.ndarray | None = None,
) -> tuple[np.ndarray, dict[str, Any]]:
    set_seed(seed)

    n, d = X.shape

    model = RootMaskedGolemModel(
        n=n,
        d=d,
        lambda_1=lambda_1,
        lambda_2=lambda_2,
        equal_variances=equal_variances,
        allowed_mask=allowed_mask,
        B_init=B_init,
        device=device,
    )

    # Match the input tensor dtype and device to the model parameter.
    X_tensor = torch.as_tensor(
        X,
        dtype=model._B.dtype,
        device=model._B.device
    )

    if X_tensor.dtype != model._B.dtype:
        raise RuntimeError(
            "GOLEM dtype mismatch: "
            f"X={X_tensor.dtype}, B={model._B.dtype}"
        )

    optimizer = torch.optim.Adam(
        model.parameters(),
        lr=learning_rate
    )

    history: list[dict[str, float]] = []

    for iteration in range(int(num_iter) + 1):
        model(X_tensor)

        if iteration > 0:
            optimizer.zero_grad()
            model.score.backward()
            optimizer.step()

        if (
            checkpoint_iter is not None
            and iteration % checkpoint_iter == 0
        ):
            history.append(
                {
                    "iteration": int(iteration),
                    "score": float(
                        model.score.detach().cpu()
                    ),
                    "likelihood": float(
                        model.likelihood.detach().cpu()
                    ),
                    "h": float(
                        model.h.detach().cpu()
                    ),
                }
            )

    model(X_tensor)

    W_raw = (
        model.B
        .detach()
        .cpu()
        .numpy()
        .astype(float)
    )

    W_raw *= allowed_mask
    np.fill_diagonal(W_raw, 0.0)

    diagnostics = {
        "score_final": float(
            model.score.detach().cpu()
        ),
        "likelihood_final": float(
            model.likelihood.detach().cpu()
        ),
        "h_final": float(
            model.h.detach().cpu()
        ),
        "dtype": str(model._B.dtype),
        "history": history,
    }

    return W_raw, diagnostics


def save_outputs(
    output_dir: str | Path,
    columns: list[str],
    W_raw: np.ndarray,
    W_dag: np.ndarray,
    metadata: dict[str, Any],
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(W_raw, index=columns, columns=columns).to_csv(
        output_dir / "golem_weighted_raw.csv",
        index_label="from",
    )
    pd.DataFrame(W_dag, index=columns, columns=columns).to_csv(
        output_dir / "golem_weighted_dag.csv",
        index_label="from",
    )
    pd.DataFrame((W_dag != 0).astype(int), index=columns, columns=columns).to_csv(
        output_dir / "golem_binary.csv",
        index_label="from",
    )

    with open(output_dir / "golem_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--root", default="A", help="Use 'none' to disable the root constraint.")
    parser.add_argument("--mode", choices=("ev", "nv", "two-stage"), default="two-stage")
    parser.add_argument("--lambda1-ev", type=float, default=2e-2)
    parser.add_argument("--lambda1-nv", type=float, default=2e-3)
    parser.add_argument("--lambda2", type=float, default=5.0)
    parser.add_argument("--learning-rate", type=float, default=1e-3)
    parser.add_argument("--num-iter-ev", type=int, default=20000)
    parser.add_argument("--num-iter-nv", type=int, default=20000)
    parser.add_argument("--checkpoint-iter", type=int, default=5000)
    parser.add_argument("--threshold", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--no-standardize", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    data = load_numeric_csv(args.input)
    columns = [str(column) for column in data.columns]

    root = None if args.root.lower() == "none" else args.root
    if root is not None and root not in columns:
        raise ValueError(f"Root node '{root}' is not present in the input columns.")
    root_index = None if root is None else columns.index(root)

    if args.no_standardize:
        X = data.to_numpy(dtype=np.float32)
        scaling = None
    else:
        X, scaling = standardize_columns(data)

    if args.device == "gpu":
        if not torch.cuda.is_available():
            raise RuntimeError("GPU was requested, but torch.cuda.is_available() is FALSE.")
        device = torch.device("cuda")
    else:
        device = torch.device("cpu")

    allowed_mask = build_allowed_mask(X.shape[1], root_index)
    checkpoint_iter = None if args.checkpoint_iter <= 0 else args.checkpoint_iter
    stage_diagnostics: dict[str, Any] = {}

    if args.mode == "ev":
        W_raw, stage_diagnostics["ev"] = train_golem_stage(
            X=X,
            allowed_mask=allowed_mask,
            lambda_1=args.lambda1_ev,
            lambda_2=args.lambda2,
            equal_variances=True,
            learning_rate=args.learning_rate,
            num_iter=args.num_iter_ev,
            checkpoint_iter=checkpoint_iter,
            seed=args.seed,
            device=device,
        )
    elif args.mode == "nv":
        W_raw, stage_diagnostics["nv"] = train_golem_stage(
            X=X,
            allowed_mask=allowed_mask,
            lambda_1=args.lambda1_nv,
            lambda_2=args.lambda2,
            equal_variances=False,
            learning_rate=args.learning_rate,
            num_iter=args.num_iter_nv,
            checkpoint_iter=checkpoint_iter,
            seed=args.seed,
            device=device,
        )
    else:
        W_ev, stage_diagnostics["ev"] = train_golem_stage(
            X=X,
            allowed_mask=allowed_mask,
            lambda_1=args.lambda1_ev,
            lambda_2=args.lambda2,
            equal_variances=True,
            learning_rate=args.learning_rate,
            num_iter=args.num_iter_ev,
            checkpoint_iter=checkpoint_iter,
            seed=args.seed,
            device=device,
        )
        W_raw, stage_diagnostics["nv"] = train_golem_stage(
            X=X,
            allowed_mask=allowed_mask,
            lambda_1=args.lambda1_nv,
            lambda_2=args.lambda2,
            equal_variances=False,
            learning_rate=args.learning_rate,
            num_iter=args.num_iter_nv,
            checkpoint_iter=checkpoint_iter,
            seed=args.seed,
            device=device,
            B_init=W_ev,
        )

    W_raw *= allowed_mask
    W_dag = postprocess(W_raw, graph_thres=args.threshold)
    W_dag *= allowed_mask

    if root_index is not None and np.any(W_dag[:, root_index] != 0):
        raise RuntimeError("GOLEM output violates the root constraint.")

    metadata = {
        "method": "GOLEM",
        "implementation": "gCastle 1.0.4 model with a hard allowed-edge mask",
        "matrix_convention": "row i, column j denotes i -> j",
        "input": str(Path(args.input).resolve()),
        "n_samples": int(X.shape[0]),
        "n_nodes": int(X.shape[1]),
        "columns": columns,
        "root": root,
        "root_constraint": root is not None,
        "mode": args.mode,
        "standardized": not args.no_standardize,
        "scaling": scaling,
        "lambda1_ev": args.lambda1_ev,
        "lambda1_nv": args.lambda1_nv,
        "lambda2": args.lambda2,
        "learning_rate": args.learning_rate,
        "num_iter_ev": args.num_iter_ev,
        "num_iter_nv": args.num_iter_nv,
        "threshold": args.threshold,
        "seed": args.seed,
        "device": str(device),
        "n_edges": int(np.count_nonzero(W_dag)),
        "diagnostics": stage_diagnostics,
    }
    save_outputs(args.output_dir, columns, W_raw, W_dag, metadata)

    final_stage = "nv" if "nv" in stage_diagnostics else "ev"
    print(
        f"GOLEM completed: mode={args.mode}, nodes={X.shape[1]}, "
        f"edges={np.count_nonzero(W_dag)}, "
        f"h={stage_diagnostics[final_stage]['h_final']:.3e}, root={root}"
    )


if __name__ == "__main__":
    main()
