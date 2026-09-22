#!/usr/bin/env python3
"""gCastle DAG-GNN with a hard exogenous-root mask.

Matrix convention:
    A[i, j] != 0 means node i -> node j.
The mask fixes A[:, root] to zero in the encoder, decoder, reconstruction loss,
sparsity loss, and acyclicity constraint throughout training.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("CASTLE_BACKEND", "pytorch")

import numpy as np
import pandas as pd
import torch
import torch.optim as optim
from torch.optim import lr_scheduler

from castle.algorithms.gradient.dag_gnn.torch.dag_gnn import DAG_GNN, set_seed
from castle.algorithms.gradient.dag_gnn.torch.models.modules import Decoder, Encoder
from castle.algorithms.gradient.dag_gnn.torch.utils import functions as func
from castle.algorithms.gradient.notears.torch.golem_utils.train import postprocess
from castle.common import Tensor


class RootMaskedEncoder(Encoder):
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        adj_A: torch.Tensor,
        allowed_mask: np.ndarray,
        device: torch.device,
        encoder_type: str = "mlp",
    ) -> None:
        super().__init__(
            input_dim=input_dim,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            adj_A=adj_A,
            device=device,
            encoder_type=encoder_type,
        )
        mask = torch.as_tensor(
            allowed_mask,
            dtype=self.adj_A.dtype,
            device=self.adj_A.device,
        )
        self.register_buffer("allowed_mask", mask)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        adj_A1 = torch.sinh(3.0 * self.adj_A) * self.allowed_mask
        adj_Aforz = torch.eye(
            adj_A1.shape[0],
            dtype=adj_A1.dtype,
            device=self.device,
        ) - adj_A1.T

        if self.encoder_type == "mlp":
            mlp_out = self.mlp(x.to(self.mlp.device))
            logits = torch.matmul(adj_Aforz, mlp_out + self.wa) - self.wa
        else:
            adj_A_inv = torch.inverse(adj_Aforz)
            meanF = torch.matmul(
                adj_A_inv,
                torch.mean(torch.matmul(adj_Aforz, x), 0),
            )
            logits = torch.matmul(adj_Aforz, x - meanF)

        return logits, adj_A1


class RootMaskedDAGGNN(DAG_GNN):
    def __init__(self, allowed_mask: np.ndarray, **kwargs: Any) -> None:
        super().__init__(**kwargs)
        self.allowed_mask = np.asarray(allowed_mask, dtype=np.float64)
        self.weight_causal_matrix = None
        self.weight_dag_matrix = None
        self.h_value = None

    def learn(self, data: np.ndarray, columns=None, **kwargs: Any) -> None:
        set_seed(self.seed)
        data = np.asarray(data, dtype=np.float64)

        if data.ndim == 2:
            data = np.expand_dims(data, axis=2)

        self.n_samples, self.n_nodes, self.input_dim = data.shape
        if self.allowed_mask.shape != (self.n_nodes, self.n_nodes):
            raise ValueError(
                f"allowed_mask has shape {self.allowed_mask.shape}, expected "
                f"{(self.n_nodes, self.n_nodes)}"
            )

        if self.latent_dim is None:
            self.latent_dim = self.input_dim

        train_loader = func.get_dataloader(
            data,
            batch_size=self.batch_size,
            device=self.device,
        )

        adj_A = torch.zeros(
            (self.n_nodes, self.n_nodes),
            dtype=torch.float64,
            requires_grad=True,
            device=self.device,
        )

        self.encoder = RootMaskedEncoder(
            input_dim=self.input_dim,
            hidden_dim=self.encoder_hidden,
            output_dim=self.latent_dim,
            adj_A=adj_A,
            allowed_mask=self.allowed_mask,
            device=self.device,
            encoder_type=self.encoder_type.lower(),
        ).double()

        self.decoder = Decoder(
            input_dim=self.latent_dim,
            hidden_dim=self.decoder_hidden,
            output_dim=self.input_dim,
            device=self.device,
            decoder_type=self.decoder_type.lower(),
        ).double()

        if self.optimizer.lower() == "adam":
            optimizer = optim.Adam(
                [
                    {"params": self.encoder.parameters()},
                    {"params": self.decoder.parameters()},
                ],
                lr=self.lr,
            )
        elif self.optimizer.lower() == "sgd":
            optimizer = optim.SGD(
                [
                    {"params": self.encoder.parameters()},
                    {"params": self.decoder.parameters()},
                ],
                lr=self.lr,
            )
        else:
            raise ValueError("optimizer must be 'adam' or 'sgd'.")

        self.scheduler = lr_scheduler.StepLR(
            optimizer,
            step_size=self.lr_decay,
            gamma=self.gamma,
        )

        c_a = self.init_c_a
        lambda_a = self.init_lambda_a
        h_a_new = torch.tensor(1.0, dtype=torch.float64, device=self.device)
        h_a_old = np.inf
        best_elbo_loss = np.inf
        origin_a = adj_A
        epoch = 0

        for step_k in range(self.k_max_iter):
            while c_a < self.c_a_thresh:
                for epoch in range(self.epochs):
                    elbo_loss, origin_a = self._train(
                        train_loader=train_loader,
                        optimizer=optimizer,
                        lambda_a=lambda_a,
                        c_a=c_a,
                    )
                    if elbo_loss < best_elbo_loss:
                        best_elbo_loss = elbo_loss

                if elbo_loss > 2.0 * best_elbo_loss:
                    break

                a_new = origin_a.detach().clone()
                h_a_new = func._h_A(a_new, self.n_nodes)

                if h_a_new.item() > self.multiply_h * h_a_old:
                    c_a *= self.eta
                else:
                    break

            h_a_old = h_a_new.item()
            logging.info("Iter: %s, epoch: %s, h_new: %s", step_k, epoch, h_a_old)
            lambda_a += c_a * h_a_new.item()

            if h_a_old <= self.h_tolerance:
                break

        W_raw = origin_a.detach().cpu().numpy().astype(float)
        W_raw *= self.allowed_mask
        np.fill_diagonal(W_raw, 0.0)

        W_dag = postprocess(W_raw, graph_thres=self.graph_threshold)
        W_dag *= self.allowed_mask

        binary = (W_dag != 0).astype(int)
        self.weight_causal_matrix = Tensor(W_raw, index=columns, columns=columns)
        self.weight_dag_matrix = Tensor(W_dag, index=columns, columns=columns)
        self.causal_matrix = Tensor(binary, index=columns, columns=columns)
        self.h_value = float(func._h_A(torch.as_tensor(W_raw), self.n_nodes).item())


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
    X = data.to_numpy(dtype=np.float64)
    means = X.mean(axis=0)
    scales = X.std(axis=0, ddof=0)
    zero_variance = np.flatnonzero(scales == 0)
    if zero_variance.size:
        names = [str(data.columns[i]) for i in zero_variance]
        raise ValueError(f"Zero-variance nodes: {', '.join(names)}")
    X = (X - means) / scales
    return X, {"mean": means.tolist(), "scale": scales.tolist()}


def build_allowed_mask(n_nodes: int, root_index: int | None) -> np.ndarray:
    mask = np.ones((n_nodes, n_nodes), dtype=np.float64)
    np.fill_diagonal(mask, 0.0)
    if root_index is not None:
        mask[:, root_index] = 0.0
    return mask


def save_outputs(
    output_dir: str | Path,
    columns: list[str],
    model: RootMaskedDAGGNN,
    metadata: dict[str, Any],
) -> None:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    W_raw = np.asarray(model.weight_causal_matrix, dtype=float)
    W_dag = np.asarray(model.weight_dag_matrix, dtype=float)
    binary = np.asarray(model.causal_matrix, dtype=int)

    pd.DataFrame(W_raw, index=columns, columns=columns).to_csv(
        output_dir / "dag_gnn_weighted_raw.csv",
        index_label="from",
    )
    pd.DataFrame(W_dag, index=columns, columns=columns).to_csv(
        output_dir / "dag_gnn_weighted_dag.csv",
        index_label="from",
    )
    pd.DataFrame(binary, index=columns, columns=columns).to_csv(
        output_dir / "dag_gnn_binary.csv",
        index_label="from",
    )

    with open(output_dir / "dag_gnn_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--root", default="A", help="Use 'none' to disable the root constraint.")
    parser.add_argument("--encoder-type", choices=("mlp", "sem"), default="mlp")
    parser.add_argument("--decoder-type", choices=("mlp", "sem"), default="mlp")
    parser.add_argument("--encoder-hidden", type=int, default=64)
    parser.add_argument("--decoder-hidden", type=int, default=64)
    parser.add_argument("--epochs", type=int, default=300)
    parser.add_argument("--k-max-iter", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--lr", type=float, default=3e-3)
    parser.add_argument("--lr-decay", type=int, default=200)
    parser.add_argument("--gamma", type=float, default=1.0)
    parser.add_argument("--tau-a", type=float, default=0.0)
    parser.add_argument("--h-tolerance", type=float, default=1e-8)
    parser.add_argument("--threshold", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--device-id", default="0")
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
        X = data.to_numpy(dtype=np.float64)
        scaling = None
    else:
        X, scaling = standardize_columns(data)

    allowed_mask = build_allowed_mask(X.shape[1], root_index)

    model = RootMaskedDAGGNN(
        allowed_mask=allowed_mask,
        encoder_type=args.encoder_type,
        decoder_type=args.decoder_type,
        encoder_hidden=args.encoder_hidden,
        decoder_hidden=args.decoder_hidden,
        epochs=args.epochs,
        k_max_iter=args.k_max_iter,
        tau_a=args.tau_a,
        batch_size=min(args.batch_size, X.shape[0]),
        lr=args.lr,
        lr_decay=args.lr_decay,
        gamma=args.gamma,
        h_tolerance=args.h_tolerance,
        graph_threshold=args.threshold,
        optimizer="adam",
        seed=args.seed,
        device_type=args.device,
        device_ids=args.device_id,
    )
    model.learn(X, columns=columns)

    W_dag = np.asarray(model.weight_dag_matrix, dtype=float)
    if root_index is not None and np.any(W_dag[:, root_index] != 0):
        raise RuntimeError("DAG-GNN output violates the root constraint.")

    metadata = {
        "method": "DAG-GNN",
        "implementation": "gCastle 1.0.4 with a hard encoder adjacency mask",
        "matrix_convention": "row i, column j denotes i -> j",
        "input": str(Path(args.input).resolve()),
        "n_samples": int(X.shape[0]),
        "n_nodes": int(X.shape[1]),
        "columns": columns,
        "root": root,
        "root_constraint": root is not None,
        "standardized": not args.no_standardize,
        "scaling": scaling,
        "encoder_type": args.encoder_type,
        "decoder_type": args.decoder_type,
        "encoder_hidden": args.encoder_hidden,
        "decoder_hidden": args.decoder_hidden,
        "epochs": args.epochs,
        "k_max_iter": args.k_max_iter,
        "batch_size": min(args.batch_size, X.shape[0]),
        "lr": args.lr,
        "lr_decay": args.lr_decay,
        "gamma": args.gamma,
        "tau_a": args.tau_a,
        "h_tolerance": args.h_tolerance,
        "threshold": args.threshold,
        "seed": args.seed,
        "device": args.device,
        "h_final": model.h_value,
        "n_edges": int(np.count_nonzero(W_dag)),
    }
    save_outputs(args.output_dir, columns, model, metadata)

    print(
        f"DAG-GNN completed: nodes={X.shape[1]}, edges={np.count_nonzero(W_dag)}, "
        f"h={model.h_value:.3e}, root={root}"
    )


if __name__ == "__main__":
    main()
