#!/usr/bin/env python3
"""Linear NOTEARS with a hard exogenous-root constraint.

Matrix convention throughout this script:
    W[i, j] != 0 means node i -> node j.
Therefore all incoming edges into root r are forbidden by fixing W[:, r] = 0
throughout optimization, not by deleting edges after fitting.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import pandas as pd
import scipy.linalg as slin
import scipy.optimize as sopt
from scipy.special import expit as sigmoid


def load_numeric_csv(path: str | Path) -> pd.DataFrame:
    data = pd.read_csv(path)

    unnamed = [column for column in data.columns if str(column).startswith("Unnamed:")]
    if unnamed:
        data = data.drop(columns=unnamed)

    if data.shape[1] < 2:
        raise ValueError("The input matrix must contain at least two node columns.")

    numeric = data.apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any():
        bad = data.columns[numeric.isna().any(axis=0)].tolist()
        raise ValueError(f"Non-numeric or missing values were found in: {bad}")

    values = numeric.to_numpy(dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("The input matrix contains NaN or infinite values.")

    return numeric


def standardize_columns(data: pd.DataFrame) -> tuple[np.ndarray, dict[str, list[float]]]:
    X = data.to_numpy(dtype=float)
    means = X.mean(axis=0)
    scales = X.std(axis=0, ddof=0)

    zero_variance = np.flatnonzero(scales == 0)
    if zero_variance.size:
        names = [str(data.columns[i]) for i in zero_variance]
        raise ValueError(f"Zero-variance nodes: {', '.join(names)}")

    X = (X - means) / scales
    return X, {"mean": means.tolist(), "scale": scales.tolist()}


def is_dag(weight_matrix: np.ndarray) -> bool:
    adjacency = np.asarray(weight_matrix != 0, dtype=int)
    np.fill_diagonal(adjacency, 0)
    return nx.is_directed_acyclic_graph(nx.DiGraph(adjacency))


def postprocess_to_dag(
    weight_matrix: np.ndarray,
    threshold: float,
) -> tuple[np.ndarray, float]:
    """Threshold weights and remove weakest remaining edges until a DAG is obtained."""
    W = np.array(weight_matrix, dtype=float, copy=True)
    W[np.abs(W) <= threshold] = 0.0
    np.fill_diagonal(W, 0.0)

    if is_dag(W):
        return W, threshold

    rows, cols = np.nonzero(W)
    weighted_edges = sorted(
        ((abs(W[i, j]), i, j) for i, j in zip(rows, cols)),
        key=lambda item: item[0],
    )

    dag_threshold = threshold
    for magnitude, i, j in weighted_edges:
        if is_dag(W):
            break
        W[i, j] = 0.0
        dag_threshold = max(dag_threshold, float(magnitude))

    if not is_dag(W):
        raise RuntimeError("Failed to obtain a DAG during post-processing.")

    return W, dag_threshold


def notears_linear_root(
    X: np.ndarray,
    lambda1: float,
    loss_type: str,
    root_index: int | None,
    max_iter: int = 100,
    h_tol: float = 1e-8,
    rho_max: float = 1e16,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Solve linear NOTEARS with all incoming root coefficients fixed to zero."""
    X = np.asarray(X, dtype=float)
    n, d = X.shape

    if root_index is not None and not 0 <= root_index < d:
        raise ValueError("root_index is outside the node range.")

    if loss_type == "l2":
        X = X - np.mean(X, axis=0, keepdims=True)

    def _loss(W: np.ndarray) -> tuple[float, np.ndarray]:
        M = X @ W

        if loss_type == "l2":
            R = X - M
            loss = 0.5 / n * (R**2).sum()
            gradient = -1.0 / n * X.T @ R
        elif loss_type == "logistic":
            loss = 1.0 / n * (np.logaddexp(0, M) - X * M).sum()
            gradient = 1.0 / n * X.T @ (sigmoid(M) - X)
        elif loss_type == "poisson":
            S = np.exp(M)
            loss = 1.0 / n * (S - X * M).sum()
            gradient = 1.0 / n * X.T @ (S - X)
        else:
            raise ValueError("loss_type must be one of: l2, logistic, poisson")

        return float(loss), gradient

    def _h(W: np.ndarray) -> tuple[float, np.ndarray]:
        E = slin.expm(W * W)
        h_value = np.trace(E) - d
        gradient = E.T * W * 2.0
        return float(h_value), gradient

    def _adj(w: np.ndarray) -> np.ndarray:
        return (w[: d * d] - w[d * d :]).reshape(d, d)

    rho = 1.0
    alpha = 0.0
    h_value = np.inf
    w_est = np.zeros(2 * d * d, dtype=float)
    history: list[dict[str, Any]] = []

    bounds = [
        (0.0, 0.0)
        if i == j or (root_index is not None and j == root_index)
        else (0.0, None)
        for _ in range(2)
        for i in range(d)
        for j in range(d)
    ]

    def _func(w: np.ndarray) -> tuple[float, np.ndarray]:
        W = _adj(w)
        loss, G_loss = _loss(W)
        h_current, G_h = _h(W)
        objective = loss + 0.5 * rho * h_current * h_current + alpha * h_current + lambda1 * w.sum()
        G_smooth = G_loss + (rho * h_current + alpha) * G_h
        gradient = np.concatenate(
            (G_smooth + lambda1, -G_smooth + lambda1),
            axis=None,
        )
        return float(objective), gradient

    for outer_iter in range(max_iter):
        w_new = None
        h_new = None
        solver_result = None

        while rho < rho_max:
            solver_result = sopt.minimize(
                _func,
                w_est,
                method="L-BFGS-B",
                jac=True,
                bounds=bounds,
            )
            w_new = solver_result.x
            h_new, _ = _h(_adj(w_new))

            if h_new > 0.25 * h_value:
                rho *= 10.0
            else:
                break

        if w_new is None or h_new is None or solver_result is None:
            raise RuntimeError("NOTEARS optimization did not produce a solution.")

        w_est = w_new
        h_value = h_new
        alpha += rho * h_value

        history.append(
            {
                "outer_iter": outer_iter + 1,
                "rho": float(rho),
                "alpha": float(alpha),
                "h": float(h_value),
                "solver_success": bool(solver_result.success),
                "solver_message": str(solver_result.message),
                "objective": float(solver_result.fun),
            }
        )

        if h_value <= h_tol or rho >= rho_max:
            break

    W_raw = _adj(w_est)
    np.fill_diagonal(W_raw, 0.0)
    if root_index is not None:
        W_raw[:, root_index] = 0.0

    final_h, _ = _h(W_raw)
    diagnostics = {
        "h_final": float(final_h),
        "rho_final": float(rho),
        "alpha_final": float(alpha),
        "outer_iterations": len(history),
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

    binary = (W_dag != 0).astype(int)
    pd.DataFrame(W_raw, index=columns, columns=columns).to_csv(
        output_dir / "notears_weighted_raw.csv",
        index_label="from",
    )
    pd.DataFrame(W_dag, index=columns, columns=columns).to_csv(
        output_dir / "notears_weighted_dag.csv",
        index_label="from",
    )
    pd.DataFrame(binary, index=columns, columns=columns).to_csv(
        output_dir / "notears_binary.csv",
        index_label="from",
    )

    with open(output_dir / "notears_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, ensure_ascii=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="CSV with observations in rows and nodes in columns.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--root", default="A", help="Known exogenous/root node; use 'none' to disable.")
    parser.add_argument("--lambda1", type=float, default=0.1)
    parser.add_argument("--loss-type", choices=("l2", "logistic", "poisson"), default="l2")
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--h-tol", type=float, default=1e-8)
    parser.add_argument("--rho-max", type=float, default=1e16)
    parser.add_argument("--threshold", type=float, default=0.3)
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
        X = data.to_numpy(dtype=float)
        scaling = None
    else:
        X, scaling = standardize_columns(data)

    W_raw, diagnostics = notears_linear_root(
        X=X,
        lambda1=args.lambda1,
        loss_type=args.loss_type,
        root_index=root_index,
        max_iter=args.max_iter,
        h_tol=args.h_tol,
        rho_max=args.rho_max,
    )
    W_dag, effective_threshold = postprocess_to_dag(W_raw, args.threshold)

    if root_index is not None and np.any(W_dag[:, root_index] != 0):
        raise RuntimeError("NOTEARS output violates the root constraint.")

    metadata = {
        "method": "NOTEARS-linear",
        "matrix_convention": "row i, column j denotes i -> j",
        "input": str(Path(args.input).resolve()),
        "n_samples": int(X.shape[0]),
        "n_nodes": int(X.shape[1]),
        "columns": columns,
        "root": root,
        "root_constraint": root is not None,
        "standardized": not args.no_standardize,
        "scaling": scaling,
        "lambda1": args.lambda1,
        "loss_type": args.loss_type,
        "max_iter": args.max_iter,
        "h_tol": args.h_tol,
        "rho_max": args.rho_max,
        "requested_threshold": args.threshold,
        "effective_dag_threshold": effective_threshold,
        "n_edges": int(np.count_nonzero(W_dag)),
        "is_dag": is_dag(W_dag),
        "diagnostics": diagnostics,
    }
    save_outputs(args.output_dir, columns, W_raw, W_dag, metadata)

    print(
        f"NOTEARS completed: nodes={X.shape[1]}, edges={np.count_nonzero(W_dag)}, "
        f"h={diagnostics['h_final']:.3e}, root={root}"
    )


if __name__ == "__main__":
    main()
