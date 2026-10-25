"""Open immune deconvolution used in place of CIBERSORTx (which needs registration).

* ``mcp_counter``: MCP-counter (Becht et al. 2016) scores, the mean log2
  expression of each population's marker genes. Scores compare samples within
  a dataset; they are not proportions.
* ``til10_nnls``: constrained least squares of TPM on the quanTIseq TIL10
  signature matrix (Finotello et al. 2019): coefficients >= 0 and summing to at
  most 1, remainder reported as "Other". This follows quanTIseq's core
  regression; it omits quanTIseq's mRNA-content rescaling and its noisy-gene
  removal step, so values are relative fractions, not quanTIseq output.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize, nnls


def load_mcp_markers(path) -> dict[str, list[str]]:
    df = pd.read_csv(path, sep="\t")
    sym_col = "HUGO symbols" if "HUGO symbols" in df.columns else df.columns[0]
    out: dict[str, list[str]] = {}
    for _, r in df.iterrows():
        out.setdefault(r["Cell population"], []).append(r[sym_col])
    return out


def mcp_counter(log_expr: pd.DataFrame, markers: dict[str, list[str]], min_genes: int = 1) -> pd.DataFrame:
    """Populations x samples; ``log_expr`` is genes x samples on log2 scale."""
    rows = {}
    for pop, genes in markers.items():
        g = [x for x in genes if x in log_expr.index]
        if len(g) >= min_genes:
            rows[pop] = log_expr.loc[g].mean(axis=0)
    return pd.DataFrame(rows).T


def load_til10(path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", index_col=0)


def _fit_one(a: np.ndarray, y: np.ndarray) -> np.ndarray:
    b, _ = nnls(a, y)
    if b.sum() <= 1:
        return b
    k = a.shape[1]
    res = minimize(lambda v: np.sum((a @ v - y) ** 2), np.full(k, 1.0 / k),
                   jac=lambda v: 2 * a.T @ (a @ v - y), method="SLSQP",
                   bounds=[(0, None)] * k, constraints=[{"type": "ineq", "fun": lambda v: 1 - v.sum()}],
                   options={"maxiter": 500, "ftol": 1e-12})
    return np.clip(res.x, 0, None)


def til10_nnls(tpm: pd.DataFrame, signature: pd.DataFrame) -> pd.DataFrame:
    """Cell types (+Other) x samples. ``tpm`` is genes x samples, linear scale."""
    genes = [g for g in signature.index if g in tpm.index]
    a = signature.loc[genes].to_numpy(float)
    y_all = tpm.loc[genes].to_numpy(float)
    # mixture and signature are both TPM, so y_g ~ sum_i f_i * A_gi with sum f_i <= 1
    out = [_fit_one(a, y_all[:, j]) for j in range(y_all.shape[1])]
    res = pd.DataFrame(np.array(out).T, index=signature.columns, columns=tpm.columns)
    res.loc["Other"] = 1 - res.sum(axis=0)
    return res
