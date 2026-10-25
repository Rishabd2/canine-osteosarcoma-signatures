"""Signature scoring: per-gene z-scores, mean signature score, and ssGSEA."""
from __future__ import annotations

from typing import Iterable, Mapping

import numpy as np
import pandas as pd


def log_zscore(expr: pd.DataFrame, pseudo: float = 2.0, log: bool = True) -> pd.DataFrame:
    """Samples x genes z-scores, as in the authors' km_clust.py.

    ``expr`` is genes x samples on a linear scale (TPM or CPM). The authors
    apply log2(x + 2) and then sklearn StandardScaler (population sd, ddof=0)
    per gene. Genes with zero variance get z = 0.
    """
    x = np.log2(expr + pseudo) if log else expr.copy()
    x = x.T
    mu = x.mean(axis=0)
    sd = x.std(axis=0, ddof=0).replace(0, np.nan)
    return ((x - mu) / sd).fillna(0.0)


def signature_genes_present(genes: Iterable[str], available: Iterable[str]) -> tuple[list[str], list[str]]:
    """Split a signature into (present, missing) relative to ``available`` gene names."""
    avail = set(available)
    genes = list(dict.fromkeys(genes))
    present = [g for g in genes if g in avail]
    missing = [g for g in genes if g not in avail]
    return present, missing


def signature_matrix(z: pd.DataFrame, genes: Iterable[str]) -> pd.DataFrame:
    """Columns of a samples x genes z matrix restricted to the signature (order preserved)."""
    present, _ = signature_genes_present(genes, z.columns)
    if not present:
        raise ValueError("no signature genes present in the matrix")
    return z[present]


def signature_score(z: pd.DataFrame, genes: Iterable[str]) -> pd.Series:
    """Mean z-score over the signature genes present (one value per sample)."""
    return signature_matrix(z, genes).mean(axis=1)


def ssgsea(expr: pd.DataFrame, gene_sets: Mapping[str, Iterable[str]], alpha: float = 0.25,
           min_size: int = 10, normalize: bool = True) -> pd.DataFrame:
    """Single-sample GSEA (Barbie et al. 2009; GSVA ``method='ssgsea'`` formulation).

    ``expr`` is genes x samples. For each sample genes are ranked by expression
    (highest first). The running sum steps up by |rank|^alpha / sum over the set
    at set members and down by 1 / (N - |S|) elsewhere; the score is the sum of
    the running-sum differences. Returns sets x samples. With ``normalize`` the
    scores are divided by (max - min) over the whole matrix, as GSVA does.
    """
    genes = np.asarray(expr.index)
    n = len(genes)
    pos = {g: i for i, g in enumerate(genes)}
    vals = expr.to_numpy(dtype=float)
    # rank 1 = lowest expression, n = highest (ties by order, as GSVA does)
    order = np.argsort(-vals, axis=0, kind="mergesort")
    rank_desc = np.empty_like(order)
    for j in range(vals.shape[1]):
        rank_desc[order[:, j], j] = np.arange(n)  # 0 = highest
    weight = (n - rank_desc).astype(float) ** alpha  # rank value n..1 for highest..lowest
    out = {}
    for name, members in gene_sets.items():
        idx = np.array(sorted({pos[g] for g in members if g in pos}), dtype=int)
        k = idx.size
        if k < min_size or k >= n:
            continue
        scores = np.empty(vals.shape[1])
        for j in range(vals.shape[1]):
            r = rank_desc[idx, j]
            w = weight[idx, j]
            o = np.argsort(r)
            r, w = r[o], w[o]
            hit_cum = np.cumsum(w) / w.sum()
            # running sums evaluated at each position; integrate exactly
            # P_hit(i) - P_miss(i) summed over i = 0..n-1
            # between consecutive hits P_hit is constant, P_miss grows linearly
            miss_step = 1.0 / (n - k)
            total = 0.0
            prev_hit_val = 0.0
            prev_pos = -1
            misses_so_far = 0
            for t in range(k):
                p = r[t]
                gap = p - prev_pos - 1  # number of non-members between hits
                if gap > 0:
                    m0 = misses_so_far + 1
                    m1 = misses_so_far + gap
                    total += gap * prev_hit_val - miss_step * (m0 + m1) * gap / 2.0
                    misses_so_far += gap
                total += hit_cum[t] - miss_step * misses_so_far
                prev_hit_val = hit_cum[t]
                prev_pos = p
            gap = n - 1 - prev_pos
            if gap > 0:
                m0 = misses_so_far + 1
                m1 = misses_so_far + gap
                total += gap * prev_hit_val - miss_step * (m0 + m1) * gap / 2.0
            scores[j] = total
        out[name] = scores
    res = pd.DataFrame(out, index=expr.columns).T
    if normalize and res.size:
        rng = res.to_numpy().max() - res.to_numpy().min()
        if rng > 0:
            res = res / rng
    return res


def read_gmt(path) -> dict[str, list[str]]:
    sets = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                sets[f[0]] = [g for g in f[2:] if g]
    return sets
