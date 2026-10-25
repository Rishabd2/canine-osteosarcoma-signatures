"""Shared analysis steps: cluster samples on a signature, then test survival."""
from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from .clustering import kmeans_labels, kmeans_silhouette, orient_two_clusters
from .scoring import log_zscore, signature_genes_present, signature_score
from .survival import split_stats

warnings.filterwarnings("ignore", module="lifelines")

ROW_COLS = ["key", "claim", "cohort", "signature", "endpoint", "method", "group_a", "group_b", "n_a", "n_b",
            "k_silhouette", "genes_used", "logrank_p", "hr_a_vs_b", "hr_lo", "hr_hi", "boot_lo", "boot_hi",
            "median_a", "median_b", "ari_vs_paper"]


def cluster_on_signature(expr_lin: pd.DataFrame, genes: list[str], seed: int = 0,
                         hi: str = "FP", lo: str = "PP") -> dict:
    """k-means (k from silhouette over 2..6) on z-scored log2(x+2) signature genes.

    Returns the silhouette choice and the two-cluster solution named by mean
    signature score (``hi`` = higher score), plus genes used.
    """
    z = log_zscore(expr_lin)
    present, missing = signature_genes_present(genes, z.columns)
    x = z[present]
    sil = kmeans_silhouette(x, seed=seed)
    two = orient_two_clusters(kmeans_labels(x, 2, seed=seed), signature_score(z, present), high=hi, low=lo)
    return {"k_silhouette": sil.k, "silhouettes": sil.silhouettes, "labels": two, "present": present,
            "missing": missing, "z": x}


def survival_row(key: str, claim: int, cohort: str, signature: str, endpoint: str, method: str,
                 groups: pd.Series, time: pd.Series, event: pd.Series, a: str, b: str,
                 k_sil=np.nan, genes_used=np.nan, paper_labels: pd.Series | None = None,
                 n_boot: int = 1000, seed: int = 0) -> dict:
    s = split_stats(groups, time, event, a, b, n_boot=n_boot, seed=seed)
    ari = np.nan
    if paper_labels is not None:
        common = groups.index.intersection(paper_labels.dropna().index)
        ari = adjusted_rand_score(paper_labels[common].astype(str), groups[common].astype(str))
    return {"key": key, "claim": claim, "cohort": cohort, "signature": signature, "endpoint": endpoint,
            "method": method, "group_a": a, "group_b": b, "n_a": int(s.n.get(a, 0)), "n_b": int(s.n.get(b, 0)),
            "k_silhouette": k_sil, "genes_used": genes_used, "logrank_p": s.logrank_p, "hr_a_vs_b": s.hr,
            "hr_lo": s.hr_lo, "hr_hi": s.hr_hi, "boot_lo": s.boot_lo, "boot_hi": s.boot_hi,
            "median_a": s.medians.get(a, np.nan), "median_b": s.medians.get(b, np.nan), "ari_vs_paper": ari}
