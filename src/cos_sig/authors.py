"""Run the authors' km_clust.py (Zenodo 10.5281/zenodo.8125077) on our inputs.

The script is executed from the downloaded copy, with four textual patches that
are applied and checked one by one:

1. ``np.argmax(sil_score)[0]`` raises ``IndexError`` on a NumPy scalar in any
   NumPy release, so the published script cannot run as written. We drop ``[0]``.
2. ``np.in1d`` was removed in NumPy 2.4; every call becomes ``np.isin``, which
   gives the same result for 1-D inputs.
3. ``KMeans(...)`` gets ``n_init=10``, the scikit-learn 0.23.1 default the
   authors ran with (current scikit-learn defaults to a single start).
4. The signature file name is pointed at the signature we pass in.

Nothing else changes. The script's own z-scoring, k range, and silhouette choice run
unchanged. It calls KMeans without a random_state, so we seed NumPy's global
generator before running it, which fixes the result.
"""
from __future__ import annotations

import os
import runpy
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

PATCHES = [
    ("idxk=np.argmax(sil_score)[0]", "idxk=np.argmax(sil_score)"),
    # scikit-learn 0.23.1 (the authors' version) defaulted to n_init=10; current releases
    # default to n_init="auto" (one k-means++ start). Restore the original default.
    ("kmeans=KMeans(n_clusters=k)", "kmeans=KMeans(n_clusters=k, n_init=10)"),
    ("km=KMeans(n_clusters=kmv[idxk])", "km=KMeans(n_clusters=kmv[idxk], n_init=10)"),
]


def patched_source(src: str, signature_file: str) -> str:
    for old, new in PATCHES:
        if src.count(old) != 1:
            raise RuntimeError(f"expected exactly one occurrence of {old!r} in km_clust.py")
        src = src.replace(old, new)
    if "np.in1d(" not in src:
        raise RuntimeError("expected np.in1d calls in km_clust.py")
    src = src.replace("np.in1d(", "np.isin(")
    old = "f3='Gardner_GS_1_genes.csv'"
    if src.count(old) != 1:
        raise RuntimeError("signature file line not found in km_clust.py")
    return src.replace(old, f"f3='{signature_file}'")


def run_km_clust(script: Path, expr_lin: pd.DataFrame, genes: list[str], time: pd.Series,
                 event: pd.Series, seed: int = 0) -> pd.Series:
    """``expr_lin``: samples x genes on a linear (TPM-like) scale. Returns cluster labels."""
    src = patched_source(Path(script).read_text(), "signature.csv")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        x = expr_lin.copy()
        x.columns = [f"{g}_x" for g in x.columns]  # the script strips a 2-character suffix
        x.index = x.index.astype(str)
        x.to_csv(td / "RNA_SEQ_FILE.csv")
        clin = pd.DataFrame({"DFI ": time.values, "DFI.Censor": event.values}, index=time.index.astype(str))
        clin.to_csv(td / "Clinical_data.csv")
        pd.DataFrame({"gene": list(genes)}).to_csv(td / "signature.csv")
        (td / "km_clust_patched.py").write_text(src)
        cwd = os.getcwd()
        try:
            os.chdir(td)
            np.random.seed(seed)
            runpy.run_path(str(td / "km_clust_patched.py"), run_name="__main__")
            out = pd.read_csv(td / "Gene_Signature_k_means_clustered_data_out.csv", index_col=0)
        finally:
            os.chdir(cwd)
    lab = out["Cluster"].astype(int)
    lab.index = lab.index.astype(str)
    return lab
