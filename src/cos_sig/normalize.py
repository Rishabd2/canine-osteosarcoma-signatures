"""Expression filtering and normalization (genes x samples matrices)."""
from __future__ import annotations

import gzip
import re
from collections import defaultdict

import numpy as np
import pandas as pd


def filter_by_expr(counts: pd.DataFrame, group=None, min_count: float = 10, min_total_count: float = 15,
                   large_n: int = 10, min_prop: float = 0.7) -> pd.Series:
    """Python port of edgeR::filterByExpr (edgeR 3.32 to 4.x default branch).

    Keeps genes with CPM >= min_count / median(lib.size) * 1e6 in at least
    ``MinSampleSize`` samples and a total count >= min_total_count.
    """
    lib = counts.sum(axis=0).astype(float)
    if group is None:
        min_n = counts.shape[1]
    else:
        min_n = int(pd.Series(group).value_counts().min())
    if min_n > large_n:
        min_n = large_n + (min_n - large_n) * min_prop
    median_lib = float(np.median(lib))
    cpm_cutoff = min_count / median_lib * 1e6
    cpm = counts.div(lib, axis=1) * 1e6
    keep_cpm = (cpm >= cpm_cutoff).sum(axis=1) >= (min_n - 1e-14)
    keep_total = counts.sum(axis=1) >= (min_total_count - 1e-14)
    return keep_cpm & keep_total


def gene_lengths_from_gtf(path) -> pd.Series:
    """Union-of-exons length per ``gene_id`` from a (gzipped) GTF."""
    exons: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    rx = re.compile(r'gene_id "([^"]+)"')
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 8)
            if len(f) < 9 or f[2] != "exon":
                continue
            m = rx.search(f[8])
            if m:
                exons[m.group(1)].append((f[0], int(f[3]), int(f[4])))
    out = {}
    for gid, iv in exons.items():
        total = 0
        for chrom in {c for c, _, _ in iv}:
            segs = sorted((s, e) for c, s, e in iv if c == chrom)
            cs, ce = segs[0]
            for s, e in segs[1:]:
                if s > ce + 1:
                    total += ce - cs + 1
                    cs, ce = s, e
                else:
                    ce = max(ce, e)
            total += ce - cs + 1
        out[gid] = total
    return pd.Series(out, name="length")


def tpm(counts: pd.DataFrame, lengths: pd.Series) -> pd.DataFrame:
    """Transcripts per million from counts and gene lengths (rows must align)."""
    rate = counts.div(lengths.reindex(counts.index).astype(float), axis=0)
    return rate.div(rate.sum(axis=0), axis=1) * 1e6


def quantile_normalize(df: pd.DataFrame) -> pd.DataFrame:
    """Quantile normalization across columns (preprocessCore::normalize.quantiles, average ties)."""
    arr = df.to_numpy(dtype=float)
    sorted_mean = np.sort(arr, axis=0).mean(axis=1)
    out = np.empty_like(arr)
    n = arr.shape[0]
    grid = np.arange(n)
    for j in range(arr.shape[1]):
        ranks = pd.Series(arr[:, j]).rank(method="average").to_numpy() - 1
        out[:, j] = np.interp(ranks, grid, sorted_mean)
    return pd.DataFrame(out, index=df.index, columns=df.columns)
