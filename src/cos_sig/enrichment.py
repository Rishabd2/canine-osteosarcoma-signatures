"""Pathway enrichment: GSEA (phenotype permutation) and hypergeometric ORA.

The paper used Broad GSEA 4.1.0 on MSigDB Hallmarks 7.5 for cluster contrasts
and Qiagen IPA for DEG lists. IPA is proprietary, so DEG lists here go through
a hypergeometric over-representation test on MSigDB Reactome 7.5.1 instead.
"""
from __future__ import annotations

from typing import Iterable, Mapping

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from .survival import bh_adjust


def gsea_phenotype(expr: pd.DataFrame, groups: pd.Series, gene_sets: Mapping[str, list[str]],
                   pos: str = "FP", neg: str = "PP", permutations: int = 1000, seed: int = 0,
                   threads: int = 4, outdir=None) -> pd.DataFrame:
    """Broad-style GSEA via gseapy: signal-to-noise ranking, phenotype permutation.

    ``expr`` is genes x samples. Returns Term, NES, FDR, NOM p, size; positive
    NES means enriched in ``pos``.
    """
    import gseapy as gp

    # gseapy takes the class of the first sample as the positive phenotype, so put ``pos`` first
    samples = [s for s in expr.columns if s in groups.index and groups[s] == pos]
    samples += [s for s in expr.columns if s in groups.index and groups[s] == neg]
    cls = [groups[s] for s in samples]
    data = expr[samples].copy()
    data.index.name = "Gene"
    res = gp.gsea(data=data.reset_index(), gene_sets=dict(gene_sets), cls=cls,
                  permutation_type="phenotype", permutation_num=permutations, method="signal_to_noise",
                  min_size=15, max_size=500, seed=seed, threads=threads, outdir=outdir, verbose=False)
    df = res.res2d.copy()
    out = pd.DataFrame({
        "term": df["Term"].astype(str),
        "nes": df["NES"].astype(float),
        "fdr": df["FDR q-val"].astype(float),
        "nom_p": df["NOM p-val"].astype(float),
    })
    return out.sort_values("nes", ascending=False).reset_index(drop=True)


def ora(genes: Iterable[str], universe: Iterable[str], gene_sets: Mapping[str, list[str]],
        min_size: int = 10, max_size: int = 500) -> pd.DataFrame:
    """Hypergeometric over-representation with BH FDR."""
    uni = set(universe)
    hits = set(genes) & uni
    n_uni, n_hit = len(uni), len(hits)
    rows = []
    for name, members in gene_sets.items():
        s = set(members) & uni
        if not (min_size <= len(s) <= max_size):
            continue
        k = len(s & hits)
        p = float(hypergeom.sf(k - 1, n_uni, len(s), n_hit)) if k > 0 else 1.0
        rows.append((name, len(s), k, p, ",".join(sorted(s & hits))))
    df = pd.DataFrame(rows, columns=["term", "set_size", "overlap", "p", "genes"])
    df["fdr"] = bh_adjust(df["p"].to_numpy()) if len(df) else []
    return df.sort_values("p").reset_index(drop=True)
