"""Canine to human gene mapping.

The canine counts in GSE238110 are keyed by NCBI CanFam4 gene symbols with a
numeric suffix (for example ``A1BG_3``). Human TARGET-OS data are keyed by
HGNC symbols. We map canine symbols to human symbols with Ensembl Compara
orthologs (BioMart, dog dataset ``clfamiliaris_gene_ensembl``), in this order:

1. ``ensembl_one2one``: the canine symbol has a one-to-one human ortholog.
2. ``ensembl_best``: one-to-many or many-to-many; we take the human ortholog
   whose symbol equals the canine symbol if present, else a high-confidence
   ortholog, else the alphabetically first one (deterministic).
3. ``identity``: no Ensembl ortholog row for that symbol, but the canine symbol
   (after ``DLA-`` to ``HLA-`` for the dog leukocyte antigen genes) is a human
   symbol in the target dataset.
4. ``unmapped``: everything else (mostly ``LOC`` identifiers).
"""
from __future__ import annotations

import re
from typing import Iterable

import numpy as np
import pandas as pd

_SUFFIX = re.compile(r"_\d+$")


def strip_suffix(symbol: str) -> str:
    """Remove the trailing ``_<n>`` that the GEO counts file appends to symbols."""
    return _SUFFIX.sub("", str(symbol))


def load_ensembl_orthologs(path) -> pd.DataFrame:
    """Read the BioMart TSV written by scripts/fetch_data.py into tidy columns."""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df.columns = ["dog_gene_id", "dog_symbol", "human_gene_id", "human_symbol", "homology_type", "confidence"]
    df = df.dropna(subset=["dog_symbol", "human_symbol"])
    df = df[(df.dog_symbol != "") & (df.human_symbol != "")]
    return df.reset_index(drop=True)


def _pick(rows: pd.DataFrame, dog_symbol: str) -> tuple[str, str]:
    one = rows[rows.homology_type == "ortholog_one2one"]
    if len(one["human_symbol"].unique()) == 1:
        return one["human_symbol"].iloc[0], "ensembl_one2one"
    cands = rows["human_symbol"].unique().tolist()
    if dog_symbol in cands:
        return dog_symbol, "ensembl_best"
    hi = sorted(rows.loc[rows.confidence == "1", "human_symbol"].unique().tolist())
    if hi:
        return hi[0], "ensembl_best"
    return sorted(cands)[0], "ensembl_best"


def build_dog_to_human_map(dog_symbols: Iterable[str], orthologs: pd.DataFrame,
                           human_symbols: Iterable[str]) -> pd.DataFrame:
    """Return one row per canine symbol: dog_symbol, human_symbol (NaN if unmapped), method."""
    human_set = set(human_symbols)
    by_dog = {k: v for k, v in orthologs.groupby("dog_symbol")}
    out = []
    for sym in pd.unique(pd.Series(list(dog_symbols), dtype=str)):
        if sym in by_dog:
            hs, method = _pick(by_dog[sym], sym)
        else:
            alt = "HLA-" + sym[4:] if sym.startswith("DLA-") else sym
            if alt in human_set:
                hs, method = alt, "identity"
            else:
                hs, method = np.nan, "unmapped"
        out.append((sym, hs, method))
    return pd.DataFrame(out, columns=["dog_symbol", "human_symbol", "method"])


def collapse_rows(expr: pd.DataFrame, keys: pd.Series, how: str = "mean") -> pd.DataFrame:
    """Aggregate rows of a genes x samples matrix that share a key (NaN keys dropped)."""
    keys = pd.Series(keys.values, index=expr.index)
    ok = keys.notna()
    grouped = expr[ok].groupby(keys[ok].values)
    return grouped.mean() if how == "mean" else grouped.sum()


def map_signature(genes: Iterable[str], mapping: pd.DataFrame) -> dict[str, str | None]:
    """Map a list of canine signature symbols to human symbols (None if unmapped)."""
    m = mapping.set_index("dog_symbol")["human_symbol"].to_dict()
    res = {}
    for g in genes:
        h = m.get(g)
        res[g] = None if (h is None or (isinstance(h, float) and np.isnan(h))) else h
    return res
