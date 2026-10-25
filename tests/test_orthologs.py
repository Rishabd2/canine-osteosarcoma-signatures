import numpy as np
import pandas as pd

from cos_sig.orthologs import build_dog_to_human_map, collapse_rows, map_signature, strip_suffix


def _orth():
    return pd.DataFrame({
        "dog_gene_id": ["d1", "d2", "d3a", "d3b", "d4", "d4"],
        "dog_symbol": ["CD8A", "MYH2", "GBP1", "GBP1", "SLC8A3", "SLC8A3"],
        "human_gene_id": ["h1", "h2", "h3", "h4", "h5", "h6"],
        "human_symbol": ["CD8A", "MYH2", "GBP3", "GBP1", "SLC8A3", "SLC8A2"],
        "homology_type": ["ortholog_one2one", "ortholog_one2one", "ortholog_one2many", "ortholog_one2many",
                          "ortholog_one2many", "ortholog_one2many"],
        "confidence": ["1", "1", "0", "1", "1", "1"],
    })


def test_strip_suffix():
    assert strip_suffix("A1BG_3") == "A1BG"
    assert strip_suffix("HSF4_12") == "HSF4"
    assert strip_suffix("DLA-DQA1") == "DLA-DQA1"


def test_one2one_and_one2many_rules():
    m = build_dog_to_human_map(["CD8A", "GBP1", "SLC8A3"], _orth(), ["CD8A", "GBP1"]).set_index("dog_symbol")
    assert m.loc["CD8A", "human_symbol"] == "CD8A"
    assert m.loc["CD8A", "method"] == "ensembl_one2one"
    # one-to-many: the same-name human ortholog wins
    assert m.loc["GBP1", "human_symbol"] == "GBP1"
    assert m.loc["GBP1", "method"] == "ensembl_best"
    assert m.loc["SLC8A3", "human_symbol"] == "SLC8A3"


def test_identity_and_dla_and_unmapped():
    m = build_dog_to_human_map(["DLA-DQA1", "IRX1", "LOC100855540"], _orth(),
                               ["HLA-DQA1", "IRX1"]).set_index("dog_symbol")
    assert m.loc["DLA-DQA1", "human_symbol"] == "HLA-DQA1"
    assert m.loc["IRX1", "method"] == "identity"
    assert m.loc["LOC100855540", "method"] == "unmapped"
    assert np.isnan(m.loc["LOC100855540", "human_symbol"])


def test_map_signature_and_collapse():
    mp = build_dog_to_human_map(["CD8A", "MYH2", "LOC1"], _orth(), ["CD8A", "MYH2"])
    assert map_signature(["CD8A", "LOC1", "NOTHERE"], mp) == {"CD8A": "CD8A", "LOC1": None, "NOTHERE": None}
    expr = pd.DataFrame({"s1": [1.0, 3.0, 5.0], "s2": [2.0, 4.0, 6.0]}, index=["a_1", "a_2", "b_1"])
    out = collapse_rows(expr, pd.Series(["A", "A", np.nan], index=expr.index))
    assert list(out.index) == ["A"]
    assert out.loc["A", "s1"] == 2.0 and out.loc["A", "s2"] == 3.0
