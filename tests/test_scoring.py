import numpy as np
import pandas as pd
import pytest

from cos_sig.normalize import filter_by_expr, quantile_normalize, tpm
from cos_sig.scoring import log_zscore, signature_genes_present, signature_score, ssgsea


def _ssgsea_bruteforce(expr, members, alpha):
    n = expr.shape[0]
    out = []
    for j in expr.columns:
        order = expr[j].sort_values(ascending=False, kind="mergesort").index
        rankval = pd.Series(np.arange(n, 0, -1), index=order).astype(float)
        inset = np.array([g in members for g in order])
        w = rankval.to_numpy() ** alpha * inset
        p_hit = np.cumsum(w) / w.sum()
        p_miss = np.cumsum(~inset) / (n - inset.sum())
        out.append(np.sum(p_hit - p_miss))
    return np.array(out)


def test_log_zscore_matches_standard_scaler():
    rng = np.random.default_rng(0)
    expr = pd.DataFrame(rng.gamma(2, 10, size=(5, 8)), index=list("abcde"), columns=[f"s{i}" for i in range(8)])
    z = log_zscore(expr)
    assert z.shape == (8, 5)
    np.testing.assert_allclose(z.mean(axis=0), 0, atol=1e-12)
    np.testing.assert_allclose(z.std(axis=0, ddof=0), 1, atol=1e-12)
    x = np.log2(expr.loc["a"] + 2)
    np.testing.assert_allclose(z["a"], (x - x.mean()) / x.std(ddof=0))


def test_signature_score_uses_present_genes_only():
    z = pd.DataFrame({"g1": [1.0, -1.0], "g2": [3.0, -3.0], "g3": [0.0, 0.0]}, index=["s1", "s2"])
    present, missing = signature_genes_present(["g1", "g2", "gX"], z.columns)
    assert present == ["g1", "g2"] and missing == ["gX"]
    np.testing.assert_allclose(signature_score(z, ["g1", "g2", "gX"]), [2.0, -2.0])


def test_ssgsea_matches_bruteforce():
    rng = np.random.default_rng(1)
    genes = [f"g{i}" for i in range(60)]
    expr = pd.DataFrame(rng.normal(size=(60, 5)), index=genes, columns=list("ABCDE"))
    members = set(genes[:12])
    for alpha in (0.25, 0.75):
        got = ssgsea(expr, {"S": members}, alpha=alpha, normalize=False).loc["S"].to_numpy()
        np.testing.assert_allclose(got, _ssgsea_bruteforce(expr, members, alpha), rtol=1e-10)


def test_ssgsea_ranks_enriched_sample_highest():
    genes = [f"g{i}" for i in range(50)]
    expr = pd.DataFrame(0.0, index=genes, columns=["hot", "cold"])
    expr.loc[genes[:10], "hot"] = 5.0
    expr.loc[genes[:10], "cold"] = -5.0
    s = ssgsea(expr, {"S": genes[:10]}, normalize=False).loc["S"]
    assert s["hot"] > 0 > s["cold"]


def test_tpm_and_quantile_normalize():
    counts = pd.DataFrame({"a": [10, 20, 30], "b": [1, 2, 3]}, index=["x", "y", "z"])
    t = tpm(counts, pd.Series({"x": 1000, "y": 2000, "z": 1000}))
    np.testing.assert_allclose(t.sum(), 1e6)
    np.testing.assert_allclose(t["a"], t["b"])
    q = quantile_normalize(pd.DataFrame({"a": [5.0, 2.0, 3.0], "b": [4.0, 1.0, 6.0]}))
    assert sorted(q["a"]) == sorted(q["b"])


def test_filter_by_expr_threshold():
    counts = pd.DataFrame(np.full((3, 12), 100), index=["keep", "low", "zero"])
    counts.loc["low"] = [0] * 11 + [500]
    counts.loc["zero"] = 0
    keep = filter_by_expr(counts, min_count=10, min_prop=0.5, large_n=5)
    assert keep.tolist() == [True, False, False]


@pytest.mark.parametrize("bad", [[], ["nope"]])
def test_signature_matrix_errors_when_no_genes(bad):
    from cos_sig.scoring import signature_matrix
    with pytest.raises(ValueError):
        signature_matrix(pd.DataFrame({"g1": [0.0]}), bad)


def test_gsea_sign_follows_pos_even_when_first_sample_is_neg():
    from cos_sig.enrichment import gsea_phenotype
    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(300)]
    x = pd.DataFrame(rng.normal(size=(300, 20)), index=genes, columns=[f"s{i}" for i in range(20)])
    grp = pd.Series(["PP", "FP"] * 10, index=x.columns)  # first sample is PP
    x.loc[genes[:30], grp.index[grp == "FP"]] += 2
    res = gsea_phenotype(x, grp, {"UP_IN_FP": genes[:30]}, pos="FP", neg="PP", permutations=50, threads=1)
    assert res.set_index("term").loc["UP_IN_FP", "nes"] > 0
