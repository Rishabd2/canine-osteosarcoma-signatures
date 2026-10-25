import numpy as np
import pandas as pd
import pytest
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

from cos_sig.clustering import orient_two_clusters
from cos_sig.survival import (align_survival, bh_adjust, bootstrap_hr, cox_beta_1d, logrank_p, median_split,
                              split_stats)


def _toy(seed=0, n=120):
    rng = np.random.default_rng(seed)
    x = rng.integers(0, 2, n)
    t = np.ceil(rng.exponential(1 / (0.4 + 0.6 * x)) * 100)
    e = (rng.random(n) < 0.8).astype(int)
    idx = [f"p{i}" for i in range(n)]
    g = pd.Series(np.where(x == 1, "PP", "FP"), index=idx)
    return g, pd.Series(t, index=idx), pd.Series(e, index=idx)


def test_median_split_is_balanced_and_strict():
    s = pd.Series([1, 2, 3, 4, 5, 6], index=list("abcdef"))
    g = median_split(s)
    assert g.value_counts().to_dict() == {"low": 3, "high": 3}
    assert g["c"] == "low" and g["d"] == "high"


def test_align_drops_missing_and_checks_event_coding():
    g = pd.Series(["FP", "PP", "PP"], index=["a", "b", "c"])
    t = pd.Series([10.0, np.nan, 5.0], index=["a", "b", "c"])
    e = pd.Series([1, 0, 1], index=["a", "b", "c"])
    df = align_survival(g, t, e)
    assert list(df.index) == ["a", "c"]
    with pytest.raises(ValueError):
        align_survival(g, t.fillna(1), pd.Series([1, 2, 0], index=["a", "b", "c"]))


def test_logrank_and_cox_match_lifelines():
    g, t, e = _toy()
    df = align_survival(g, t, e)
    a, b = df[df.group == "PP"], df[df.group == "FP"]
    assert logrank_p(df) == pytest.approx(logrank_test(a.time, b.time, a.event, b.event).p_value)
    x = (df.group == "PP").astype(float)
    ref = CoxPHFitter().fit(pd.DataFrame({"t": df.time, "e": df.event, "x": x}), "t", "e").params_["x"]
    assert cox_beta_1d(df.time.to_numpy(), df.event.to_numpy(), x.to_numpy()) == pytest.approx(ref, abs=1e-5)


def test_split_stats_direction_and_bootstrap():
    g, t, e = _toy(seed=3)
    s = split_stats(g, t, e, "PP", "FP", n_boot=300, seed=0)
    assert s.hr > 1  # PP simulated with the higher hazard
    assert s.boot_lo < s.hr < s.boot_hi
    assert s.n["PP"] + s.n["FP"] == len(g)


def test_zero_event_group_gives_nan_hr():
    g = pd.Series(["PP"] * 6 + ["FP"] * 6, index=range(12))
    t = pd.Series(np.arange(1, 13) * 10.0)  # PP die early, FP censored late
    e = pd.Series([1] * 6 + [0] * 6)
    s = split_stats(g, t, e, "PP", "FP", n_boot=50)
    assert np.isnan(s.hr) and s.logrank_p < 0.05


def test_bootstrap_is_seeded():
    g, t, e = _toy(seed=5)
    df = align_survival(g, t, e)
    assert bootstrap_hr(df, "PP", "FP", n_boot=200, seed=1) == bootstrap_hr(df, "PP", "FP", n_boot=200, seed=1)


def test_orient_two_clusters_uses_score_not_outcome():
    labels = pd.Series([0, 0, 1, 1], index=list("abcd"))
    score = pd.Series([2.0, 1.5, -1.0, -2.0], index=list("abcd"))
    out = orient_two_clusters(labels, score)
    assert out.to_dict() == {"a": "FP", "b": "FP", "c": "PP", "d": "PP"}


def test_bh_adjust_known_values():
    p = np.array([0.01, 0.04, 0.03, np.nan, 0.5])
    q = bh_adjust(p)
    np.testing.assert_allclose(q[[0, 1, 2, 4]], [0.04, 0.05333333, 0.05333333, 0.5])
    assert np.isnan(q[3])


def test_fast_logrank_matches_lifelines():
    from cos_sig.survival import logrank_2group_fast
    for seed in range(3):
        g, t, e = _toy(seed=seed)
        a, b = g == "PP", g == "FP"
        ref = logrank_test(t[a], t[b], e[a], e[b]).p_value
        assert logrank_2group_fast(t.to_numpy(), e.to_numpy(), a.to_numpy().astype(int)) == pytest.approx(ref, rel=1e-9)


def test_bootstrap_counts_zero_event_resamples_instead_of_dropping():
    from cos_sig.survival import bootstrap_hr_detail
    # 4 FP patients, only one FP event: many resamples leave FP with no events (HR = inf)
    g = pd.Series(["PP"] * 16 + ["FP"] * 4)
    t = pd.Series(list(range(10, 170, 10)) + [500, 600, 700, 800], dtype=float)
    e = pd.Series([1] * 16 + [1, 0, 0, 0])
    df = align_survival(g, t, e)
    r = bootstrap_hr_detail(df, "PP", "FP", n_boot=500, seed=0)
    assert r["n_inf"] > 0.05 * 500
    # the upper bound is reported as not estimable (inf, or NaN when too many resamples
    # are degenerate), never as a finite number computed from the surviving resamples
    assert not np.isfinite(r["hi"])
    assert r["n_inf"] + r["n_zero"] + r["n_degenerate"] <= 500
