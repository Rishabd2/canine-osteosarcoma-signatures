"""Survival splits and statistics: log-rank, Cox HR, bootstrap CI, KM medians."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test


def median_split(score: pd.Series, high: str = "high", low: str = "low") -> pd.Series:
    """Split at the median: values strictly above the median are ``high``."""
    med = score.median()
    return score.map(lambda v: high if v > med else low).rename("group")


def align_survival(groups: pd.Series, time: pd.Series, event: pd.Series) -> pd.DataFrame:
    """Inner-join group labels with time/event on the index; drop missing rows."""
    df = pd.concat({"group": groups, "time": time, "event": event}, axis=1, join="inner").dropna()
    df["time"] = df["time"].astype(float)
    df["event"] = df["event"].astype(int)
    if not set(df["event"].unique()) <= {0, 1}:
        raise ValueError("event must be coded 0/1")
    return df


def logrank_p(df: pd.DataFrame) -> float:
    g = df["group"].unique()
    if len(g) < 2:
        return float("nan")
    if len(g) == 2:
        a, b = (df[df.group == x] for x in g)
        return float(logrank_test(a.time, b.time, a.event, b.event).p_value)
    return float(multivariate_logrank_test(df.time, df.group, df.event).p_value)


def cox_hr(df: pd.DataFrame, exposed: str, reference: str) -> tuple[float, float, float, float]:
    """HR of ``exposed`` vs ``reference`` (two groups). Returns (hr, lo95, hi95, wald p)."""
    d = df[df.group.isin([exposed, reference])].copy()
    d["x"] = (d.group == exposed).astype(float)
    cph = CoxPHFitter()
    cph.fit(d[["time", "event", "x"]], duration_col="time", event_col="event")
    s = cph.summary.loc["x"]
    return (float(s["exp(coef)"]), float(s["exp(coef) lower 95%"]), float(s["exp(coef) upper 95%"]), float(s["p"]))


def cox_beta_1d(time: np.ndarray, event: np.ndarray, x: np.ndarray, max_iter: int = 50,
                tol: float = 1e-9) -> float:
    """Cox log-HR for a single covariate, Efron ties, Newton-Raphson (matches lifelines)."""
    t = np.asarray(time, float)
    e = np.asarray(event, int)
    x = np.asarray(x, float)
    order = np.argsort(t, kind="mergesort")
    t, e, x = t[order], e[order], x[order]
    ev_times, d = np.unique(t[e == 1], return_counts=True)
    if ev_times.size == 0:
        return float("nan")
    start = np.searchsorted(t, ev_times, side="left")          # risk set = t >= u
    k_of = np.searchsorted(ev_times, t[e == 1])                 # event -> event-time index
    rep = np.repeat(np.arange(ev_times.size), d)                # one row per (time, l)
    frac = np.concatenate([np.arange(n) / n for n in d])
    xd_sum = np.bincount(k_of, weights=x[e == 1], minlength=ev_times.size)
    beta = 0.0
    for _ in range(max_iter):
        w = np.exp(beta * x)
        c0 = np.concatenate([np.cumsum(w[::-1])[::-1], [0.0]])
        c1 = np.concatenate([np.cumsum((w * x)[::-1])[::-1], [0.0]])
        s0, s1 = c0[start], c1[start]
        d0 = np.bincount(k_of, weights=w[e == 1], minlength=ev_times.size)
        d1 = np.bincount(k_of, weights=(w * x)[e == 1], minlength=ev_times.size)
        a0 = s0[rep] - frac * d0[rep]
        a1 = s1[rep] - frac * d1[rep]
        # x is binary or general: s2 needed for general x
        c2 = np.concatenate([np.cumsum((w * x * x)[::-1])[::-1], [0.0]])
        d2 = np.bincount(k_of, weights=(w * x * x)[e == 1], minlength=ev_times.size)
        a2 = c2[start][rep] - frac * d2[rep]
        grad = xd_sum.sum() - np.sum(a1 / a0)
        hess = -np.sum(a2 / a0 - (a1 / a0) ** 2)
        if hess == 0 or not np.isfinite(hess):
            return float("nan")
        step = float(np.clip(grad / hess, -1.0, 1.0))
        beta -= step
        if abs(step) < tol:
            break
    return float(beta)


def bootstrap_hr(df: pd.DataFrame, exposed: str, reference: str, n_boot: int = 1000,
                 seed: int = 0) -> tuple[float, float]:
    """Percentile 95% CI for the HR, resampling patients with labels held fixed.

    Replicates where one group has no events or vanishes are skipped.
    """
    rng = np.random.default_rng(seed)
    d = df[df.group.isin([exposed, reference])].reset_index(drop=True)
    x = (d.group == exposed).to_numpy(float)
    t = d.time.to_numpy(float)
    e = d.event.to_numpy(int)
    vals = []
    for _ in range(n_boot):
        i = rng.integers(0, len(d), len(d))
        if x[i].min() == x[i].max() or e[i][x[i] == 1].sum() == 0 or e[i][x[i] == 0].sum() == 0:
            continue
        b = cox_beta_1d(t[i], e[i], x[i])
        if np.isfinite(b) and abs(b) < 20:
            vals.append(float(np.exp(b)))
    if len(vals) < 0.5 * n_boot:
        return float("nan"), float("nan")
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(lo), float(hi)


def km_medians(df: pd.DataFrame) -> dict[str, float]:
    out = {}
    for g, sub in df.groupby("group"):
        kmf = KaplanMeierFitter().fit(sub.time, sub.event)
        m = kmf.median_survival_time_
        out[str(g)] = float(m) if np.isfinite(m) else float("nan")
    return out


@dataclass
class SplitStats:
    n: dict
    logrank_p: float
    hr: float
    hr_lo: float
    hr_hi: float
    cox_p: float
    boot_lo: float
    boot_hi: float
    medians: dict


def split_stats(groups: pd.Series, time: pd.Series, event: pd.Series, exposed: str, reference: str,
                n_boot: int = 1000, seed: int = 0) -> SplitStats:
    df = align_survival(groups, time, event)
    n = df.group.value_counts().to_dict()
    p = logrank_p(df)
    ev = df.groupby("group")["event"].sum()
    if df.group.nunique() == 2 and {exposed, reference} <= set(df.group) and ev.min() == 0:
        # a group with no events: the Cox HR is infinite or zero, report it as not estimable
        hr = lo = hi = cp = blo = bhi = float("nan")
    elif df.group.nunique() == 2 and {exposed, reference} <= set(df.group):
        hr, lo, hi, cp = cox_hr(df, exposed, reference)
        blo, bhi = bootstrap_hr(df, exposed, reference, n_boot=n_boot, seed=seed) if n_boot else (np.nan, np.nan)
    else:
        hr = lo = hi = cp = blo = bhi = float("nan")
    return SplitStats(n, p, hr, lo, hi, cp, blo, bhi, km_medians(df))


def bh_adjust(p) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN passed through)."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    q = p[ok]
    n = q.size
    if n == 0:
        return out
    order = np.argsort(q)
    ranked = q[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    res = np.empty(n)
    res[order] = np.minimum(ranked, 1.0)
    out[ok] = res
    return out
