"""k-means with silhouette-chosen k, as described in the paper's Methods."""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, silhouette_score


@dataclass
class KMeansResult:
    labels: pd.Series
    k: int
    silhouettes: dict[int, float] = field(default_factory=dict)


def kmeans_labels(x: pd.DataFrame, k: int, seed: int = 0, n_init: int = 10) -> pd.Series:
    km = KMeans(n_clusters=k, random_state=seed, n_init=n_init)
    return pd.Series(km.fit_predict(x.to_numpy()), index=x.index, name="cluster")


def kmeans_silhouette(x: pd.DataFrame, ks=(2, 3, 4, 5, 6), seed: int = 0, n_init: int = 10,
                      force_k: int | None = None) -> KMeansResult:
    """Fit k-means for each k, keep the k with the largest silhouette (or ``force_k``)."""
    sils = {}
    fits = {}
    for k in ks:
        lab = kmeans_labels(x, k, seed=seed, n_init=n_init)
        fits[k] = lab
        sils[k] = float(silhouette_score(x.to_numpy(), lab.to_numpy())) if lab.nunique() > 1 else float("nan")
    k_best = force_k if force_k is not None else max(sils, key=lambda k: sils[k])
    if k_best not in fits:
        fits[k_best] = kmeans_labels(x, k_best, seed=seed, n_init=n_init)
    return KMeansResult(fits[k_best], k_best, sils)


def orient_two_clusters(labels: pd.Series, score: pd.Series, high: str = "FP", low: str = "PP") -> pd.Series:
    """Name two clusters by their mean signature score, without using outcome.

    The paper calls the GS-1/GS-2 high-enrichment cluster "favorable prognosis"
    (FP). We fix that rule before looking at survival: the cluster with the
    higher mean score is FP.
    """
    if labels.nunique() != 2:
        raise ValueError("orient_two_clusters needs exactly two clusters")
    means = score.groupby(labels).mean()
    hi = means.idxmax()
    return labels.map(lambda c: high if c == hi else low).rename("group")


def seed_stability(x: pd.DataFrame, k: int, seeds=range(50), n_init: int = 1) -> dict[str, float]:
    """Adjusted Rand index of single-init k-means runs against the seed-0 run."""
    ref = kmeans_labels(x, k, seed=list(seeds)[0], n_init=n_init)
    aris = [adjusted_rand_score(ref, kmeans_labels(x, k, seed=s, n_init=n_init)) for s in list(seeds)[1:]]
    return {"ari_mean": float(np.mean(aris)), "ari_min": float(np.min(aris))}
