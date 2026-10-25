"""Checks the paper did not run.

1. Choice of k: for k = 2..6, silhouette and multi-group log-rank p.
2. k-means seed stability: adjusted Rand index of 50 single-start runs.
3. Whole-pipeline bootstrap: resample patients, re-cluster (k = 2), re-test; the
   share of resamples with log-rank p < 0.05 says how often the finding recurs.
4. Random-signature null in TARGET: 1000 random gene sets of the same size,
   clustered the same way; empirical p = share of random sets whose log-rank p is
   at or below the observed one.
5. Leakage: the frozen de novo signatures verify against their sha256, and every
   human split uses orientation by signature score only.
Writes results/sensitivity_*.csv.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.clustering import kmeans_labels, seed_stability  # noqa: E402
from cos_sig.orthologs import map_signature  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.scoring import log_zscore, signature_genes_present  # noqa: E402
from cos_sig.signatures import PUBLISHED_ISA, load_frozen, published_gs  # noqa: E402
from cos_sig.survival import align_survival, logrank_2group_fast, logrank_p  # noqa: E402
from sklearn.metrics import silhouette_score  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
SEED = 0
N_PIPE_BOOT = 200
N_RANDOM = 1000


def lr(labels, time, event) -> float:
    return logrank_p(align_survival(labels, time, event))


def lr2(labels: np.ndarray, time: np.ndarray, event: np.ndarray) -> float:
    """Fast two-group log-rank for resampling loops (labels 0/1, arrays aligned)."""
    return logrank_2group_fast(time, event, labels)


def main() -> None:
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    gs = published_gs(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx")
    frozen = load_frozen(RESULTS / "frozen_signatures.json")
    dt = pd.read_pickle(W / "dog_tpm_qn.pkl")
    dt.columns = dt.columns.astype(str)
    dm = pd.read_pickle(W / "dog_meta.pkl")
    dm.index = dm.index.astype(str)
    ht = pd.read_pickle(W / "hum_tpm_qn.pkl")
    hm = pd.read_pickle(W / "hum_meta.pkl")
    hg = lambda g: [h for h in map_signature(g, omap).values() if h]  # noqa: E731
    nonmet = hm.index[hm.metastatic == False]  # noqa: E712
    met = hm.index[hm.metastatic == True]  # noqa: E712
    cases = [
        ("DOG2", "GS-1", dt, gs["GS-1"], dm, ["dfi", "os"]),
        ("DOG2", "GS-2", dt, gs["GS-2"], dm, ["dfi", "os"]),
        ("TARGET all", "GS-1", ht, hg(gs["GS-1-human"]), hm, ["pfs", "os"]),
        ("TARGET all", "GS-2", ht, hg(gs["GS-2"]), hm, ["pfs", "os"]),
        ("TARGET metastatic", "GS-1", ht[met], hg(gs["GS-1-human"]), hm.loc[met], ["pfs", "os"]),
        ("TARGET metastatic", "GS-2", ht[met], hg(gs["GS-2"]), hm.loc[met], ["pfs", "os"]),
        ("TARGET all", "ISA-3", ht, hg(PUBLISHED_ISA["ISA-3"]), hm, ["pfs"]),
        ("TARGET non-metastatic", "ISA-1", ht[nonmet], hg(PUBLISHED_ISA["ISA-1"]), hm.loc[nonmet], ["pfs"]),
        ("TARGET non-metastatic", "ISA-3", ht[nonmet], hg(PUBLISHED_ISA["ISA-3"]), hm.loc[nonmet], ["pfs"]),
        ("TARGET all", "ISA-4", ht, hg(PUBLISHED_ISA["ISA-4"]), hm, ["pfs"]),
    ]
    n_primary = len(cases)  # resampling checks run on the paper's signatures only
    for name, genes in frozen.items():
        cases.append(("TARGET all", name, ht, hg(genes), hm, ["pfs", "os"]))
        cases.append(("TARGET non-metastatic", name, ht[nonmet], hg(genes), hm.loc[nonmet], ["pfs"]))

    k_rows, stab_rows, boot_rows, null_rows = [], [], [], []
    rng = np.random.default_rng(SEED)
    for ci, (cohort, sname, expr, genes, meta, eps) in enumerate(cases):
        z = log_zscore(expr)
        present, _ = signature_genes_present(genes, z.columns)
        x = z[present]
        # 1. k sensitivity
        for k in range(2, 7):
            lab = kmeans_labels(x, k, seed=SEED)
            row = {"cohort": cohort, "signature": sname, "k": k,
                   "silhouette": silhouette_score(x.to_numpy(), lab.to_numpy()),
                   "min_cluster": int(lab.value_counts().min())}
            for ep in eps:
                row[f"logrank_p_{ep}"] = lr(lab.astype(str), meta[ep], meta[f"{ep}_event"])
            k_rows.append(row)
        # 2. seed stability
        stab_rows.append({"cohort": cohort, "signature": sname, "genes": len(present), **seed_stability(x, 2)})
        if ci >= n_primary:
            print("done", cohort, sname, flush=True)
            continue
        # 3. whole-pipeline bootstrap (re-z-score and re-cluster each resample)
        ids = np.asarray(expr.columns)
        hits = {ep: [] for ep in eps}
        for b in range(N_PIPE_BOOT):
            pick = rng.choice(ids, size=len(ids), replace=True)
            e = expr[pick]
            e.columns = [f"{c}#{i}" for i, c in enumerate(pick)]
            zb = log_zscore(e)[present]
            lab = kmeans_labels(zb, 2, seed=b, n_init=5).to_numpy()
            m = meta.loc[pick]
            for ep in eps:
                hits[ep].append(lr2(lab, m[ep].to_numpy(float), m[f"{ep}_event"].to_numpy(int)))
        for ep in eps:
            p = np.array(hits[ep])
            boot_rows.append({"cohort": cohort, "signature": sname, "endpoint": ep.upper(),
                              "share_p_lt_0.05": float(np.nanmean(p < 0.05)),
                              "median_p": float(np.nanmedian(p)), "n_resamples": N_PIPE_BOOT})
        # 4. random-signature null (TARGET only)
        if cohort.startswith("TARGET"):
            pool = z.columns[(expr.mean(axis=1) > 1).reindex(z.columns).fillna(False).to_numpy()]
            obs_lab = kmeans_labels(x, 2, seed=SEED)
            for ep in eps:
                obs = lr(obs_lab.astype(str), meta[ep], meta[f"{ep}_event"])
                null = []
                tt = meta.loc[z.index, ep].to_numpy(float)
                ee = meta.loc[z.index, f"{ep}_event"].to_numpy(int)
                for r in range(N_RANDOM):
                    g = rng.choice(pool, size=len(present), replace=False)
                    null.append(lr2(kmeans_labels(z[g], 2, seed=r, n_init=3).to_numpy(), tt, ee))
                null = np.array(null)
                null_rows.append({"cohort": cohort, "signature": sname, "endpoint": ep.upper(),
                                  "genes": len(present), "observed_p": obs,
                                  "empirical_p": float((np.sum(null <= obs) + 1) / (len(null) + 1)),
                                  "share_random_p_lt_0.05": float(np.mean(null < 0.05)),
                                  "n_random": N_RANDOM})
        print("done", cohort, sname, flush=True)
    pd.DataFrame(k_rows).to_csv(RESULTS / "sensitivity_k.csv", index=False)
    pd.DataFrame(stab_rows).to_csv(RESULTS / "sensitivity_seed_stability.csv", index=False)
    pd.DataFrame(boot_rows).to_csv(RESULTS / "sensitivity_pipeline_bootstrap.csv", index=False)
    pd.DataFrame(null_rows).to_csv(RESULTS / "sensitivity_random_signature_null.csv", index=False)
    leak = {"frozen_sha256_verified": True,
            "frozen_signatures": {k: len(v) for k, v in frozen.items()},
            "orientation_rule": "FP/high = cluster with the higher mean signature z-score; outcome not used",
            # static check: the derivation script never names a human (hum_*) work file
            "human_data_read_by_derivation_script": "hum_" in (ROOT / "scripts" / "03_denovo_dog.py").read_text()}
    if leak["human_data_read_by_derivation_script"]:
        raise RuntimeError("scripts/03_denovo_dog.py references human data")
    (RESULTS / "leakage_checks.json").write_text(json.dumps(leak, indent=1))
    print(pd.DataFrame(boot_rows).to_string())
    print(pd.DataFrame(null_rows).to_string())


if __name__ == "__main__":
    main()
