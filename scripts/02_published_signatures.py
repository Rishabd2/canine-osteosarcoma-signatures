"""Claims 1 and 3 (GS part): published GS-1 / GS-2 in DOG2 and TARGET-OS.

For each cohort and signature we report three versions of the survival split:
  paper_labels   the paper's own cluster labels (Supplementary Data 1) with the
                 paper's outcomes, re-tested here (checks the figures' arithmetic)
  rederived      our clustering from public expression data (pipeline in Methods)
  authors_code   the authors' km_clust.py (patched only to run) on the same input
Writes results/claims_published_signatures.csv and results/cluster_assignments_*.csv.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.authors import run_km_clust  # noqa: E402
from cos_sig.clustering import orient_two_clusters  # noqa: E402
from cos_sig.orthologs import map_signature  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.pipeline import ROW_COLS, cluster_on_signature, survival_row  # noqa: E402
from cos_sig.scoring import log_zscore, signature_score  # noqa: E402
from cos_sig.signatures import published_gs  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
SEED = 0
N_BOOT = 1000


def main() -> None:
    sigs = published_gs(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx")
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    rows, notes = [], {}

    # ---------------- DOG2 (claim 1) ----------------
    dt = pd.read_pickle(W / "dog_tpm_qn.pkl")
    dt.columns = dt.columns.astype(str)
    dm = pd.read_pickle(W / "dog_meta.pkl")
    dm.index = dm.index.astype(str)
    dog_assign = pd.DataFrame(index=dm.index)
    for sname, col in [("GS-1", "paper_gs1"), ("GS-2", "paper_gs2")]:
        genes = sigs[sname]
        res = cluster_on_signature(dt, genes, seed=SEED)
        lab = res["labels"]
        dog_assign[f"{sname}_rederived"] = lab
        dog_assign[f"{sname}_paper"] = dm[col]
        # authors' script on the same input (samples x genes, linear scale)
        al = run_km_clust(D / "zenodo" / "km_clust.py", dt.T, res["present"], dm["dfi"], dm["dfi_event"], seed=SEED)
        z = log_zscore(dt)
        al = orient_two_clusters(al.reindex(dm.index), signature_score(z, res["present"])) if al.nunique() == 2 else al
        dog_assign[f"{sname}_authors_code"] = al
        notes[f"dog_{sname}"] = {"genes_listed": len(genes), "genes_present": len(res["present"]),
                                 "missing": res["missing"], "silhouettes": res["silhouettes"],
                                 "authors_code_k": int(al.nunique()),
                                 "ari_rederived_vs_authors_code": adjusted_rand_score(lab, al.reindex(lab.index))}
        for ep, tcol in [("DFI", "dfi"), ("OS", "os")]:
            key = f"dog_{sname.lower().replace('-', '')}_{ep.lower()}"
            for method, groups in [("paper_labels", dm[col]), ("rederived", lab), ("authors_code", al)]:
                rows.append(survival_row(key, 1, "DOG2", sname, ep, method, groups, dm[tcol], dm[f"{tcol}_event"],
                                         "PP", "FP", k_sil=res["k_silhouette"] if method == "rederived" else np.nan,
                                         genes_used=len(res["present"]),
                                         paper_labels=None if method == "paper_labels" else dm[col],
                                         n_boot=N_BOOT, seed=SEED))
    dog_assign.index.name = "dog_pid"
    dog_assign.to_csv(RESULTS / "cluster_assignments_dog.csv")

    # ---------------- TARGET-OS (claim 3) ----------------
    ht = pd.read_pickle(W / "hum_tpm_qn.pkl")
    hm = pd.read_pickle(W / "hum_meta.pkl")
    hum_assign = pd.DataFrame(index=hm.index)
    cohorts = {"TARGET all": hm.index, "TARGET non-metastatic": hm.index[hm.metastatic == False],  # noqa: E712
               "TARGET metastatic": hm.index[hm.metastatic == True]}  # noqa: E712
    prefix = {"TARGET all": "hum", "TARGET non-metastatic": "hum_nonmet", "TARGET metastatic": "hum_met"}
    for sname, dog_key in [("GS-1", "GS-1-human"), ("GS-2", "GS-2")]:
        genes = [h for h in map_signature(sigs[dog_key], omap).values() if h]
        for cname, idx in cohorts.items():
            res = cluster_on_signature(ht[idx], genes, seed=SEED)
            lab = res["labels"]
            tag = f"{sname}_{prefix[cname]}"
            hum_assign[f"{tag}_rederived"] = lab
            if cname == "TARGET all":
                paper = hm[f"paper_{sname.lower().replace('-', '')}"]
            elif cname == "TARGET metastatic":
                paper = hm.loc[idx, f"paper_{sname.lower().replace('-', '')}_met"]
            else:
                paper = None  # the paper's non-metastatic labels are not in the source data
            if paper is not None:
                hum_assign[f"{tag}_paper"] = paper
            al = None
            if cname == "TARGET all":
                al = run_km_clust(D / "zenodo" / "km_clust.py", ht[idx].T, res["present"], hm.loc[idx, "pfs"],
                                  hm.loc[idx, "pfs_event"], seed=SEED)
                z = log_zscore(ht[idx])
                al = orient_two_clusters(al.reindex(idx), signature_score(z, res["present"])) if al.nunique() == 2 else al
                hum_assign[f"{tag}_authors_code"] = al
            notes[f"{tag}"] = {"genes_mapped": len(genes), "genes_present": len(res["present"]),
                               "missing": res["missing"], "silhouettes": res["silhouettes"]}
            for ep, tcol in [("PFS", "pfs"), ("OS", "os")]:
                key = f"{prefix[cname]}_{sname.lower().replace('-', '')}_{ep.lower()}"
                variants = [("rederived", lab)]
                if paper is not None:
                    variants.insert(0, ("paper_labels", paper))
                if al is not None:
                    variants.append(("authors_code", al))
                for method, groups in variants:
                    rows.append(survival_row(key, 3, cname, sname, ep, method, groups, hm.loc[idx, tcol],
                                             hm.loc[idx, f"{tcol}_event"], "PP", "FP",
                                             k_sil=res["k_silhouette"] if method == "rederived" else np.nan,
                                             genes_used=len(res["present"]),
                                             paper_labels=None if method == "paper_labels" else paper,
                                             n_boot=N_BOOT, seed=SEED))
    hum_assign.index.name = "target_case"
    hum_assign.to_csv(RESULTS / "cluster_assignments_target.csv")
    pd.DataFrame(rows, columns=ROW_COLS).to_csv(RESULTS / "claims_published_signatures.csv", index=False)
    (RESULTS / "notes_published_signatures.json").write_text(json.dumps(notes, indent=1, default=float))
    print(pd.DataFrame(rows)[["key", "method", "n_a", "n_b", "logrank_p", "hr_a_vs_b", "ari_vs_paper"]].to_string())


if __name__ == "__main__":
    main()
