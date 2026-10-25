"""Claim 3 (ISA part): apply canine-derived signatures to TARGET-OS.

Signatures applied, all fixed before this script reads any human data:
  published ISA-1..4    Supplementary Table 3 (canine symbols -> human orthologs)
  denovo-ISA-*          results/frozen_signatures.json from scripts/03_denovo_dog.py;
                        the file's sha256 is checked and the run stops if it changed
Each signature: k-means (silhouette over k = 2..6, plus the 2-cluster solution) in
TARGET all and TARGET non-metastatic, then log-rank on PFS and OS.
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
from cos_sig.clustering import seed_stability  # noqa: E402
from cos_sig.orthologs import map_signature  # noqa: E402
from cos_sig.paths import RESULTS, work_dir  # noqa: E402
from cos_sig.pipeline import ROW_COLS, cluster_on_signature, survival_row  # noqa: E402
from cos_sig.signatures import PUBLISHED_ISA, load_frozen  # noqa: E402
from cos_sig.survival import median_split  # noqa: E402

warnings.filterwarnings("ignore")
W = work_dir()
SEED = 0
N_BOOT = 1000


def main() -> None:
    frozen = load_frozen(RESULTS / "frozen_signatures.json")  # raises if edited after freezing
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    ht = pd.read_pickle(W / "hum_tpm_qn.pkl")
    hm = pd.read_pickle(W / "hum_meta.pkl")
    cohorts = {"TARGET all": ("hum", hm.index),
               "TARGET non-metastatic": ("hum_nonmet", hm.index[hm.metastatic == False])}  # noqa: E712
    sigs = {**PUBLISHED_ISA, **frozen}
    rows, notes, assign = [], {}, pd.DataFrame(index=hm.index)
    for sname, dog_genes in sigs.items():
        hmap = map_signature(dog_genes, omap)
        genes = [h for h in hmap.values() if h]
        for cname, (pref, idx) in cohorts.items():
            res = cluster_on_signature(ht[idx], genes, seed=SEED, hi="high", lo="low")
            lab = res["labels"]
            tag = f"{pref}_{sname.lower().replace('-', '')}"
            assign[f"{sname}_{pref}"] = lab
            stab = seed_stability(res["z"], 2)
            notes[tag] = {"canine_genes": len(dog_genes), "mapped": len(genes),
                          "present": len(res["present"]), "unmapped": [g for g, h in hmap.items() if not h],
                          "k_silhouette": res["k_silhouette"], "silhouettes": res["silhouettes"],
                          "sizes_k2": lab.value_counts().to_dict(), "seed_stability_k2": stab}
            pi = sname.split("-")[-1] if sname.startswith("ISA-") else None
            for ep, tcol in [("PFS", "pfs"), ("OS", "os")]:
                key = f"{tag}_{ep.lower()}"
                variants = [("rederived", lab, "low", "high")]
                if pi in {"1", "2", "3"} and ep == "PFS":
                    col = f"paper_isa{pi}_all" if cname == "TARGET all" else f"paper_isa{pi}_nonmet"
                    pl = hm.loc[idx, col].astype(int).astype(str)
                    variants.insert(0, ("paper_labels", pl, "1", "0"))
                for method, groups, a, b in variants:
                    paper = None
                    if method == "rederived" and pi in {"1", "2", "3"}:
                        col = f"paper_isa{pi}_all" if cname == "TARGET all" else f"paper_isa{pi}_nonmet"
                        paper = hm.loc[idx, col]
                    rows.append(survival_row(key, 3, cname, sname, ep, method, groups, hm.loc[idx, tcol],
                                             hm.loc[idx, f"{tcol}_event"], a, b,
                                             k_sil=res["k_silhouette"] if method == "rederived" else np.nan,
                                             genes_used=len(res["present"]), paper_labels=paper,
                                             n_boot=N_BOOT, seed=SEED))
    # IRX1 median split in non-metastatic TARGET (Supplementary Fig. 9a)
    idx = cohorts["TARGET non-metastatic"][1]
    irx = median_split(np.log2(ht.loc["IRX1", idx] + 2))
    rows.append(survival_row("hum_nonmet_irx1_pfs", 3, "TARGET non-metastatic", "IRX1 median split", "PFS",
                             "rederived", irx, hm.loc[idx, "pfs"], hm.loc[idx, "pfs_event"], "high", "low",
                             genes_used=1, n_boot=N_BOOT, seed=SEED))
    assign.index.name = "target_case"
    assign.to_csv(RESULTS / "cluster_assignments_target_isa.csv")
    pd.DataFrame(rows, columns=ROW_COLS).to_csv(RESULTS / "claims_transfer_target.csv", index=False)
    (RESULTS / "notes_transfer_target.json").write_text(json.dumps(notes, indent=1, default=float))
    print(pd.DataFrame(rows)[["key", "method", "n_a", "n_b", "k_silhouette", "logrank_p", "hr_a_vs_b",
                              "ari_vs_paper"]].to_string())


if __name__ == "__main__":
    main()
