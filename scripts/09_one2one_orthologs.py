"""Sensitivity: use only Ensembl one-to-one orthologs for every cross-species step.

The main run also maps one-to-many orthologs (by rule) and identical symbols with no
Ensembl row. Here those genes are dropped, and two things are re-run:
  1. every signature applied to TARGET (published GS-1, GS-2, ISA-1..4 and the frozen
     canine signatures), in TARGET all, non-metastatic and (GS only) metastatic;
  2. the canine ssGSEA k-means step, whose Hallmark scoring needs human symbols.
Writes results/sensitivity_one2one_orthologs.csv.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.clustering import kmeans_labels  # noqa: E402
from cos_sig.orthologs import collapse_rows, map_signature, strip_suffix  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.pipeline import cluster_on_signature  # noqa: E402
from cos_sig.scoring import read_gmt, ssgsea  # noqa: E402
from cos_sig.signatures import PUBLISHED_ISA, load_frozen, published_gs  # noqa: E402
from cos_sig.survival import align_survival, logrank_p  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
SEED = 0
IMMUNE_SETS = ["HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_INFLAMMATORY_RESPONSE"]


def lr(groups, time, event) -> float:
    return logrank_p(align_survival(groups, time, event))


def main() -> None:
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    o2o = omap[omap.method == "ensembl_one2one"]
    gs = published_gs(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx")
    sigs = {"GS-1": gs["GS-1-human"], "GS-2": gs["GS-2"], **PUBLISHED_ISA,
            **load_frozen(RESULTS / "frozen_signatures.json")}
    ht = pd.read_pickle(W / "hum_tpm_qn.pkl")
    hm = pd.read_pickle(W / "hum_meta.pkl")
    cohorts = {"TARGET all": hm.index, "TARGET non-metastatic": hm.index[hm.metastatic == False],  # noqa: E712
               "TARGET metastatic": hm.index[hm.metastatic == True]}  # noqa: E712
    rows = []
    for sname, genes in sigs.items():
        full = [h for h in map_signature(genes, omap).values() if h]
        strict = [h for h in map_signature(genes, o2o).values() if h]
        for cname, idx in cohorts.items():
            if cname == "TARGET metastatic" and not sname.startswith("GS"):
                continue
            a = cluster_on_signature(ht[idx], full, seed=SEED, hi="high", lo="low")
            b = cluster_on_signature(ht[idx], strict, seed=SEED, hi="high", lo="low")
            for ep in ["pfs", "os"]:
                rows.append({"step": "transfer to TARGET", "signature": sname, "cohort": cname,
                             "endpoint": ep.upper(), "genes_all_mappings": len(a["present"]),
                             "genes_one2one": len(b["present"]),
                             "p_all_mappings": lr(a["labels"], hm.loc[idx, ep], hm.loc[idx, f"{ep}_event"]),
                             "p_one2one": lr(b["labels"], hm.loc[idx, ep], hm.loc[idx, f"{ep}_event"]),
                             "ari_between_clusterings": adjusted_rand_score(a["labels"], b["labels"])})

    # canine ssGSEA step with one-to-one orthologs only
    meta = pd.read_pickle(W / "dog_meta.pkl")
    meta.index = meta.index.astype(str)
    lc = pd.read_csv(W / "dog_logcpm_combat.csv", index_col=0)
    lc.columns = lc.columns.astype(str)
    lc_sym = lc.groupby(lc.index.map(strip_suffix)).mean()
    lc_sym = lc_sym[~lc_sym.index.str.startswith("LOC")]
    hall = read_gmt(D / "ref" / "h.all.v7.5.1.symbols.gmt")
    labs, ngenes = {}, {}
    for tag, m in [("all mappings", omap), ("one2one", o2o)]:
        hmap = m.set_index("dog_symbol")["human_symbol"]
        lc_h = collapse_rows(lc_sym, pd.Series(hmap.reindex(lc_sym.index).values, index=lc_sym.index))
        ss = ssgsea(lc_h, hall, alpha=0.75)
        zs = (ss.T - ss.T.mean()) / ss.T.std(ddof=0)
        lab = kmeans_labels(zs, 2, seed=SEED)
        fp = zs[IMMUNE_SETS].mean(axis=1).groupby(lab).mean().idxmax()
        labs[tag] = lab.map(lambda c: "FP" if c == fp else "PP")
        ngenes[tag] = lc_h.shape[0]
    for ep in ["dfi", "os"]:
        rows.append({"step": "DOG2 ssGSEA k-means", "signature": "Hallmark ssGSEA", "cohort": "DOG2",
                     "endpoint": ep.upper(), "genes_all_mappings": ngenes["all mappings"],
                     "genes_one2one": ngenes["one2one"],
                     "p_all_mappings": lr(labs["all mappings"], meta[ep], meta[f"{ep}_event"]),
                     "p_one2one": lr(labs["one2one"], meta[ep], meta[f"{ep}_event"]),
                     "ari_between_clusterings": adjusted_rand_score(labs["all mappings"], labs["one2one"])})
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS / "sensitivity_one2one_orthologs.csv", index=False)
    changed = df[(df.p_all_mappings < 0.05) != (df.p_one2one < 0.05)]
    print(df.round(4).to_string())
    print(f"\n{len(changed)} of {len(df)} tests change their call at p < 0.05")


if __name__ == "__main__":
    main()
