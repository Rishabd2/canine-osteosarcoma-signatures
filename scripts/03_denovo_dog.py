"""Claim 2: de novo signatures in DOG2 (ssGSEA -> k-means -> limma-voom DEGs -> ISA).

This script reads canine data only. It ends by freezing the de novo signatures
(results/frozen_signatures.json, with a sha256 of the gene content) before any
human data is touched; scripts/04_transfer_target.py refuses to run if that file
was edited afterwards.
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
from cos_sig.clustering import kmeans_labels, kmeans_silhouette  # noqa: E402
from cos_sig.orthologs import collapse_rows, strip_suffix  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.pipeline import ROW_COLS, survival_row  # noqa: E402
from cos_sig.rrun import run_r  # noqa: E402
from cos_sig.scoring import read_gmt, ssgsea  # noqa: E402
from cos_sig.signatures import PUBLISHED_ISA, freeze, published_gs  # noqa: E402
from cos_sig.survival import median_split  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
SEED = 0
N_BOOT = 1000
ISA_RUNS = 25
IMMUNE_SETS = ["HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
               "HALLMARK_INFLAMMATORY_RESPONSE"]


def jaccard(a, b) -> float:
    a, b = set(a), set(b)
    return len(a & b) / len(a | b) if a | b else 0.0


def consensus_biclusters(runs: pd.DataFrame, n_runs: int, min_frac: float = 0.75, jac: float = 0.5) -> list[dict]:
    """Group biclusters from repeated ISA runs by gene Jaccard >= ``jac`` (single linkage).

    A group is kept if it appears in more than ``min_frac`` of runs (the paper's 75%
    rule). Consensus genes / samples are those present in at least half of the
    group's members.
    """
    items = [(int(r.run), set(str(r.genes).split(";")), set(str(r.samples).split(";"))) for r in runs.itertuples()]
    n = len(items)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(n):
        for j in range(i + 1, n):
            if jaccard(items[i][1], items[j][1]) >= jac:
                parent[find(i)] = find(j)
    groups: dict[int, list[int]] = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    out = []
    for members in groups.values():
        run_ids = {items[i][0] for i in members}
        frac = len(run_ids) / n_runs
        gcount = pd.Series([g for i in members for g in items[i][1]]).value_counts()
        scount = pd.Series([s for i in members for s in items[i][2]]).value_counts()
        out.append({"runs": len(run_ids), "frac_runs": frac, "kept": frac > min_frac,
                    "genes": sorted(gcount[gcount >= len(members) / 2].index),
                    "samples": sorted(scount[scount >= len(members) / 2].index)})
    return sorted(out, key=lambda d: (-d["frac_runs"], -len(d["genes"])))


def main() -> None:
    meta = pd.read_pickle(W / "dog_meta.pkl")
    meta.index = meta.index.astype(str)
    counts = pd.read_pickle(W / "dog_counts.pkl")
    counts.columns = counts.columns.astype(str)
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    rows, notes = [], {}

    # 1. authors' TMM/ComBat log-CPM (R/authors_tmm_combat.R)
    cpath, bpath, opath = W / "dog_counts.csv", W / "dog_batch.csv", W / "dog_logcpm_combat.csv"
    counts.to_csv(cpath)
    pd.DataFrame({"Batch.ID": meta.trial}).to_csv(bpath)
    if not opath.exists():
        print(run_r("authors_tmm_combat.R", cpath, bpath, opath).strip())
    lc = pd.read_csv(opath, index_col=0)
    lc.columns = lc.columns.astype(str)
    lc_sym = lc.groupby(lc.index.map(strip_suffix)).mean()
    lc_sym = lc_sym[~lc_sym.index.str.startswith("LOC")]
    notes["logcpm_probes"] = int(lc.shape[0])
    notes["logcpm_genes_after_collapse"] = int(lc_sym.shape[0])

    # 2. ssGSEA on Hallmarks (canine symbols mapped to human orthologs)
    hmap = omap.set_index("dog_symbol")["human_symbol"]
    lc_h = collapse_rows(lc_sym, pd.Series(hmap.reindex(lc_sym.index).values, index=lc_sym.index))
    hall = read_gmt(D / "ref" / "h.all.v7.5.1.symbols.gmt")
    ss = ssgsea(lc_h, hall, alpha=0.75)
    ss.to_csv(RESULTS / "ssgsea_hallmark_dog.csv.gz", float_format="%.5f",
              compression={"method": "gzip", "mtime": 0})
    notes["ssgsea_sets_scored"] = int(ss.shape[0])
    zs = ((ss.T - ss.T.mean()) / ss.T.std(ddof=0))
    sil = kmeans_silhouette(zs, seed=SEED)
    lab2 = kmeans_labels(zs, 2, seed=SEED)
    immune = zs[IMMUNE_SETS].mean(axis=1).groupby(lab2).mean()
    fp_cluster = immune.idxmax()
    grp = lab2.map(lambda c: "FP" if c == fp_cluster else "PP")
    notes["ssgsea_kmeans"] = {"k_silhouette": sil.k, "silhouettes": sil.silhouettes,
                              "sizes": grp.value_counts().to_dict(),
                              "ari_vs_paper": adjusted_rand_score(meta.paper_ssgsea, grp)}
    for ep, tcol in [("DFI", "dfi"), ("OS", "os")]:
        key = f"dog_ssgsea_{ep.lower()}"
        rows.append(survival_row(key, 2, "DOG2", "ssGSEA hallmark k-means", ep, "paper_labels", meta.paper_ssgsea,
                                 meta[tcol], meta[f"{tcol}_event"], "PP", "FP", n_boot=N_BOOT, seed=SEED))
        rows.append(survival_row(key, 2, "DOG2", "ssGSEA hallmark k-means", ep, "rederived", grp, meta[tcol],
                                 meta[f"{tcol}_event"], "PP", "FP", k_sil=sil.k, genes_used=ss.shape[0],
                                 paper_labels=meta.paper_ssgsea, n_boot=N_BOOT, seed=SEED))

    # 3-4. limma-voom DEGs between ssGSEA clusters (trial as batch), then ISA on the DEGs.
    # Run twice: on our ssGSEA clusters ("denovo") and on the paper's published ssGSEA
    # clusters ("paperlab"), which separates clustering differences from DEG/ISA differences.
    gs = published_gs(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx")
    denovo = {}
    for tag, labels in [("denovo", grp), ("paperlab", meta.paper_ssgsea)]:
        lpath, dpath = W / f"dog_ssgsea_labels_{tag}.csv", W / f"dog_ssgsea_limma_{tag}.csv"
        pd.DataFrame({"Cluster": (labels == "FP").astype(int), "Batch.ID": meta.trial}).to_csv(lpath)
        print(run_r("limma_voom.R", cpath, lpath, dpath, "batch").strip())
        tt = pd.read_csv(dpath, index_col=0)
        # The paper reports 259 DEGs for this contrast; with n = 186 an adj.P-only cutoff returns
        # thousands, so we apply the Methods' DEG threshold (BH adj.P < 0.05 and fold change > 3).
        degs_probe = tt.index[(tt["adj.P.Val"] < 0.05) & (tt.logFC.abs() > np.log2(3))]
        degs = sorted({g for g in map(strip_suffix, degs_probe) if not g.startswith("LOC")})
        notes[f"ssgsea_degs_{tag}"] = {"n_adjP_lt_0.05_probes": int((tt["adj.P.Val"] < 0.05).sum()),
                                       "n_adjP_and_absFC_gt_3_probes": int(len(degs_probe)), "n_genes": len(degs),
                                       "overlap_GS1": len(set(degs) & set(gs["GS-1"])),
                                       "overlap_GS2": len(set(degs) & set(gs["GS-2"]))}
        pd.Series(degs, name="gene").to_csv(RESULTS / f"dog_ssgsea_degs_{tag}.csv", index=False)
        xpath = W / f"dog_isa_input_{tag}.csv"
        lc_sym.loc[[g for g in degs if g in lc_sym.index]].T.to_csv(xpath)
        print(run_r("isa_runs.R", xpath, W / f"dog_{tag}", ISA_RUNS, SEED).strip().splitlines()[-1])
        runs = pd.read_csv(W / f"dog_{tag}_isa_runs.csv")
        cons = consensus_biclusters(runs, ISA_RUNS)
        kept = [c for c in cons if c["kept"]]
        notes[f"isa_{tag}"] = {"runs": ISA_RUNS, "biclusters_total": int(len(runs)), "groups": len(cons),
                               "groups_kept_gt_75pct": len(kept), "kept": []}
        for i, c in enumerate(kept, start=1):
            name = f"{tag}-ISA-{i}"
            denovo[name] = c["genes"]
            best = max(PUBLISHED_ISA, key=lambda p: len(set(c["genes"]) & set(PUBLISHED_ISA[p])))
            notes[f"isa_{tag}"]["kept"].append({
                "name": name, "genes": len(c["genes"]), "samples": len(c["samples"]), "frac_runs": c["frac_runs"],
                "closest_published": best, "overlap_with_closest": len(set(c["genes"]) & set(PUBLISHED_ISA[best])),
                "jaccard_with_closest": jaccard(c["genes"], PUBLISHED_ISA[best]), "genes_list": c["genes"]})
            grp_isa = pd.Series("REST", index=meta.index)
            grp_isa[[s for s in c["samples"] if s in grp_isa.index]] = "ISA"
            for ep, tcol in [("DFI", "dfi"), ("OS", "os")]:
                rows.append(survival_row(f"dog_{name.lower().replace('-', '_')}_{ep.lower()}", 2, "DOG2",
                                         f"{name} bicluster vs rest", ep, "rederived", grp_isa, meta[tcol],
                                         meta[f"{tcol}_event"], "ISA", "REST", genes_used=len(c["genes"]),
                                         n_boot=N_BOOT, seed=SEED))

    # 5. the paper's own ISA sample biclusters (Supplementary Data 1, Fig. 7), re-tested
    for i in range(1, 5):
        for ep, tcol in [("DFI", "dfi"), ("OS", "os")]:
            key = f"dog_isa{i}" if ep == "DFI" else f"dog_isa{i}_os"
            rows.append(survival_row(key, 2, "DOG2", f"ISA-{i} bicluster vs rest", ep, "paper_labels",
                                     meta[f"paper_isa{i}"], meta[tcol], meta[f"{tcol}_event"], "ISA", "REST",
                                     n_boot=N_BOOT, seed=SEED))

    # 6. IRX1 median split (Supplementary Fig. 9b)
    tq = pd.read_pickle(W / "dog_tpm_qn.pkl")
    tq.columns = tq.columns.astype(str)
    irx = median_split(np.log2(tq.loc["IRX1"] + 2))
    rows.append(survival_row("dog_irx1_dfi", 2, "DOG2", "IRX1 median split", "DFI", "rederived", irx,
                             meta.dfi, meta.dfi_event, "high", "low", genes_used=1, n_boot=N_BOOT, seed=SEED))

    # 7. freeze: de novo signatures (canine symbols) before any human data is read
    digest = freeze(denovo, RESULTS / "frozen_signatures.json",
                    {"script": "scripts/03_denovo_dog.py", "inputs": ["GSE238110 counts", "DOG2 outcomes"],
                     "seed": SEED, "isa_runs": ISA_RUNS, "note": "derived from canine data only"})
    notes["frozen_sha256"] = digest
    pd.DataFrame(rows, columns=ROW_COLS).to_csv(RESULTS / "claims_denovo_dog.csv", index=False)
    (RESULTS / "notes_denovo_dog.json").write_text(json.dumps(notes, indent=1, default=float))
    print(json.dumps(notes, indent=1, default=float)[:6000])
    print(pd.DataFrame(rows)[["key", "method", "n_a", "n_b", "logrank_p", "hr_a_vs_b", "ari_vs_paper"]].to_string())


if __name__ == "__main__":
    main()
