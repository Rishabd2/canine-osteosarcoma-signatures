"""Claim 4: shared pathway and immune programs between GS-1 FP and PP groups.

* Hallmark GSEA (phenotype permutation, signal-to-noise) FP vs PP in DOG2 and
  TARGET all / non-metastatic / metastatic, compared with the paper's Fig. 4 lists.
* limma-voom DEGs FP vs PP (paper threshold: BH adj.P < 0.05 and |FC| > 3) in both
  species, shared DEGs, and Reactome over-representation on the shared set
  (open substitute for IPA).
* Immune deconvolution with MCP-counter and TIL10 constrained least squares
  (open substitutes for CIBERSORTx/LM22); FP vs PP by two-sided Mann-Whitney U,
  BH within each dataset and method, compared with Table 1 directions.
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.deconv import load_mcp_markers, load_til10, mcp_counter, til10_nnls  # noqa: E402
from cos_sig.enrichment import gsea_phenotype, ora  # noqa: E402
from cos_sig.orthologs import collapse_rows, strip_suffix  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402
from cos_sig.rrun import run_r  # noqa: E402
from cos_sig.scoring import read_gmt  # noqa: E402
from cos_sig.survival import bh_adjust  # noqa: E402

warnings.filterwarnings("ignore")
D, W = data_dir(), work_dir()
SEED = 0
PERM = 1000
FC = np.log2(3)

# Paper Table 1 (CIBERSORTx LM22): direction of the FP vs PP difference and BH p in each species
PAPER_TABLE1 = {
    "T cells CD8": {"dog": ("FP>PP", 2.66e-08), "human": ("FP>PP", 0.011562)},
    "Macrophages M0": {"dog": ("PP>FP", 6.63e-08), "human": ("PP>FP", 1.20e-06)},
    "Macrophages M1": {"dog": ("FP>PP", 4.97e-06), "human": ("FP>PP", 6.69e-05)},
    "Macrophages M2": {"dog": ("FP>PP", 1.47e-06), "human": ("FP>PP", 6.69e-05)},
}
# Our cell types that correspond to the Table 1 rows
CELL_MATCH = {
    ("TIL10", "T.cells.CD8"): "T cells CD8",
    ("TIL10", "Macrophages.M1"): "Macrophages M1",
    ("TIL10", "Macrophages.M2"): "Macrophages M2",
    ("MCP", "CD8 T cells"): "T cells CD8",
    ("MCP", "Monocytic lineage"): "Macrophages (all)",
}


def paper_gsea_sets(xlsx: Path) -> dict[str, set[str]]:
    df = pd.read_excel(xlsx, sheet_name="Figure 4", header=None)
    def norm(s):
        s = str(s).upper().strip().replace(" ", "_").replace("IS2", "IL2").replace("PHYORYLATION", "PHORYLATION")
        s = s.replace("IL6_JAK_STAT5", "IL6_JAK_STAT3").replace("KRAS_SIGNALING_UP", "KRAS_SIGNALING_UP")
        return s
    out = {}
    for name, (r0, r1, c) in {"DOG2": (3, 23, 1), "TARGET all": (3, 23, 8), "TARGET non-metastatic": (26, 46, 1),
                              "TARGET metastatic": (26, 46, 8)}.items():
        block = df.iloc[r0:r1, [c, c + 3]].dropna()
        out[name] = {norm(t) for t, q in block.itertuples(index=False) if float(q) < 0.05}
    return out


def canon(term: str) -> str:
    t = term.replace("HALLMARK_", "")
    for a, b in [("_RESPONSE", ""), ("_SIGNALING", ""), ("_VIA_NFKB", "")]:
        t = t.replace(a, b)
    return t


def main() -> None:
    omap = pd.read_pickle(W / "ortholog_map.pkl")
    hmap = omap.set_index("dog_symbol")["human_symbol"]
    da = pd.read_csv(RESULTS / "cluster_assignments_dog.csv", index_col=0)
    da.index = da.index.astype(str)
    ta = pd.read_csv(RESULTS / "cluster_assignments_target.csv", index_col=0)
    hm = pd.read_pickle(W / "hum_meta.pkl")
    dog_grp = da["GS-1_rederived"]
    hum_grp = {"TARGET all": ta["GS-1_hum_rederived"],
               "TARGET non-metastatic": ta["GS-1_hum_nonmet_rederived"].dropna(),
               "TARGET metastatic": ta["GS-1_hum_met_rederived"].dropna()}
    hall = read_gmt(D / "ref" / "h.all.v7.5.1.symbols.gmt")
    notes: dict = {}

    # expression on human symbols, log2(TPM-qn + 1)
    dq = pd.read_pickle(W / "dog_tpm_qn.pkl")
    dq.columns = dq.columns.astype(str)
    dq_h = collapse_rows(dq, pd.Series(hmap.reindex(dq.index).values, index=dq.index))
    hq = pd.read_pickle(W / "hum_tpm_qn.pkl")

    # ---------- GSEA ----------
    gsea_rows = []
    sets = {"DOG2": (np.log2(dq_h + 1), dog_grp)}
    sets.update({k: (np.log2(hq[g.index] + 1), g) for k, g in hum_grp.items()})
    for cname, (expr, grp) in sets.items():
        res = gsea_phenotype(expr, grp, hall, permutations=PERM, seed=SEED)
        res.insert(0, "cohort", cname)
        gsea_rows.append(res)
    gsea = pd.concat(gsea_rows)
    gsea.to_csv(RESULTS / "gsea_hallmark_fp_vs_pp.csv", index=False)
    paper = paper_gsea_sets(D / "supp" / "42003_2023_5208_MOESM4_ESM.xlsx")
    gs_notes = {}
    for cname in sets:
        ours = {canon(t) for t in gsea[(gsea.cohort == cname) & (gsea.nes > 0) & (gsea.fdr < 0.05)].term}
        theirs = {canon("HALLMARK_" + t) for t in paper[cname]}
        gs_notes[cname] = {"ours_fp_sig": len(ours), "paper_fp_sig": len(theirs), "overlap": len(ours & theirs),
                           "ours": sorted(ours), "paper": sorted(theirs),
                           "ours_pp_sig": int(((gsea.cohort == cname) & (gsea.nes < 0) & (gsea.fdr < 0.05)).sum())}
    dog_sig = set(gs_notes["DOG2"]["ours"])
    gs_notes["dog_fp_sets_also_sig_in_TARGET_all"] = len(dog_sig & set(gs_notes["TARGET all"]["ours"]))
    notes["gsea"] = gs_notes

    # ---------- DEGs and shared pathways ----------
    dc = pd.read_pickle(W / "dog_counts.pkl")
    dc.columns = dc.columns.astype(str)
    dm = pd.read_pickle(W / "dog_meta.pkl")
    dm.index = dm.index.astype(str)
    cpath, lpath, opath = W / "dog_counts.csv", W / "dog_gs1_labels.csv", W / "dog_gs1_limma.csv"
    if not cpath.exists():
        dc.to_csv(cpath)
    pd.DataFrame({"Cluster": (dog_grp == "FP").astype(int), "Batch.ID": dm.trial}).to_csv(lpath)
    run_r("limma_voom.R", cpath, lpath, opath, "batch")
    dt = pd.read_csv(opath, index_col=0)
    hc = pd.read_pickle(W / "hum_counts.pkl")
    deg_sets = {}
    dsig = dt[(dt["adj.P.Val"] < 0.05) & (dt.logFC.abs() > FC)]
    deg_sets["DOG2"] = {h for h in hmap.reindex([strip_suffix(g) for g in dsig.index]).dropna()}
    for cname, g in hum_grp.items():
        lp, op = W / f"hum_gs1_labels_{cname[7:10]}.csv", W / f"hum_gs1_limma_{cname[7:10]}.csv"
        pd.DataFrame({"Cluster": (g == "FP").astype(int)}).to_csv(lp)
        sub = W / f"hum_counts_{cname[7:10]}.csv"
        hc[g.index].to_csv(sub)
        run_r("limma_voom.R", sub, lp, op, "nobatch")
        tt = pd.read_csv(op, index_col=0)
        deg_sets[cname] = set(tt[(tt["adj.P.Val"] < 0.05) & (tt.logFC.abs() > FC)].index)
    shared = deg_sets["DOG2"] & deg_sets["TARGET all"] & deg_sets["TARGET non-metastatic"]
    universe = set(hq.index) & set(dq_h.index)
    reactome = read_gmt(D / "ref" / "c2.cp.reactome.v7.5.1.symbols.gmt")
    ora_res = ora(shared, universe, reactome)
    ora_res.to_csv(RESULTS / "reactome_ora_shared_degs.csv", index=False)
    pd.Series(sorted(shared), name="gene").to_csv(RESULTS / "shared_degs_gs1.csv", index=False)
    notes["degs"] = {k: len(v) for k, v in deg_sets.items()}
    notes["degs"].update({"shared_three_sets": len(shared), "paper_dog": 317, "paper_target_all": 281,
                          "paper_target_nonmet": 321, "paper_shared": 84,
                          "reactome_fdr_lt_0.05": int((ora_res.fdr < 0.05).sum()),
                          "reactome_top": ora_res.head(10)[["term", "overlap", "fdr"]].to_dict("records")})

    # ---------- immune deconvolution ----------
    mcp = load_mcp_markers(D / "ref" / "MCPcounter_genes.txt")
    til = load_til10(D / "ref" / "TIL10_signature.txt")
    dtpm = pd.read_pickle(W / "dog_tpm.pkl")
    dtpm.columns = dtpm.columns.astype(str)
    dtpm_h = collapse_rows(dtpm, pd.Series(hmap.reindex(dtpm.index).values, index=dtpm.index), how="sum")
    dtpm_h = dtpm_h / dtpm_h.sum() * 1e6
    htpm = pd.read_pickle(W / "hum_tpm.pkl")
    decon = {"dog": (dtpm_h, dog_grp), "human": (htpm, hum_grp["TARGET all"])}
    drows = []
    for sp, (tp, grp) in decon.items():
        for meth, est in [("MCP", mcp_counter(np.log2(tp + 1), mcp)), ("TIL10", til10_nnls(tp, til))]:
            est.to_csv(RESULTS / f"deconv_{meth.lower()}_{sp}.csv.gz", float_format="%.6f",
                       compression={"method": "gzip", "mtime": 0})
            sub = []
            for cell in est.index:
                fp = est.loc[cell, grp.index[grp == "FP"]]
                pp = est.loc[cell, grp.index[grp == "PP"]]
                p = mannwhitneyu(fp, pp, alternative="two-sided").pvalue if fp.nunique() + pp.nunique() > 2 else np.nan
                sub.append({"species": sp, "method": meth, "cell": cell, "mean_fp": fp.mean(), "mean_pp": pp.mean(),
                            "direction": "FP>PP" if fp.median() > pp.median() else "PP>FP", "p": p})
            sub = pd.DataFrame(sub)
            sub["q_bh"] = bh_adjust(sub.p.to_numpy())
            drows.append(sub)
    dec = pd.concat(drows)
    dec["paper_row"] = [CELL_MATCH.get((m, c), "") for m, c in zip(dec.method, dec.cell)]
    dec.to_csv(RESULTS / "deconvolution_fp_vs_pp.csv", index=False)
    notes["deconvolution_genes"] = {
        "til10_genes_dog": int(len(set(til.index) & set(dtpm_h.index))),
        "til10_genes_human": int(len(set(til.index) & set(htpm.index))), "til10_total": int(til.shape[0])}
    (RESULTS / "notes_immune_pathways.json").write_text(json.dumps(notes, indent=1, default=float))
    print(json.dumps(notes, indent=1, default=float)[:4000])
    print(dec[dec.paper_row != ""].to_string())


if __name__ == "__main__":
    main()
