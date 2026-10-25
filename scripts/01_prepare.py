"""Build analysis-ready matrices and clinical tables for DOG2 and TARGET-OS.

Outputs (uncommitted, in $COS_WORK_DIR):
  dog_counts.pkl        raw counts, 186 primary tumours, canine symbols (suffix kept)
  dog_tpm_qn.pkl        filterByExpr genes, TPM, duplicate symbols averaged, quantile normalized
  dog_meta.pkl          pid, trial (batch), DFI/OS from the paper's source data
  hum_counts.pkl        TARGET-OS STAR unstranded counts, protein-coding, 79 samples
  hum_tpm_qn.pkl        TPM (protein-coding), quantile normalized
  hum_tpm.pkl           TPM before quantile normalization (for deconvolution)
  hum_meta.pkl          PFS/OS from the paper's source data, metastatic status, GDC cross-check
  ortholog_map.pkl      canine symbol -> human symbol
Committed summaries: results/data_summary.json, results/ortholog_mapping_summary.csv
"""
from __future__ import annotations

import gzip
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.normalize import filter_by_expr, gene_lengths_from_gtf, quantile_normalize, tpm  # noqa: E402
from cos_sig.orthologs import build_dog_to_human_map, load_ensembl_orthologs, strip_suffix  # noqa: E402
from cos_sig.paths import RESULTS, data_dir, work_dir  # noqa: E402

D = data_dir()
W = work_dir()
SRC = D / "supp" / "42003_2023_5208_MOESM4_ESM.xlsx"


def dog() -> dict:
    counts = pd.read_csv(D / "geo" / "GSE238110_RawCountFile_combined.csv.gz", index_col=0, encoding="utf-8-sig")
    prim = [c for c in counts.columns if re.match(r"^\d+_\d{4}_.*[Tt]umor$", c) and "_Me" not in c]
    pid = {c: int(re.match(r"^\d+_(\d+)", c).group(1)) for c in prim}
    clin = pd.read_excel(SRC, sheet_name="Figure 2", index_col=0)
    clin.index = clin.index.astype(int)
    keep = [c for c in prim if pid[c] in clin.index]
    counts = counts[keep].rename(columns=pid)
    counts = counts[sorted(counts.columns)]
    # trial (batch) from GEO sample titles, e.g. "10_0512_Tumor [COTC022-0512-T]"
    titles = [ln for ln in gzip.open(D / "geo" / "GSE238110_series_matrix.txt.gz", "rt")
              if ln.startswith("!Sample_title")][0].rstrip("\n").split("\t")[1:]
    trial = {}
    for t in titles:
        m = re.search(r"(COTC\d+)-(\d+)", t)
        if m:
            trial[int(m.group(2))] = m.group(1)
    arm = pd.read_excel(SRC, sheet_name="Figure S1", index_col=0)["Study"]
    meta = pd.DataFrame({
        "trial": pd.Series(trial).reindex(counts.columns),
        "arm": arm.reindex(counts.columns),
        "dfi": clin["DFI"].reindex(counts.columns).astype(float),
        "dfi_event": clin["DFI.Censor"].reindex(counts.columns).astype(int),
        "os": clin["OSV"].reindex(counts.columns).astype(float),
        "os_event": clin["Surv.censor"].reindex(counts.columns).astype(int),
        "paper_gs1": clin["PROGNOSIS 2b"].reindex(counts.columns).map({"FAVORABLE": "FP", "POOR": "PP"}),
        "paper_gs2": clin["PROGNOSIS 2e"].reindex(counts.columns).map({"FAVORABLE": "FP", "POOR": "PP"}),
    })
    meta.index.name = "pid"
    s11 = pd.read_excel(SRC, sheet_name="Figure S11", index_col=0)
    meta["paper_ssgsea"] = s11["Prognosis"].reindex(meta.index).map({"Favorable": "FP", "Poor": "PP"})
    f7 = pd.read_excel(SRC, sheet_name="Figure 7", index_col=0)
    for i, col in enumerate(["7-e", "7-f", "7-g", "7-h"], start=1):
        meta[f"paper_isa{i}"] = f7[col].reindex(meta.index).map(lambda v: "ISA" if str(v).startswith("ISA") else "REST")

    keep_rows = filter_by_expr(counts, min_count=10, min_prop=0.5, large_n=5)
    lengths = gene_lengths_from_gtf(D / "ref" / "GCF_011100685.1_UU_Cfam_GSD_1.0_genomic.gtf.gz")
    filt = counts[keep_rows]
    sym = filt.index.map(strip_suffix)
    has_len = sym.isin(lengths.index)
    filt_len = filt[has_len]
    lens = pd.Series(lengths.reindex(sym[has_len]).to_numpy(), index=filt_len.index)
    t = tpm(filt_len, lens)
    t = t.groupby(t.index.map(strip_suffix)).mean()
    t = t[~t.index.str.startswith("LOC")]
    tq = quantile_normalize(t)
    counts.to_pickle(W / "dog_counts.pkl")
    tq.to_pickle(W / "dog_tpm_qn.pkl")
    t.to_pickle(W / "dog_tpm.pkl")
    meta.to_pickle(W / "dog_meta.pkl")
    return {
        "dog_primary_tumour_columns_in_geo": len(prim),
        "dog_samples_matched_to_source_data": int(counts.shape[1]),
        "dog_probes_total": int(counts.shape[0]),
        "dog_probes_after_filterByExpr": int(keep_rows.sum()),
        "dog_probes_without_gtf_length": int((~has_len).sum()),
        "dog_genes_after_collapse_and_LOC_removal": int(tq.shape[0]),
        "dog_trial_counts": meta.trial.value_counts().to_dict(),
        "dog_trial_equals_arm": bool((meta.trial.map({"COTC021": "Rapa", "COTC022": "SOC"}) == meta.arm).all()),
        "dog_dfi_events": int(meta.dfi_event.sum()),
        "dog_os_events": int(meta.os_event.sum()),
    }


def human() -> dict:
    man = pd.DataFrame(json.loads((D / "gdc" / "file_manifest.json").read_text()))
    man["short"] = man["case"].str.replace("TARGET-40-", "", regex=False)
    f3 = pd.read_excel(SRC, sheet_name="Figure 3", index_col=0)
    tpms, cnts = {}, {}
    gene_info = None
    for _, r in man.iterrows():
        d = pd.read_csv(D / "gdc" / "star_counts" / f"{r.file_id}.tsv", sep="\t", comment="#")
        d = d[d.gene_id.str.startswith("ENSG")]
        if gene_info is None:
            gene_info = d[["gene_id", "gene_name", "gene_type"]].set_index("gene_id")
        tpms[r.short] = d.set_index("gene_id")["tpm_unstranded"]
        cnts[r.short] = d.set_index("gene_id")["unstranded"]
    tp = pd.DataFrame(tpms)
    ct = pd.DataFrame(cnts)
    pc = gene_info.gene_type.isin(["protein_coding", "polymorphic_pseudogene"])
    names = gene_info.gene_name
    tp_pc = tp[pc].groupby(names[pc]).sum()
    ct_pc = ct[pc].groupby(names[pc]).sum()
    samples = [s for s in f3.index if s in tp_pc.columns]
    tp_pc, ct_pc = tp_pc[samples], ct_pc[samples]
    tq = quantile_normalize(tp_pc)
    s4 = pd.read_excel(SRC, sheet_name="Figure S4", index_col=0)
    f9 = pd.read_excel(SRC, sheet_name="Figure 9", index_col=0)
    meta = pd.DataFrame({
        "pfs": f3["PFS"].astype(float), "pfs_event": f3["PFS.censor"].astype(int),
        "os": f3["Survival"].astype(float), "os_event": f3["Surv.censor"].astype(int),
        "paper_gs1": f3["PROGNOSIS 3b"].map({"FAVORABLE": "FP", "POOR": "PP"}),
        "paper_gs2": f3["PROGNOSIS 3e"].map({"FAVORABLE": "FP", "POOR": "PP"}),
    }).reindex(samples)
    meta["metastatic"] = [True if s in s4.index else (False if s in f9.index else np.nan) for s in samples]
    meta["paper_gs1_met"] = s4["S4-b"].reindex(samples)
    meta["paper_gs2_met"] = s4["s4-e"].reindex(samples)
    f8 = pd.read_excel(SRC, sheet_name="Figure 8", index_col=0)
    for i, col in enumerate([" 8-d", "8-e", "8-f"], start=1):
        meta[f"paper_isa{i}_all"] = f8[col].reindex(samples)
    for i, col in enumerate(["9-d", "9-e", "9-f"], start=1):
        meta[f"paper_isa{i}_nonmet"] = f9[col].reindex(samples)
    # independent check of metastatic status and OS against the GDC clinical records
    cases = json.loads((D / "gdc" / "clinical_cases.json").read_text())
    gmet, gvital = {}, {}
    for c in cases:
        s = c["submitter_id"].replace("TARGET-40-", "")
        dx = (c.get("diagnoses") or [{}])[0]
        gmet[s] = dx.get("metastasis_at_diagnosis")
        gvital[s] = (c.get("demographic") or {}).get("vital_status")
    meta["gdc_metastasis_at_diagnosis"] = pd.Series(gmet).reindex(samples)
    meta["gdc_vital_status"] = pd.Series(gvital).reindex(samples)
    meta.index.name = "sample"
    ct_pc.to_pickle(W / "hum_counts.pkl")
    tp_pc.to_pickle(W / "hum_tpm.pkl")
    tq.to_pickle(W / "hum_tpm_qn.pkl")
    meta.to_pickle(W / "hum_meta.pkl")
    gm = meta.gdc_metastasis_at_diagnosis.map({"Metastatic": True, "Non-metastatic (confirmed)": False,
                                               "Non-metastatic (Confirmed)": False, "No Metastasis": False,
                                               "Metastasis, NOS": True})
    agree_met = int((gm == meta.metastatic).sum())
    dead = meta.gdc_vital_status.eq("Dead").astype(int)
    return {
        "target_star_files": int(len(man)),
        "target_samples_used": len(samples),
        "target_samples_missing_from_gdc": [s for s in f3.index if s not in tp_pc.columns],
        "target_protein_coding_genes": int(tp_pc.shape[0]),
        "target_metastatic": int(meta.metastatic.sum()),
        "target_nonmetastatic": int((meta.metastatic == False).sum()),  # noqa: E712
        "target_gdc_metastasis_values": meta.gdc_metastasis_at_diagnosis.value_counts(dropna=False).to_dict(),
        "target_metastatic_agrees_with_gdc": agree_met,
        "target_os_event_agrees_with_gdc_vital_status": int((dead == meta.os_event).sum()),
        "target_pfs_events": int(meta.pfs_event.sum()),
        "target_os_events": int(meta.os_event.sum()),
    }


def orthologs(summary: dict) -> None:
    dog_counts = pd.read_pickle(W / "dog_counts.pkl")
    hum = pd.read_pickle(W / "hum_tpm_qn.pkl")
    orth = load_ensembl_orthologs(D / "ref" / "ensembl_dog_human_orthologs.tsv")
    sig = pd.read_excel(D / "supp" / "42003_2023_5208_MOESM5_ESM.xlsx", sheet_name="Gene Signature")
    dog_syms = sorted(set(dog_counts.index.map(strip_suffix)) | set(sig.stack().dropna().astype(str)))
    mp = build_dog_to_human_map(dog_syms, orth, hum.index)
    mp.to_pickle(W / "ortholog_map.pkl")
    tab = mp.groupby("method").size().rename("n_canine_symbols").reset_index()
    tab["note"] = tab["method"].map({
        "ensembl_one2one": "Ensembl one-to-one ortholog",
        "ensembl_best": "Ensembl one-to-many or many-to-many, rule in orthologs.py",
        "identity": "no Ensembl row, canine symbol is a TARGET symbol (DLA- read as HLA-)",
        "unmapped": "no ortholog found (mostly LOC ids)"})
    tab.to_csv(RESULTS / "ortholog_mapping_summary.csv", index=False)
    summary["ensembl_release"] = json.loads((D / "ref" / "ensembl_release.json").read_text())["releases"]
    summary["ortholog_methods"] = tab.set_index("method")["n_canine_symbols"].to_dict()


def main() -> None:
    RESULTS.mkdir(exist_ok=True)
    s = {}
    s.update(dog())
    s.update(human())
    orthologs(s)
    (RESULTS / "data_summary.json").write_text(json.dumps(s, indent=1, default=str))
    print(json.dumps(s, indent=1, default=str))


if __name__ == "__main__":
    main()
