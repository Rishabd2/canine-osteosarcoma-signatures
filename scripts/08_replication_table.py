"""Assemble results/replication_table.md (and .csv): paper value vs our value vs agree.

Agreement rule for survival tests, fixed before reading our numbers:
  yes      same call at p < 0.05 (both significant in the same direction, or both not)
  partial  same direction of effect, but one p is below 0.05 and the other is not
  no       opposite direction with at least one p < 0.05
  different partition   our clusters overlap the paper's with ARI < 0.5, so the two
                        p-values test different patient groupings
Direction is the sign of log HR. The paper's direction is the HR from re-testing its
own labels (Supplementary Data 1); our groups are matched to the paper's groups by
majority overlap before comparing signs. Where the paper's labels are not public
(non-metastatic GS splits), the paper's direction is FP better, as its text states.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from cos_sig.paths import RESULTS  # noqa: E402
from cos_sig.survival import bh_adjust  # noqa: E402


def fmt_p(p) -> str:
    if p is None or (isinstance(p, float) and np.isnan(p)):
        return "n/a"
    return f"{p:.2g}" if p >= 0.001 else f"{p:.1e}"


def fmt_bound(v, side: str) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "not estimable"
    if side == "hi" and np.isinf(v):
        return "not estimable"
    if side == "lo" and v == 0:
        return "not estimable"
    return f"{v:.2f}"


def fmt_hr(r) -> str:
    if r is None or np.isnan(r.hr_a_vs_b):
        return "not estimable (a group has no events)" if r is not None else "n/a"
    txt = f"{r.hr_a_vs_b:.2f} [{fmt_bound(r.boot_lo, 'lo')}, {fmt_bound(r.boot_hi, 'hi')}]"
    nf = int(r.get("boot_nonfinite", 0) or 0)
    if nf:
        txt += f" ({nf}/1000 resamples with a zero-event group)"
    return txt


def norm_label(v) -> str:
    s = str(v)
    return s[:-2] if s.endswith(".0") else s


def agree(paper_p, ours_p, same_dir) -> str:
    ps, os_ = paper_p < 0.05, ours_p < 0.05
    if ps == os_:
        return "yes" if (same_dir or not ps) else "no"
    return "partial" if same_dir else "no"


def frozen_section(notes: dict) -> list[str]:
    """Canine-derived signatures frozen before TARGET was read, applied to TARGET."""
    ft = pd.read_csv(RESULTS / "claims_transfer_target.csv")
    ft = ft[ft.signature.str.startswith(("denovo-", "paperlab-")) & (ft.method == "rederived")]
    info = {k["name"]: k for tag in ("isa_denovo", "isa_paperlab") for k in notes.get(tag, {}).get("kept", [])}
    n, hits = len(ft), int((ft.logrank_p < 0.05).sum())
    out = ["", "## Frozen canine signatures in TARGET (leakage-safe test of claim 3)", "",
           f"{len(info)} signatures derived from canine data only, hashed in `results/frozen_signatures.json`",
           "before any human file was read, applied to TARGET all and non-metastatic (k-means, k = 2,",
           f"high vs low by mean signature score), PFS and OS: {n} tests. {hits} reach p < 0.05, against",
           f"{0.05 * n:.0f} expected by chance if none carried signal. `denovo` signatures come from our ssGSEA",
           "clusters; `paperlab` signatures come from the paper's published ssGSEA clusters.", "",
           "| Signature | Genes (closest published, shared) | Cohort | Endpoint | low / high n | p | HR low vs high [boot CI] |",
           "|---|---|---|---|---|---|---|"]
    for _, r in ft.iterrows():
        k = info.get(r.signature, {})
        out.append(f"| {r.signature} | {k.get('genes', '')} ({k.get('closest_published', '')}, "
                   f"{k.get('overlap_with_closest', '')}) | {r.cohort} | {r.endpoint} | {r.n_a} / {r.n_b} | "
                   f"{fmt_p(r.logrank_p)} | {fmt_hr(r)} |")
    return out


def main() -> None:
    paper = pd.read_csv(ROOT / "data" / "paper_reported_values.csv")
    claims = pd.concat([pd.read_csv(RESULTS / f) for f in ["claims_published_signatures.csv",
                                                            "claims_denovo_dog.csv", "claims_transfer_target.csv"]])
    dn_notes = json.loads((RESULTS / "notes_denovo_dog.json").read_text())
    out = []
    for _, pr in paper.iterrows():
        sub = claims[claims.key == pr.key]
        pl = sub[sub.method == "paper_labels"]
        rd = sub[sub.method == "rederived"]
        ac = sub[sub.method == "authors_code"]
        pl = pl.iloc[0] if len(pl) else None
        rd = rd.iloc[0] if len(rd) else None
        ac = ac.iloc[0] if len(ac) else None
        ours = rd if rd is not None else pl
        ours_kind = "re-derived" if rd is not None else "paper labels re-tested"
        if pr.key.startswith("dog_isa") and pr.key[7:8].isdigit():
            # match the paper's ISA-k to the ISA bicluster we re-derived from the paper's ssGSEA clusters
            k = pr.key[7]
            hit = [b for b in dn_notes.get("isa_paperlab", {}).get("kept", [])
                   if b["closest_published"] == f"ISA-{k}" and b["overlap_with_closest"] > 0]
            if hit:
                nm = hit[0]["name"].lower().replace("-", "_")
                cand = claims[claims.key == f"dog_{nm}_dfi"]
                if len(cand):
                    ours = cand.iloc[0]
                    ours_kind = f"re-derived ({hit[0]['name']}, {hit[0]['overlap_with_closest']} shared genes)"
        if ours is None:
            continue
        # paper direction from re-testing the paper's labels; ours matched to the paper's groups
        if pl is not None and np.isfinite(pl.hr_a_vs_b):
            paper_dir = pl.hr_a_vs_b > 1
        elif pl is not None:  # zero-event group: direction from medians (FP never reaches its median)
            paper_dir = True
        else:
            paper_dir = True  # paper text: higher signature enrichment = favorable
        if np.isfinite(ours.hr_a_vs_b):
            ours_dir = ours.hr_a_vs_b > 1
            lab_a = norm_label(ours.get("paper_label_of_a", np.nan))
            if pl is not None and lab_a not in ("nan", "") and lab_a != norm_label(pl.group_a):
                ours_dir = not ours_dir  # our group a corresponds to the paper's group b
            same_dir = ours_dir == paper_dir
        else:
            same_dir = True
        call = agree(pr.paper_p, ours.logrank_p, same_dir)
        if np.isfinite(ours.ari_vs_paper) and ours.ari_vs_paper < 0.5:
            call = "different partition"
        out.append({
            "claim": pr.claim, "key": pr.key, "test": f"{pr.cohort}: {pr.signature}, {pr.endpoint}",
            "paper_groups": pr.n_groups, "paper_p": pr.paper_p, "paper_source": pr.source,
            "paper_labels_p": pl.logrank_p if pl is not None else np.nan,
            "ours_kind": ours_kind, "ours_groups": f"{ours.group_a} {ours.n_a} / {ours.group_b} {ours.n_b}",
            "ours_p": ours.logrank_p, "ours_hr": fmt_hr(ours), "hr_groups": f"{ours.group_a} vs {ours.group_b}",
            "ari_vs_paper": ours.ari_vs_paper,
            "authors_code_p": ac.logrank_p if ac is not None else np.nan,
            "agree": call})
    t = pd.DataFrame(out)
    t["ours_q_bh"] = bh_adjust(t.ours_p.to_numpy())
    t.to_csv(RESULTS / "replication_table.csv", index=False)

    lines = ["# Replication table", "",
             "Generated by `scripts/08_replication_table.py` from `results/claims_*.csv` and",
             "`data/paper_reported_values.csv` (values read from the paper's figures). All p-values are",
             "two-sided log-rank. HR is Cox (Efron ties) for the first group vs the second, with a percentile",
             "95% bootstrap CI over 1000 patient resamples (labels held fixed). `q` is Benjamini-Hochberg",
             "over every p-value in the 'ours' column. 'Paper labels' re-tests the paper's own cluster",
             "assignments from Supplementary Data 1 against the paper's outcomes. ARI compares our",
             "clusters with the paper's.", "",
             "Agreement: **yes** = same call at p < 0.05; **partial** = same direction, call differs;",
             "**no** = opposite direction; **different partition** = our clusters overlap the paper's",
             "with ARI < 0.5. Bootstrap bounds that land on a resample where one group has no events",
             "(HR infinite or zero) are reported as not estimable, with the count of such resamples.", ""]
    hdr = ("| Claim | Test | Paper groups | Paper p | Paper labels re-tested p | Ours | Our p (q) | HR [95% boot CI] |"
           " ARI vs paper | Agree |")
    lines += [hdr, "|" + "---|" * 10]
    for _, r in t.iterrows():
        ari = "" if np.isnan(r.ari_vs_paper) else f"{r.ari_vs_paper:.2f}"
        lines.append(f"| {r.claim} | {r.test} ({r.paper_source}) | {r.paper_groups} | {fmt_p(r.paper_p)} | "
                     f"{fmt_p(r.paper_labels_p)} | {r.ours_kind}: {r.ours_groups} | {fmt_p(r.ours_p)} "
                     f"({fmt_p(r.ours_q_bh)}) | {r.ours_hr} ({r.hr_groups}) | {ari} | {r.agree} |")
    ac = t[t.authors_code_p.notna()]
    if len(ac):
        lines += ["", "## Authors' km_clust.py on the same inputs", "",
                  "| Test | Re-derived p | Authors' code p |", "|---|---|---|"]
        for _, r in ac.iterrows():
            lines.append(f"| {r.test} | {fmt_p(r.ours_p)} | {fmt_p(r.authors_code_p)} |")

    # claim 2: de novo biclusters (no one-to-one paper counterpart)
    dn = pd.read_csv(RESULTS / "claims_denovo_dog.csv")
    dn = dn[dn.key.str.contains("denovo_isa|paperlab_isa")]
    notes = json.loads((RESULTS / "notes_denovo_dog.json").read_text())
    lines += ["", "## Claim 2: biclusters re-derived in DOG2", "",
              "Paper: 4 ISA biclusters; only ISA-4 separated DFI (p = 0.026). Ours: ISA consensus biclusters from",
              "25 seeded runs (kept if present in more than 75% of runs), from our ssGSEA clusters (`denovo`) and from",
              "the paper's published ssGSEA clusters (`paperlab`).", "",
              "| Bicluster | Genes | Closest published (shared genes) | Endpoint | ISA / rest n | p | HR ISA vs rest [boot CI] |",
              "|---|---|---|---|---|---|---|"]
    info = {k["name"]: k for tag in ("isa_denovo", "isa_paperlab") for k in notes.get(tag, {}).get("kept", [])}
    for _, r in dn.iterrows():
        nm = r.signature.split(" ")[0]
        k = info.get(nm, {})
        lines.append(f"| {nm} | {k.get('genes', '')} | {k.get('closest_published', '')} "
                     f"({k.get('overlap_with_closest', '')}) | {r.endpoint} | {r.n_a} / {r.n_b} | {fmt_p(r.logrank_p)} | "
                     f"{fmt_hr(r)} |")

    # claim 4: deconvolution and GSEA
    dec = pd.read_csv(RESULTS / "deconvolution_fp_vs_pp.csv")
    dec = dec[dec.paper_row.notna() & (dec.paper_row != "")]
    paper_dir = {"T cells CD8": "FP>PP", "Macrophages M1": "FP>PP", "Macrophages M2": "FP>PP",
                 "Macrophages (all)": "n/a (paper: M0 PP>FP, M1 and M2 FP>PP)"}
    lines += ["", "## Claim 4: immune deconvolution, GS-1 FP vs PP", "",
              "Paper: CIBERSORTx LM22 (Table 1). Ours: MCP-counter and TIL10 constrained least squares",
              "(open substitutes; see docs/DEVIATIONS.md). Mann-Whitney U, BH within dataset and method.", "",
              "| Species | Method | Cell type | Paper direction | Our direction | Our q | Agree |", "|---|---|---|---|---|---|---|"]
    for _, r in dec.iterrows():
        pdx = paper_dir.get(r.paper_row, "")
        ag = "n/a" if pdx.startswith("n/a") else ("yes" if (r.direction == pdx and r.q_bh < 0.05) else
                                                  ("partial" if r.direction == pdx else "no"))
        lines.append(f"| {r.species} | {r.method} | {r.cell} | {pdx} | {r.direction} | {fmt_p(r.q_bh)} | {ag} |")
    gs = json.loads((RESULTS / "notes_immune_pathways.json").read_text())["gsea"]
    lines += ["", "## Claim 4: Hallmark GSEA, FP vs PP (FDR < 0.05, enriched in FP)", "",
              "| Cohort | Paper sets | Our sets | Shared |", "|---|---|---|---|"]
    for c in ["DOG2", "TARGET all", "TARGET non-metastatic", "TARGET metastatic"]:
        g = gs[c]
        lines.append(f"| {c} | {g['paper_fp_sig']} | {g['ours_fp_sig']} | {g['overlap']} |")
    lines.append(f"\nDOG2 FP-enriched sets that are also FP-enriched in TARGET all: "
                 f"{gs['dog_fp_sets_also_sig_in_TARGET_all']} of {gs['DOG2']['ours_fp_sig']} "
                 f"(paper: 10 of 10). Sets enriched in PP at FDR < 0.05: "
                 + ", ".join(f"{c} {gs[c]['ours_pp_sig']}" for c in ["DOG2", "TARGET all"]) + " (paper: none).")
    lines += frozen_section(notes)
    (RESULTS / "replication_table.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
