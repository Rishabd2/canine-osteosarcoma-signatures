# Deviations from the paper

Every place where this replication does something other than what Mannheimer et al. (2023)
describe, and why. Items marked **(forced)** were not a choice: the original input or tool is not
public, is proprietary, or no longer installs.

## Inputs

1. **Canine outcomes come from the paper's Supplementary Data 1 (forced).** GEO GSE238110 has raw
   counts and sample titles but no clinical fields. The source-data workbook
   (`42003_2023_5208_MOESM4_ESM.xlsx`, sheet "Figure 2") has DFI, DFI censoring, overall survival
   and survival censoring for all 186 dogs, keyed by COTC patient number, which matches the GEO
   column names (for example `10_0512__Tumor` is dog 512). One GEO column is named
   `30_0814_RNAla_Ter__tumor`; it is dog 814 and is included. Four primary-tumour columns in GEO
   (902, 1701, 1702, 1706) are not in the 186 and are dropped, as the paper did.
2. **TARGET-OS outcomes come from the same workbook (sheet "Figure 3").** PFS and OS for the 79
   analysed cases. We cross-checked against GDC case records: OS event status matches
   `vital_status` for 79/79 cases, and the paper's metastatic split (sheets "Figure S4" and
   "Figure 9") matches GDC `metastasis_at_diagnosis` for 79/79. We did not recompute the paper's
   RNA QC exclusion (6 of 85 specimens with exon/coding fraction below 0.5 and intergenic rate
   at least 0.3); we used its list of 79 cases.
3. **TARGET expression is the current GDC STAR-Counts release (forced).** The GDC re-processed
   TARGET-OS with the STAR two-pass / GENCODE v36 pipeline; the paper used an earlier release.
   We use `tpm_unstranded` and `unstranded` counts for protein-coding genes (19,986 symbols).
4. **Canine TPM is computed from counts and CanFam4 gene lengths.** GEO has counts only. We take
   union-of-exon lengths per gene from the NCBI UU_Cfam_GSD_1.0 annotation (the reference the
   authors aligned to). RSEM's effective-length TPM, which the authors used, is not public.
5. **Signature lists.** GS-1 and GS-2 come from Supplementary Data 2 ("Gene Signature"): 25 GS-1,
   156 GS-2 and 2 "DOG Only" genes. The paper says 27 GS-1 genes were present in DOG2 and that two
   canine-specific genes were absent from TARGET, and Fig. 2a lists DLA-DQA1 and DLA-DQB1 among
   the GS-1 rows, so canine GS-1 = 25 + the 2 DLA genes, human GS-1 = 25. ISA signatures are
   transcribed from Supplementary Table 3. The main text gives ISA-3 39 genes and ISA-4 15 genes;
   the table (and Fig. 7c,d) give 15 and 39. We use the table.

## Methods

6. **Authors' `km_clust.py` needed four patches to run.** `np.argmax(sil_score)[0]` raises
   `IndexError` in every NumPy version, so the published script cannot run as written; `np.in1d`
   no longer exists in NumPy 2.4; and scikit-learn changed the `KMeans` default from 10 starts to
   one, so we restore `n_init=10`. The signature file name is pointed at our list. The patched
   source is built at run time from the downloaded file (`src/cos_sig/authors.py`); no copy of the
   authors' code is committed.
7. **Authors' `TMM_Normalization.R` is run as written, including a no-op line.** It calls
   `calcNormFactors(Xv, method='TMM')` without assigning the result, so the TMM factors are never
   applied and the log-CPM uses raw library sizes. We keep that behaviour for the ssGSEA input
   (`R/authors_tmm_combat.R`) because it is what the code does; ComBat batch correction on trial
   (COTC021 vs COTC022) follows.
8. **Trial, arm and batch coincide.** In GEO, COTC021 = the rapamycin arm (93 dogs) and COTC022 =
   standard of care (93 dogs), so the "batch" the authors corrected for is also the treatment arm.
   The paper reports no arm effect on outcome (Supp. Fig. 11), so this matters little here, but
   ComBat on this factor would also remove any true treatment-associated expression.
9. **ssGSEA is our implementation, not GenePattern 10.1.0.** `cos_sig.scoring.ssgsea` follows the
   Barbie et al. / GSVA formulation with weight 0.75 (GenePattern's default exponent), on Hallmark
   7.5.1 after mapping canine symbols to human orthologs. GenePattern's own rank handling and
   score normalisation differ in detail. The k-means step z-scores the ssGSEA scores; the authors'
   script would have taken log2(x + 2) of them, which is undefined for negative scores, so we read
   "same format" as standardise only.
10. **The ssGSEA FP/PP naming rule is ours.** The paper named clusters by their survival. To keep
    outcome out of the labelling, we call FP the cluster with the higher mean ssGSEA z-score over
    three immune Hallmarks (allograft rejection, interferon-gamma response, inflammatory
    response). Log-rank p-values do not depend on the naming; HR direction does. For GS-1 / GS-2,
    FP is the cluster with the higher mean signature z-score, which is how the paper describes
    them ("increased enrichment of GS-1 was associated with improved DFI").
11. **DEG threshold before ISA.** The paper reports 259 DEGs between its ssGSEA clusters; the
    limma code keeps BH adj.P < 0.05 only. With 186 dogs that cutoff gives 7,646 probes for our
    clusters and 7,325 for the paper's clusters, so we apply the Methods' stated DEG threshold for
    pathway input (adj.P < 0.05 and fold change > 3) before ISA. That leaves 109 probes (94 genes)
    for our clusters and 234 probes (190 genes) for the paper's clusters (paper: 259).
12. **ISA uses CRAN `isa2` instead of Bioconductor `eisa` (forced).** `eisa` is deprecated and does
    not install on current R/Bioconductor for arm64; its `ISANormalize`, `ISAIterate`,
    `ISAUnique` and `ISAFilterRobust` are thin wrappers over `isa2`'s `isa.normalize`,
    `isa.iterate`, `isa.unique` and isa2's robustness filter, which we call with the authors'
    thresholds (1.5 genes, 0.7 samples). The seeding (Ward hierarchical clustering into 10 gene and
    7 sample groups, 3 random members each) follows `GET_ISA_SEEDs.py`. The paper does not say how
    biclusters were matched across its 25 runs; we join biclusters with gene Jaccard of at least
    0.5 and keep groups found in more than 75% of runs.
13. **Pathway analysis: Reactome over-representation instead of IPA (forced).** IPA is proprietary.
    Shared DEGs go through a hypergeometric test on MSigDB Reactome 7.5.1 with BH FDR. The
    Hallmark GSEA uses gseapy (signal-to-noise, 1000 phenotype permutations) in place of Broad GSEA
    4.1.0 desktop.
14. **Deconvolution: MCP-counter and TIL10 least squares instead of CIBERSORTx (forced).**
    CIBERSORTx needs registration and LM22 is distributed under its terms. We use MCP-counter
    marker means and a constrained least-squares fit on the quanTIseq TIL10 matrix (fractions
    non-negative, summing to at most 1). The TIL10 fit omits quanTIseq's mRNA-content rescaling and
    noisy-gene removal, so its numbers are relative fractions, not quanTIseq output. Neither
    method has an "M0 macrophage" state, so the paper's strongest cell-type result (M0 higher in
    PP) has no direct counterpart. Canine TPM is mapped to human symbols with Ensembl orthologs
    before deconvolution, which is as species-naive as the paper's use of LM22 on dog tissue.
15. **Statistics in Python.** Log-rank and Kaplan-Meier use lifelines rather than survminer. For
    the paper's own labels, our log-rank p-values equal the published ones at the printed
    precision in 23 of 24 re-tested comparisons; the exception is DOG2 GS-2 DFI (0.026 here,
    0.028 printed).
16. **Fig. 7 endpoint.** Fig. 7's y-axis and legend say overall survival; the text says DFI.
    Re-testing the paper's own ISA labels reproduces the printed p-values on DFI (0.81, 0.32, 0.16,
    0.026) and not on OS (0.93, 0.24, 0.48, 0.014), so we treat Fig. 7 as DFI.

## Added checks (not in the paper)

17. Bootstrap CIs on every HR; BH across all replication p-values; k = 2..6 sensitivity;
    seed stability of k-means; whole-pipeline bootstrap (re-cluster each resample);
    random-signature null in TARGET; a frozen, hashed copy of the canine-derived signatures that
    the human-transfer script verifies before running; and a static check that the canine
    derivation script never opens human data.
