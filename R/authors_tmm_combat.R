# Adapted from TMM_Normalization.R in Mannheimer et al., Zenodo 10.5281/zenodo.8125077 (CC BY 4.0).
# Changes: file paths come from the command line. The analysis lines are unchanged,
# including calcNormFactors() whose result is not assigned in the original script,
# so the TMM factors are computed but not applied (cpm() then uses library sizes only).
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages({library(sva); library(edgeR)})
Xdata <- read.csv(file = args[1], header = TRUE, row.names = 1, check.names = FALSE)
mdata <- read.csv(file = args[2], header = TRUE, row.names = 1, check.names = FALSE)
batch_labels <- as.factor(mdata$Batch.ID)
Xcount <- DGEList(Xdata)
keep <- filterByExpr(Xcount, min.count = 10, min.prop = 0.5, large.n = 5)
Xv <- Xcount[keep, , keep.lib.sizes = FALSE]
calcNormFactors(Xv, method = 'TMM')  # as in the original: result discarded
Xcpm <- cpm(Xv, normalized.lib.sizes = TRUE, log = TRUE)
Xcombat <- ComBat(Xcpm, batch = batch_labels)
write.csv(Xcombat, file = args[3])
cat("authors_tmm_combat: kept", sum(keep), "of", length(keep), "probes\n")
