# Adapted from LIMMA_DOG2.R and TARGETS_Limma.R in Mannheimer et al., Zenodo 10.5281/zenodo.8125077 (CC BY 4.0).
# Usage: Rscript limma_voom.R counts.csv labels.csv out.csv [batch|nobatch]
# labels.csv: row names = sample ids, columns Cluster (0/1) and, for dogs, Batch.ID.
# Changes: paths from the command line; the full table is written (not only adj.P < 0.05)
# so callers can apply the paper's thresholds; filterByExpr uses min.count (the DOG2 script
# passed "min_count", which filterByExpr ignores; the default is also 10).
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(edgeR))
Xdata1 <- read.csv(args[1], header = TRUE, row.names = 1, check.names = FALSE)
lab <- read.csv(args[2], header = TRUE, row.names = 1, check.names = FALSE)
lab <- lab[colnames(Xdata1), , drop = FALSE]
cluster <- as.factor(lab$Cluster)
Xd <- DGEList(Xdata1)
keep <- filterByExpr(Xd, min.count = 10, min.prop = 0.5, large.n = 5)
X0 <- Xd[keep, , keep.lib.sizes = FALSE]
X0 <- calcNormFactors(X0, method = 'TMM')
if (length(args) >= 4 && args[4] == "batch") {
  batch <- as.factor(lab$Batch.ID)
  mm <- model.matrix(~batch + cluster)
} else {
  mm <- model.matrix(~cluster)
}
y <- voom(X0, mm, plot = FALSE)
fit <- eBayes(lmFit(y, mm))
tt <- topTable(fit, coef = 'cluster1', sort.by = 'P', number = Inf, adjust = "BH")
write.csv(tt, file = args[3])
cat("limma_voom:", sum(tt$adj.P.Val < 0.05), "genes with BH adj.P < 0.05 of", nrow(tt), "\n")
