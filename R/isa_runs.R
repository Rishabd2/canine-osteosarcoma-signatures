# ISA biclustering, following ISA_Cluster.R and GET_ISA_SEEDs.py from Mannheimer et al.,
# Zenodo 10.5281/zenodo.8125077 (CC BY 4.0). The original calls the Bioconductor 'eisa'
# wrappers (ISANormalize, ISAIterate, ISAUnique, ISAFilterRobust), which are thin wrappers over
# the CRAN 'isa2' functions used here (isa.normalize, isa.iterate, isa.unique, isa.filter.robust).
# eisa is not installable for current R/Bioconductor on arm64.
# Seeds: hclust(method = "ward.D") on genes (cut to 10) and samples (cut to 7), then 3 random
# genes / samples per hierarchical cluster, as in GET_ISA_SEEDs.py. Repeated n_runs times.
# Usage: Rscript isa_runs.R expr.csv out_prefix n_runs seed
#   expr.csv: samples x genes (log, normalized), as in the original RNA_SEQ_DATA.csv
args <- commandArgs(trailingOnly = TRUE)
suppressPackageStartupMessages(library(isa2))
X <- as.matrix(read.csv(args[1], row.names = 1, check.names = FALSE))
X <- t(X)                     # genes x samples
n_runs <- as.integer(args[3])
set.seed(as.integer(args[4]))
gcl <- cutree(hclust(dist(X), method = "ward.D"), k = 10)
scl <- cutree(hclust(dist(t(X)), method = "ward.D"), k = 7)
nm <- isa.normalize(X)
rows <- list()
for (r in seq_len(n_runs)) {
  fs <- matrix(0, nrow(X), 10); ss <- matrix(0, ncol(X), 7)
  for (i in 1:10) { idx <- which(gcl == i); fs[if (length(idx) > 3) sample(idx, 3) else idx, i] <- 1 }
  for (i in 1:7)  { idx <- which(scl == i); ss[if (length(idx) > 3) sample(idx, 3) else idx, i] <- 1 }
  bc <- isa.iterate(nm, row.seeds = fs, col.seeds = ss, thr.row = 1.5, thr.col = 0.7)
  bc <- isa.unique(nm, bc)
  bc <- isa.filter.robust(X, nm, bc)
  nb <- ncol(bc$rows)
  if (nb > 0) for (b in seq_len(nb)) {
    g <- rownames(X)[bc$rows[, b] != 0]
    s <- colnames(X)[bc$columns[, b] != 0]
    rows[[length(rows) + 1]] <- data.frame(run = r, bicluster = b,
      genes = paste(g, collapse = ";"), samples = paste(s, collapse = ";"))
  }
  cat("run", r, ":", nb, "robust biclusters\n")
}
out <- do.call(rbind, rows)
write.csv(out, file = paste0(args[2], "_isa_runs.csv"), row.names = FALSE)
