# R 4.4 with Bioconductor edgeR, limma, sva, and CRAN survival, isa2.
# Versions used for the committed results: R 4.4.3, edgeR 4.4.0, limma 3.62.1, sva 3.54.0, isa2 0.3.6.
# One way to get them: conda create -p <env> -c conda-forge -c bioconda r-base=4.4 bioconductor-edger
#   bioconductor-limma bioconductor-sva r-survival ; then run this file to add isa2 (needs a C compiler).
if (!requireNamespace("isa2", quietly = TRUE)) install.packages("isa2", repos = "https://cloud.r-project.org")
for (p in c("edgeR", "limma", "sva", "isa2")) cat(p, as.character(packageVersion(p)), "\n")
