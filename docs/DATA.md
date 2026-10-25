# Data sources

Nothing below is committed except small derived tables in `results/` and the two hand-made files
in `data/` (`paper_reported_values.csv`, values read from the paper's figures, and
`checksums.sha256`). `python scripts/fetch_data.py` downloads everything into `$COS_DATA_DIR`
(default `data/raw/`), checks upstream md5 sums where the source publishes them (Zenodo, GDC), and
checks sha256 against `data/checksums.sha256`.

| Dataset | What we use | Source | License / terms | Citation |
|---|---|---|---|---|
| DOG2 canine osteosarcoma RNA-seq | raw counts, 186 primary tumours; sample titles (trial) | GEO GSE238110, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238110 | NCBI GEO public data (no use restrictions stated by NCBI) | Mannheimer et al. 2023 |
| Paper supplementary data | Supplementary Data 1 (per-patient outcomes and cluster labels), Supplementary Data 2 (signature gene lists), Supplementary Information PDF (Supp. Table 3 and figures) | https://doi.org/10.1038/s42003-023-05208-z | CC BY 4.0 | Mannheimer et al. 2023 |
| Authors' code | 9 scripts (Python, R) | Zenodo, https://doi.org/10.5281/zenodo.8125077 | CC BY 4.0 | Mannheimer and LeBlanc 2023 |
| TARGET-OS RNA-seq | 88 STAR-Counts files (open access), 79 used | NCI GDC, https://portal.gdc.cancer.gov/projects/TARGET-OS | NCI GDC open-access data, subject to the GDC Data Use Agreement and TARGET publication policy | TARGET Osteosarcoma project, NCI (phs000468) |
| TARGET-OS clinical | GDC case records (vital status, metastasis at diagnosis) for cross-checks | GDC API `/cases` | as above | as above |
| Ensembl Compara orthologs | dog (ROS_Cfam_1.0) to human gene pairs | Ensembl BioMart, release 116 | Ensembl data are freely available without restriction | Harrison et al., Ensembl 2024, Nucleic Acids Res. |
| CanFam4 annotation | GTF for gene lengths | NCBI GCF_011100685.1 (UU_Cfam_GSD_1.0) | NCBI public data | Wang et al. 2021, Commun. Biol. 4:185 |
| MSigDB Hallmark and Reactome 7.5.1 | gene sets (human symbols) | https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/ | MSigDB terms (CC BY 4.0 for these collections); Reactome CC0 | Liberzon et al. 2015, Cell Syst. 1:417 |
| quanTIseq TIL10 | signature matrix | github.com/icbi-lab/quanTIseq @ 039eb16 | BSD-style (see repository) | Finotello et al. 2019, Genome Med. 11:34 |
| MCP-counter markers | marker gene table | github.com/ebecht/MCPcounter @ b6eac73 | GPL-3 (marker table used as data, not redistributed) | Becht et al. 2016, Genome Biol. 17:218 |

## Disk layout used for the committed run

The internal disk of the machine that produced `results/` had little free space, so raw data,
intermediate matrices, the Python venv and the R environment lived on an external volume:

```
COS_DATA_DIR=/Volumes/T9/claudia/github-data/canine-osteosarcoma-signatures/data
COS_WORK_DIR=/Volumes/T9/claudia/github-data/canine-osteosarcoma-signatures/work
```

with `data/raw` and `data/work` symlinked there (both are gitignored). Any location works.
