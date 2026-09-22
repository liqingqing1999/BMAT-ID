# External inputs

This repository contains no sequencing data. Two public datasets are needed, and
both are downloaded by hand before running the pipeline.

## 0. Project root

Every script locates the project root from its own position (the parent of
`scripts/`) and needs no editing. If you keep the scripts somewhere else, point
the environment variable **`BMAT_ID_DIR`** at the root:

```bash
export BMAT_ID_DIR=/path/to/BMAT-ID
```

All 41 scripts read this variable, so it is the only one needed for Stages 1–4
and 6–7. The two variables below are needed only by the reference-atlas stage.

## 1. Bulk RNA-seq (the data analysed here)

**GSE291355** — human bone marrow adipocytes, 24 libraries (4 donors × 3
anatomical sites × 2 culture states).

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291355>
- Download the deposited counts matrix and place it at:

```
data/raw/GSE291355_counts.tsv.gz
```

The file must be a gzipped, tab-separated gene-by-sample matrix:

| property | expected |
|---|---|
| columns | 25: the Ensembl gene id, then one column per library |
| header row | the 24 library ids, in the form `H410_epi_primary_adip` — the gene-id column is left unnamed, which is what `read.table(..., row.names = 1)` expects |
| rows | 46,206 genes (plus the header) |
| gene ids | Ensembl, version suffixes allowed (the scripts strip `.<n>`) |
| file size | ≈ 1.20 MB gzipped |
| file md5 | `d3d585fa5cd75bf2cdb211df22651cd2` |

The library ids are decoded, never hard-coded: the scripts parse the donor
(`H410`), the site (`epi` / `meta` / `scat`) and the state (`primary` / `diff`)
out of the column names, so the column order does not matter but the naming
convention does. No preprocessing is applied beyond a CPM filter inside the
scripts: every script starts from these raw counts.

*Confirmed on the copy used for the manuscript; if your download differs, check
the md5 first — a different counts matrix will not reproduce the reported
numbers.*

## 2. Single-cell reference for the composition signatures

**GSE169396** — human femoral-head single-cell atlas (4 donors, 18,205 cells,
26 clusters). Only its non-adipocyte compartments are used, and only to build
the composition signatures of Methods 2.6.

- GEO series: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169396>

Two files are needed and both are referenced through environment variables so
that nothing has to be edited:

| environment variable | what to point it at | used by |
|---|---|---|
| `BMAT_REF_RDS` | a Seurat object (`.rds`) holding the full atlas, with a `seurat_clusters` column and a `group` column in its metadata | `30_*` only |
| `BMAT_REF_FEATURES` | the features/genes file shipped with the GSE169396 supplementary archive (`GSM5201883_S1_features.tsv.gz`), used only to bridge Ensembl ids to gene symbols offline | `31_*` only |

`32_*` needs neither: it reads the score table written by `31_*`.

Build the atlas with Seurat from the GEO supplementary matrices, e.g.

```r
# sketch -- see the GSE169396 series for the exact files
library(Seurat)
counts <- Read10X("GSE169396/")
obj    <- CreateSeuratObject(counts)
obj    <- NormalizeData(obj)
obj    <- FindVariableFeatures(obj)
obj    <- ScaleData(obj)
obj    <- RunPCA(obj)
obj    <- FindNeighbors(obj, dims = 1:20)
obj    <- FindClusters(obj, resolution = 0.6)
obj$group <- ...      # the donor / clinical group labels of the series
saveRDS(obj, "gse169396_all.rds")
```

then

```bash
export BMAT_REF_RDS=/path/to/gse169396_all.rds
export BMAT_REF_FEATURES=/path/to/GSM5201883_S1_features.tsv.gz
```

Reference processing needs `Seurat >= 5` (the scripts call `JoinLayers()`).
