## Export a reference signature from the local human bone-marrow scRNA-seq atlas
## (GSE169396, femoral head, 4 donors, 18,205 cells) for bulk deconvolution of GSE291355.
# ---- portable project root -------------------------------------------------
# Nothing to edit: the root is the parent of this script's own directory.
# Override with the environment variable BMAT_ID_DIR if you keep the scripts
# somewhere else.  See README.md, section "Running the pipeline".
base_dir <- local({
  env <- Sys.getenv("BMAT_ID_DIR", unset = "")
  if (nzchar(env)) return(normalizePath(env, mustWork = FALSE))
  a <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(a)) return(normalizePath(file.path(dirname(a[1]), ".."), mustWork = FALSE))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(file.path(dirname(of), ".."), mustWork = FALSE))
  }
  normalizePath(".", mustWork = FALSE)
})
for (.d in c("logs", "results/tables", "results/figures"))
  dir.create(file.path(base_dir, .d), recursive = TRUE, showWarnings = FALSE)
BASE <- base_dir
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

outdir <- file.path(BASE, "results", "tables")
logf   <- file.path(BASE, "logs", "ref_export_out.txt")
con <- file(logf, open = "wt", encoding = "UTF-8")
say <- function(...) { writeLines(paste0(...), con); flush(con) }

p <- Sys.getenv("BMAT_REF_RDS", unset = "")
if (!nzchar(p) || !file.exists(p))
  stop("Set BMAT_REF_RDS to the GSE169396 Seurat object (.rds); see README.md, ", 
       "section \"External inputs\".")
say("== reading GSE169396 atlas from ", p)
say("== reading ", p)
o <- readRDS(p)
say("dim ", nrow(o), " x ", ncol(o))
say("assays: ", paste(Assays(o), collapse = ","), " | default: ", DefaultAssay(o))
say("rownames head: ", paste(head(rownames(o), 5), collapse = ", "))
say("meta: ", paste(colnames(o@meta.data), collapse = " | "))

DefaultAssay(o) <- "RNA"
say("RNA layers: ", paste(Layers(o[["RNA"]]), collapse = ","))
if (length(Layers(o[["RNA"]])) > 1) {
  o[["RNA"]] <- tryCatch(JoinLayers(o[["RNA"]]),
                         error = function(e) { say("JoinLayers err: ", conditionMessage(e)); o[["RNA"]] })
  say("after JoinLayers: ", paste(Layers(o[["RNA"]]), collapse = ","))
}
m <- NULL
for (ly in c("counts", "data", "scale.data")) {
  m <- tryCatch(GetAssayData(o, assay = "RNA", layer = ly),
                error = function(e) { say("  layer '", ly, "' err: ", conditionMessage(e)); NULL })
  if (!is.null(m) && prod(dim(m)) > 0) { say("  using RNA layer: ", ly); break }
  m <- NULL
}
if (is.null(m)) stop("no RNA matrix available")
say("counts matrix: ", class(m)[1], " ", nrow(m), " x ", ncol(m))
say("counts total: ", sum(m))

ct  <- as.character(o$seurat_clusters)
grp <- as.character(o$group)
cl  <- sort(unique(ct))
say("n clusters: ", length(cl))

## pseudo-bulk per cluster (sum of raw counts)
pb <- sapply(cl, function(k) Matrix::rowSums(m[, ct == k, drop = FALSE]))
colnames(pb) <- cl
say("pseudo-bulk dim: ", nrow(pb), " x ", ncol(pb))
write.csv(data.frame(ensg = rownames(pb), pb, check.names = FALSE),
          file.path(outdir, "BMATID_ref_pseudobulk_clusters.csv"), row.names = FALSE)

## cluster x donor-group composition
tb <- table(cluster = ct, group = grp)
write.csv(as.data.frame.matrix(tb), file.path(outdir, "BMATID_ref_cluster_by_group.csv"))
say("\n== cluster x group ==")
print(tb)

## marker scoring on pseudo-bulk (log2 CPM, z-scored across clusters)
cpm <- t(t(pb) / colSums(pb) * 1e6)
lg  <- log2(cpm + 1)
z   <- t(scale(t(lg)))          # z across clusters per gene
z[is.na(z)] <- 0
markers <- list(
  haematopoietic = c("PTPRC","CD14","LYZ","CD68","CSF1R","SPI1"),
  T_cell         = c("CD3D","CD3E","CD2","IL7R","CCL5"),
  B_cell         = c("MS4A1","CD79A","CD19","BANK1","TCL1A"),
  plasma_cell    = c("MZB1","JCHAIN","IGHG1","SDC1","DERL3","IGKC"),
  osteoblast     = c("RUNX2","SP7","ALPL","IBSP","COL1A1","BGLAP","SPP1"),
  osteocyte      = c("DMP1","PHEX","SOST","MEPE","ENPP1"),
  MSC_stroma     = c("LEPR","PDGFRB","NT5E","ENG","CXCL12","ADIPOR1"),
  adipocyte      = c("PLIN1","PLIN4","ADIPOQ","FABP4","LPL","CIDEC","CFD"),
  endothelial    = c("PECAM1","VWF","CDH5","CLDN5"),
  erythroid      = c("HBB","HBA1","ALAS2","AHSP"),
  chondrocyte    = c("ACAN","COL2A1","SOX9","COMP"),
  mural_muscle   = c("ACTA2","TAGLN","MYH11","NOTCH3")
)
present <- rownames(z)
sc <- sapply(markers, function(g) {
  g <- intersect(g, present)
  if (length(g) == 0) return(rep(NA_real_, ncol(z)))
  colMeans(z[g, , drop = FALSE])
})
rownames(sc) <- colnames(z)
say("\n== marker score per cluster (mean z; NA = no marker present) ==")
print(round(sc, 2))
say("\n== marker coverage ==")
for (k in names(markers)) {
  g <- markers[[k]]; say("  ", k, ": ", length(intersect(g, present)), "/", length(g),
                         " present -> ", paste(intersect(g, present), collapse = ","))
}
write.csv(data.frame(cluster = rownames(sc), sc, check.names = FALSE),
          file.path(outdir, "BMATID_ref_cluster_markerscores.csv"), row.names = FALSE)

## also export: mean expression of the 10 GSE291355 composition-axis markers by cluster
axis_markers <- c("PTPRC","CD14","CD68","LYZ","CD3E","MS4A1","MZB1","JCHAIN","CSF1R","CD19",
                  "CD79A","COL1A1","RUNX2","SP7","ALPL","IBSP","SPP1","MEPE","DMP1","PHEX",
                  "BGLAP","TMEM119","ENPP1","ADIPOQ","PLIN1","PLIN4","FABP4","GPAM","CFD","LEP","LPL","PLIN2","CIDEC")
am <- intersect(axis_markers, present)
say("\n== axis markers present in reference: ", length(am), "/", length(axis_markers))
say("   ", paste(am, collapse = ","))
if (length(am) > 0) {
  sub <- as.data.frame(round(lg[am, , drop = FALSE], 2))
  sub <- data.frame(gene = rownames(sub), sub, check.names = FALSE)
  write.csv(sub, file.path(outdir, "BMATID_ref_axismarkers_log2cpm.csv"), row.names = FALSE)
}

say("\nDONE")
close(con)
cat("done\n")
