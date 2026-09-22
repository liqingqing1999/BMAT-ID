## Composition-axis attribution: refit the LRT after removing genes that co-vary
## with an axis, keeping the rest (keep = setdiff(all, to_remove)).
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

out <- file.path(BASE, "logs", "deseq2_axis_out.txt")
con <- file(out, "w"); w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }
suppressMessages(library(DESeq2))

dd <- read.delim(gzfile(file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")),
                 header = TRUE, row.names = 1, check.names = FALSE)
rownames(dd) <- sub("\\..*$", "", rownames(dd)); dd <- dd[!duplicated(rownames(dd)), ]
cn <- colnames(dd)
donor <- sub("^H([0-9]+)_.*$", "\\1", cn)
location <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", cn)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", cn)
CONTRASTS <- list(c("meta", "epi"), c("scat", "epi"), c("scat", "meta"))

quantify <- function(idx, keep = NULL, do_lrt = FALSE) {
  cnt <- round(as.matrix(dd[, idx]))
  if (!is.null(keep)) { k <- intersect(rownames(cnt), keep); cnt <- cnt[k, , drop = FALSE] }
  cnt <- cnt[rowSums(cnt) >= 10, , drop = FALSE]
  m <- data.frame(donor = factor(donor[idx]),
                  location = factor(location[idx], levels = c("epi", "meta", "scat")))
  dds <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  dds <- DESeq(dds, quiet = TRUE)
  res <- list()
  for (ct in CONTRASTS) {
    r <- as.data.frame(lfcShrink(dds, contrast = c("location", ct[1], ct[2]),
                                 type = "normal", quiet = TRUE))
    r$gene <- rownames(r); res[[sprintf("%s_vs_%s", ct[1], ct[2])]] <- r
  }
  allg <- unique(unlist(lapply(res, function(x) x$gene)))
  nFDR <- 0; nLFC1 <- 0; nLFC05 <- 0
  for (g in allg) {
    mf <- NA_real_; ml <- 0
    for (nm in names(res)) {
      x <- res[[nm]]; j <- which(x$gene == g); if (!length(j)) next
      p <- x$padj[j]; if (is.na(p)) next
      mf <- if (is.na(mf)) p else min(mf, p)
      if (p < 0.05) ml <- max(ml, abs(x$log2FoldChange[j]))
    }
    if (!is.na(mf) && mf < 0.05) { nFDR <- nFDR + 1; if (ml > 1) nLFC1 <- nLFC1 + 1; if (ml > 0.5) nLFC05 <- nLFC05 + 1 }
  }
  LRT05 <- NA_integer_; ntest <- NA_integer_
  if (do_lrt) {
    d2 <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
    d2 <- DESeq(d2, test = "LRT", reduced = ~ donor, quiet = TRUE)
    r2 <- results(d2); LRT05 <- sum(!is.na(r2$padj) & r2$padj < 0.05); ntest <- sum(!is.na(r2$padj))
  }
  list(ngene = nrow(cnt), ntest = ntest, LRT05 = LRT05,
       wald05 = nFDR, wald_lfc1 = nLFC1, wald_lfc05 = nLFC05)
}

rc <- readRDS(file.path(BASE, "results/tables/BMATID_axis_correlations.rds"))
rownames(rc) <- sub("\\..*$", "", rownames(rc))
H <- names(which(abs(rc[, "Hemato"]) > 0.7))
B <- names(which(abs(rc[, "Bone"])   > 0.7))
A <- names(which(abs(rc[, "Adipo"])  > 0.7))
SCHEMES <- list(
  "基准：全基因"        = NULL,
  "剔 Hemato"           = H,
  "剔 Bone"             = B,
  "剔 Hemato+Bone"      = union(H, B),
  "剔 Hemato+Bone+Adipo" = union(union(H, B), A))

w("==================================================================")
w("  组成轴归因 v2（DESeq2 + lfcShrink normal + 效应量阈值）  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("  keep = setdiff(全部基因, 该轴 |r|>0.7 的基因)；Wald 三对比并集")
w("==================================================================")
tab <- data.frame()
for (nm in names(SCHEMES)) {
  rm_g <- SCHEMES[[nm]]
  keep <- if (is.null(rm_g)) NULL else setdiff(rownames(dd), rm_g)
  rem <- if (is.null(rm_g)) 0 else length(rm_g)
  w(""); w(sprintf("### %s   （剔除 %d 个基因）", nm, rem))
  P <- quantify(which(state == "primary"), keep, do_lrt = TRUE)
  D <- quantify(which(state == "diff"), keep, do_lrt = FALSE)
  w(sprintf("  体内 primary  tested=%d  LRT FDR<.05=%d  并集 FDR<.05=%d  &|LFC|>.5=%d  &|LFC|>1=%d",
            P$ntest, P$LRT05, P$wald05, P$wald_lfc05, P$wald_lfc1))
  w(sprintf("  体外 diff     并集 FDR<.05=%d  &|LFC|>.5=%d  &|LFC|>1=%d", D$wald05, D$wald_lfc05, D$wald_lfc1))
  tab <- rbind(tab, data.frame(scheme = nm, removed = rem,
                               prim_ntest = P$ntest, prim_LRT05 = P$LRT05,
                               prim_wald05 = P$wald05, prim_lfc1 = P$wald_lfc1,
                               diff_wald05 = D$wald05, diff_lfc1 = D$wald_lfc1))
}
write.csv(tab, file.path(BASE, "results/tables/BMATID_deseq2_axis_attribution.csv"), row.names = FALSE)

w(""); w("=================================================================="); w("  汇总表"); w("==================================================================")
w(sprintf("%-22s %8s %12s %14s %12s %10s", "方案", "removed", "体内 LRT", "体内并集|LFC|>1", "体外|LFC|>1", "体内保留"))
w(paste(rep("-", 84), collapse = ""))
b <- tab$prim_lfc1[1]
for (i in seq_len(nrow(tab))) {
  w(sprintf("%-22s %8d %12s %14d %12d %9.1f%%", tab$scheme[i], tab$removed[i],
            ifelse(is.na(tab$prim_LRT05[i]), "-", as.character(tab$prim_LRT05[i])),
            tab$prim_lfc1[i], tab$diff_lfc1[i], 100 * tab$prim_lfc1[i] / b))
}
w("")
w(sprintf("  体内基准（全基因，并集&|LFC|>1）= %d ; 体外基准 = %d  → 体外保留 %.1f%%",
          tab$prim_lfc1[1], tab$diff_lfc1[1], 100 * tab$diff_lfc1[1] / tab$prim_lfc1[1]))
close(con); cat("done\n")
