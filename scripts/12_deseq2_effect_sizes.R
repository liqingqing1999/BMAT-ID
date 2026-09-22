## 最终定量：DESeq2 + Wald + LFC 收缩，加效应量阈值，体内 vs 体外
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

out <- file.path(BASE, "logs", "deseq2_final_out.txt")
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

analyze <- function(idx, tag, save_prefix = NULL) {
  cnt <- round(as.matrix(dd[, idx])); cnt <- cnt[rowSums(cnt) >= 10, , drop = FALSE]
  m <- data.frame(donor = factor(donor[idx]),
                  location = factor(location[idx], levels = c("epi", "meta", "scat")))
  dds <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  dds <- DESeq(dds, quiet = TRUE)
  res <- list()
  for (ct in CONTRASTS) {
    nm <- sprintf("%s_vs_%s", ct[1], ct[2])
    r <- lfcShrink(dds, contrast = c("location", ct[1], ct[2]), type = "normal", quiet = TRUE)
    r <- as.data.frame(r)
    r$gene <- rownames(r); r$contrast <- nm
    res[[nm]] <- r
  }
  # 并集显著（任一对比 FDR<0.05）
  allg <- unique(unlist(lapply(res, function(x) x$gene)))
  best <- data.frame(gene = allg, baseMean = NA_real_, maxabsLFC = 0,
                     minFDR = NA_real_, nB = 0, stringsAsFactors = FALSE)
  for (i in seq_along(allg)) {
    g <- allg[i]; mf <- NA_real_; ml <- 0; nb <- 0; bm <- NA_real_
    for (nm in names(res)) {
      x <- res[[nm]]; j <- which(x$gene == g)
      if (!length(j)) next
      if (is.na(bm)) bm <- x$baseMean[j]
      if (!is.na(x$padj[j])) {
        mf <- if (is.na(mf)) x$padj[j] else min(mf, x$padj[j])
        if (x$padj[j] < 0.05) { nb <- nb + 1; ml <- max(ml, abs(x$log2FoldChange[j])) }
      }
    }
    best$baseMean[i] <- bm; best$minFDR[i] <- mf; best$maxabsLFC[i] <- ml; best$nB[i] <- nb
  }
  w(""); w("================================================================"); w("  ", tag)
  w("================================================================")
  w(sprintf("  基因数（prefilter 后）= %d", nrow(cnt)))
  w(sprintf("  %-46s %8s", "判据", "基因数"))
  w(paste(rep("-", 58), collapse = ""))
  for (th in c(0, 0.25, 0.5, 1, 1.5, 2)) {
    n <- sum(best$nB > 0 & best$maxabsLFC > th, na.rm = TRUE)
    lab <- if (th == 0) "FDR<0.05（任一对比，不看效应量）" else sprintf("FDR<0.05 且 |log2FC| > %.2f", th)
    w(sprintf("  %-46s %8d", lab, n))
  }
  n_fdr10 <- sum(best$minFDR < 0.10 & best$maxabsLFC > 1, na.rm = TRUE)
  w(sprintf("  %-46s %8d", "FDR<0.10 且 |log2FC| > 1", n_fdr10))
  w("")
  sel <- best[!is.na(best$minFDR) & best$minFDR < 0.05 & best$maxabsLFC > 1, ]
  sel <- sel[order(sel$minFDR), ]
  w(sprintf("  === FDR<0.05 且 |log2FC|>1 的基因（n=%d）===", nrow(sel)))
  show <- head(sel, 30)
  for (i in seq_len(nrow(show)))
    w(sprintf("     %-19s baseMean=%9.1f  max|LFC|=%5.2f  minFDR=%.3g  nContrast=%d",
              show$gene[i], show$baseMean[i], show$maxabsLFC[i], show$minFDR[i], show$nB[i]))
  if (!is.null(save_prefix)) {
    write.csv(best, file.path(BASE, sprintf("results/tables/%s.csv", save_prefix)), row.names = FALSE)
    for (nm in names(res))
      write.csv(res[[nm]], file.path(BASE, sprintf("results/tables/%s_%s.csv", save_prefix, nm)), row.names = FALSE)
    w("  -> 已导出 results/tables/", save_prefix, "*.csv")
  }
  best
}

w("==================================================================")
w("  DESeq2 + LFC 收缩 最终定量  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("  model ~donor+location; Wald + lfcShrink(type=normal); prefilter rowSums>=10")
w("==================================================================")
P <- analyze(which(state == "primary"), "体内 primary (n=12)", "BMATID_deseq2_primary_shrunk")
D <- analyze(which(state == "diff"),    "体外 diff    (n=12)", "BMATID_deseq2_diff_shrunk")

w(""); w("=================================================================="); w("  直接对照"); w("==================================================================")
getn <- function(b, th) sum(b$nB > 0 & b$maxabsLFC > th, na.rm = TRUE)
w(sprintf("%-34s %10s %10s %10s", "判据", "体内", "体外", "体外/体内"))
for (th in c(0, 0.5, 1, 1.5, 2)) {
  lab <- if (th == 0) "FDR<0.05（不看效应量）" else sprintf("FDR<0.05 & |log2FC|>%.2f", th)
  a <- getn(P, th); b <- getn(D, th)
  w(sprintf("%-34s %10d %10d %9.1f%%", lab, a, b, 100 * b / max(1, a)))
}
close(con); cat("done\n")
