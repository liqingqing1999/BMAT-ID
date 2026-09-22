## 体外残余基因：预处理 + LFC 收缩 + 离群样本稳健性
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

out <- file.path(BASE, "logs", "deseq2_robust_out.txt")
con <- file(out, "w"); w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }
suppressMessages(library(DESeq2))

dd <- read.delim(gzfile(file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")),
                 header = TRUE, row.names = 1, check.names = FALSE)
rownames(dd) <- sub("\\..*$", "", rownames(dd)); dd <- dd[!duplicated(rownames(dd)), ]
cn <- colnames(dd)
donor <- sub("^H([0-9]+)_.*$", "\\1", cn)
location <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", cn)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", cn)
di <- which(state == "diff")

TARGET <- c(HOXB3 = "ENSG00000120093", HOXB4 = "ENSG00000182742", HOXB7 = "ENSG00000260027",
            ITGA8 = "ENSG00000077943", GPAT3 = "ENSG00000138678", CDKN2C = "ENSG00000123080",
            AHNAK2 = "ENSG00000185567", TLL1 = "ENSG00000038295")

run <- function(idx, tag) {
  cnt <- round(as.matrix(dd[, idx]))
  cnt <- cnt[rowSums(cnt) >= 10, , drop = FALSE]
  m <- data.frame(donor = factor(donor[idx]),
                  location = factor(location[idx], levels = c("epi", "meta", "scat")))
  dds <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  dds <- DESeq(dds, test = "LRT", reduced = ~ donor, quiet = TRUE)
  r <- results(dds)
  w(""); w("================================================================"); w("  ", tag)
  w("================================================================")
  w(sprintf("  LRT: tested=%d  FDR<.05=%d  FDR<.10=%d  raw p<.05=%d (%.1f%%)",
            sum(!is.na(r$padj)), sum(!is.na(r$padj) & r$padj < 0.05),
            sum(!is.na(r$padj) & r$padj < 0.10),
            sum(!is.na(r$pvalue) & r$pvalue < 0.05),
            100 * sum(!is.na(r$pvalue) & r$pvalue < 0.05) / max(1, sum(!is.na(r$padj)))))
  q <- quantile(r$pvalue, c(0, .25, .5, .75, 1), na.rm = TRUE)
  w(sprintf("  p 分位: %s", paste(sprintf("%.3f", q), collapse = " / ")))

  # Wald + 收缩（normal 型，不需要额外包）
  d2 <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  d2 <- DESeq(d2, quiet = TRUE)
  w("")
  w("  === 目标基因：Wald + LFC 收缩（type=normal） ===")
  w(sprintf("  %-8s %-19s %13s %13s %13s", "symbol", "ensembl", "meta_vs_epi", "scat_vs_epi", "scat_vs_meta"))
  for (nm in names(TARGET)) {
    e <- TARGET[[nm]]; cells <- c()
    for (cf in c("location_meta_vs_epi", "location_scat_vs_epi", "location_scat_vs_meta")) {
      if (!cf %in% resultsNames(d2)) { cells <- c(cells, "na"); next }
      rs <- lfcShrink(d2, coef = cf, type = "normal", quiet = TRUE)
      if (e %in% rownames(rs)) {
        cells <- c(cells, sprintf("%+.2f/%.2g", rs[e, "log2FoldChange"], rs[e, "padj"]))
      } else cells <- c(cells, "na")
    }
    w(sprintf("  %-8s %-19s %13s %13s %13s", nm, e, cells[1], cells[2], cells[3]))
  }
  w("  （格式：收缩后 log2FC / FDR）")

  # 收缩后按 FDR 排序的残余基因（体内最想看的）
  rs <- lfcShrink(d2, coef = "location_scat_vs_epi", type = "normal", quiet = TRUE)
  rs2 <- lfcShrink(d2, coef = "location_meta_vs_epi", type = "normal", quiet = TRUE)
  minsig <- pmin(rs$padj, rs2$padj, na.rm = TRUE)
  top <- order(minsig)[1:15]
  w("")
  w("  === 收缩后最强部位基因 top15（取两对比 FDR 较小者） ===")
  for (i in top) {
    if (is.na(minsig[i])) next
    w(sprintf("     %-19s baseMean=%9.1f  LFC(scat/epi)=%+6.2f  LFC(meta/epi)=%+6.2f  minFDR=%.3g",
              rownames(rs)[i], rs$baseMean[i], rs$log2FoldChange[i], rs2$log2FoldChange[i], minsig[i]))
  }
  invisible(list(r = r, n05 = sum(!is.na(r$padj) & r$padj < 0.05)))
}

w("==================================================================")
w("  体外残余部位基因稳健性（DESeq2, prefilter rowSums>=10）  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("==================================================================")
w("体外 12 样本：", paste(cn[di], collapse = ", "))
S1 <- run(di, "S1 全部 12 个体外样本")
S2 <- run(di[!(cn[di] %in% "H410_epi_diff_adip")], "S2 剔除离群样本 H410_epi_diff_adip（11 个）")
w("")
w("==================================================================")
w("  小结")
w("==================================================================")
w(sprintf("  S1 (n=12)  LRT FDR<0.05 = %d", S1$n05))
w(sprintf("  S2 (n=11)  LRT FDR<0.05 = %d", S2$n05))
close(con); cat("done\n")
