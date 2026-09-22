## DESeq2 复核（精简版，后台跑）：体内 vs 体外部位效应 + 组成轴归因
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

out <- file.path(BASE, "logs", "deseq2_out.txt")
con <- file(out, "w"); w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }
suppressMessages({ library(DESeq2) })

w("==================================================================")
w("  DESeq2 复核（精简单轮）  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
w("  DESeq2 ", as.character(packageVersion("DESeq2")), "   model: full ~ donor+location / reduced ~ donor (LRT)")
w("==================================================================")

f <- file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")
dd <- read.delim(gzfile(f), header = TRUE, row.names = 1, check.names = FALSE)
rownames(dd) <- sub("\\..*$", "", rownames(dd))
dd <- dd[!duplicated(rownames(dd)), ]
cn <- colnames(dd)
donor <- sub("^H([0-9]+)_.*$", "\\1", cn)
location <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", cn)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", cn)
prim_i <- which(state == "primary"); diff_i <- which(state == "diff")
w(sprintf("counts: %d genes x %d samples | primary %d, diff %d",
          nrow(dd), ncol(dd), length(prim_i), length(diff_i)))

fit_one <- function(idx, tag, gene_subset = NULL, do_pair = FALSE) {
  cnt <- round(as.matrix(dd[, idx])); m <- data.frame(
    donor = factor(donor[idx]),
    location = factor(location[idx], levels = c("epi", "meta", "scat")))
  if (!is.null(gene_subset)) cnt <- cnt[intersect(rownames(cnt), gene_subset), , drop = FALSE]
  cnt <- cnt[rowSums(cnt) >= 10, , drop = FALSE]
  dds <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  dds <- DESeq(dds, test = "LRT", reduced = ~ donor, quiet = TRUE)
  r <- results(dds); padj <- r$padj
  n05 <- sum(!is.na(padj) & padj < 0.05); n01 <- sum(!is.na(padj) & padj < 0.01)
  nraw <- sum(!is.na(r$pvalue) & r$pvalue < 0.05); nt <- sum(!is.na(padj))
  w("")
  w(sprintf("--- %-40s tested=%6d  FDR<.05=%5d  FDR<.01=%5d  raw p<.05=%5d (%.1f%%)",
            tag, nt, n05, n01, nraw, 100 * nraw / max(1, nt)))
  q <- quantile(r$pvalue, c(0, .25, .5, .75, 1), na.rm = TRUE)
  w(sprintf("    p 分位(0/.25/.5/.75/1): %s", paste(sprintf("%.3f", q), collapse = " / ")))
  sig <- r[!is.na(r$padj) & r$padj < 0.05, , drop = FALSE]
  if (nrow(sig) > 0 && nrow(sig) <= 40) {
    w("    全部显著基因："); o <- order(sig$padj)
    for (i in o) w(sprintf("      %-18s baseMean=%9.1f  log2FC=%+7.2f  FDR=%.3g",
                           rownames(sig)[i], sig$baseMean[i], sig$log2FoldChange[i], sig$padj[i]))
  }
  if (do_pair) {
    for (cf in c("location_meta_vs_epi", "location_scat_vs_epi", "location_scat_vs_meta")) {
      if (!cf %in% resultsNames(dds)) next
      rp <- results(dds, name = cf)
      w(sprintf("    %-26s FDR<.05=%5d  FDR<.05 & |FC|>2=%5d  (NA=%d)",
                cf, sum(rp$padj < 0.05, na.rm = TRUE),
                sum(rp$padj < 0.05 & abs(rp$log2FoldChange) > 2, na.rm = TRUE), sum(is.na(rp$padj))))
    }
  }
  invisible(list(dds = dds, res = r, n05 = n05, ngene = nt))
}

w(""); w("############ 1. 体内 primary（全基因）")
P <- fit_one(prim_i, "1 体内 primary  LRT", do_pair = TRUE)
w(""); w("############ 2. 体外 diff（全基因）")
D <- fit_one(diff_i, "2 体外 diff     LRT", do_pair = TRUE)

w(""); w("############ 3. 组成轴归因（剔除与组成轴共变的基因后重算）")
## The correlation matrix is the one written by 09_composition_axes.R, whose axes
## are Hemato / Bone / Adipo. Version suffixes are stripped from its rownames
## before gene sets are compared: the counts matrix above is read version-stripped
## (sub("\\..*$", "", ...)), so a setdiff() on the raw rownames would match nothing.
rcf <- file.path(BASE, "results/tables/BMATID_axis_correlations.rds")
res_summary <- data.frame()
no_version <- function(x) sub("\\..*$", "", x)
if (file.exists(rcf)) {
  rc <- readRDS(rcf)
  w("   rc 维度: ", nrow(rc), " x ", ncol(rc), " | 列: ", paste(colnames(rc), collapse = ", "))
  axes <- intersect(c("Hemato", "Bone", "Adipo"), colnames(rc))
  if (!length(axes)) {
    w("   [warn] 未识别到组成轴列，本节跳过")
  } else {
    removed_for <- function(cols) {
      g <- unique(unlist(lapply(cols, function(a) {
        # NB: rc is a matrix, so a column must be taken with rc[, a] -- double-bracket
        # indexing would be linear and raise "subscript out of bounds".
        v <- rc[, a]; names(v) <- no_version(rownames(rc))
        names(v)[!is.na(v) & abs(v) > 0.7]
      })))
      g[!is.na(g)]
    }
    combos <- lapply(axes, function(a) a); names(combos) <- axes
    if (all(c("Hemato", "Bone") %in% axes)) combos[["Hemato+Bone"]] <- c("Hemato", "Bone")
    for (nm in names(combos)) {
      rm_genes <- removed_for(combos[[nm]])
      keep <- setdiff(rownames(dd), rm_genes)
      resP <- fit_one(prim_i, sprintf("3 体内 primary  剔 %s (n=%d)", nm, length(rm_genes)), gene_subset = keep)
      resD <- fit_one(diff_i, sprintf("3 体外 diff     剔 %s (n=%d)", nm, length(rm_genes)), gene_subset = keep)
      res_summary <- rbind(res_summary, data.frame(axis = nm, removed = length(rm_genes),
                                                   n_prim = resP$n05, n_prim_test = resP$ngene,
                                                   n_diff = resD$n05, n_diff_test = resD$ngene))
    }
  }
} else w("   [warn] 未找到 ", rcf)

w(""); w("==================================================================")
w("  汇总：DESeq2 口径下的部位效应")
w("==================================================================")
w(sprintf("%-14s %8s %12s %14s %12s %14s", "axis", "removed", "体内 FDR<.05", "/tested", "体外 FDR<.05", "/tested"))
w(paste(rep("-", 78), collapse = ""))
w(sprintf("%-14s %8s %12d %14s %12d %14s", "（全基因）", 0, P$n05, sprintf("(%d)", P$ngene), D$n05, sprintf("(%d)", D$ngene)))
if (nrow(res_summary) > 0)
  for (i in seq_len(nrow(res_summary))) {
    s <- res_summary[i, ]
    w(sprintf("%-14s %8d %12d %14s %12d %14s", s$axis, s$removed, s$n_prim,
              sprintf("(%d)", s$n_prim_test), s$n_diff, sprintf("(%d)", s$n_diff_test)))
  }

## 导出体外那 14 个残余基因的完整结果（含 pairwise）
w(""); w("############ 4. 体外残余显著基因的逐对比 log2FC")
sigD <- D$res[!is.na(D$res$padj) & D$res$padj < 0.05, , drop = FALSE]
if (nrow(sigD) > 0) {
  cnt <- round(as.matrix(dd[, diff_i])); cnt <- cnt[rownames(sigD), , drop = FALSE]
  m <- data.frame(donor = factor(donor[diff_i]),
                  location = factor(location[diff_i], levels = c("epi", "meta", "scat")))
  d2 <- DESeqDataSetFromMatrix(cnt, m, design = ~ donor + location)
  d2 <- DESeq(d2, quiet = TRUE)
  for (cf in c("location_meta_vs_epi", "location_scat_vs_epi", "location_scat_vs_meta")) {
    if (!cf %in% resultsNames(d2)) next
    rr <- results(d2, name = cf)
    w("  -- ", cf)
    for (i in seq_len(nrow(rr)))
      w(sprintf("      %-18s log2FC=%+7.2f  FDR=%.3g", rownames(rr)[i], rr$log2FoldChange[i], rr$padj[i]))
  }
  write.csv(as.data.frame(sigD), file.path(BASE, "results/tables/BMATID_deseq2_diff_residual_genes.csv"))
  w("  -> 已导出 results/tables/BMATID_deseq2_diff_residual_genes.csv")
}
close(con); cat("done\n")
