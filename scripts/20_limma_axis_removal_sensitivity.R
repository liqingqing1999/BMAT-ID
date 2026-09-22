
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

suppressMessages(library(limma))
logf <- file.path(BASE, "logs", "sens_out.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con) }

p <- file.path(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)

cpm0 <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep <- rowSums(cpm0 >= 1) >= 4
dd <- d0[keep, ]
w("genes kept: ", nrow(dd))

w(""); w("===== library size 明细（升序）=====")
ls <- colSums(d0)
o <- order(ls)
for (i in o) w(sprintf("  %-26s %-8s %12s", samp[i], state[i], format(ls[i], big.mark=",")))

locF <- function(use, tag) {
  idx <- which(use)
  m <- as.matrix(dd[, idx]); m[is.na(m)] <- 0
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  ce <- grep("^lc", colnames(des), value=TRUE)
  tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
  # pairwise 最大 |log2FC|（EPI-META / EPI-SCAT / META-SCAT）
  cm <- makeContrasts(EM = lcepi, ES = lcepi - lcscat, MS = lcscat, levels=des)
  f2 <- eBayes(contrasts.fit(lmFit(v, des), cm))
  lfc <- sapply(colnames(cm), function(cc) abs(topTable(f2, coef=cc, number=Inf, sort.by="none")$logFC))
  maxlfc <- apply(lfc, 1, max)
  res <- list(n05 = sum(tt$adj.P.Val < 0.05),
              n05fc2 = sum(tt$adj.P.Val < 0.05 & maxlfc > 2),
              medF = median(tt$F), q90F = quantile(tt$F, 0.9), maxF = max(tt$F),
              nP05 = sum(tt$P.Value < 0.05), n = length(idx))
  w(sprintf("  %-40s n=%2d | FDR<.05=%5d | FDR<.05&|lfc|>2=%5d | p<.05=%5d | F: med=%.2f q90=%.2f max=%.1f",
            tag, res$n, res$n05, res$n05fc2, res$nP05, res$medF, res$q90F, res$maxF))
  invisible(res)
}

w(""); w("===== 敏感性分析 =====")
w("  [S1 基准：全部 4 供体]")
locF(state == "primary", "S1 primary (all 4 donors)")
locF(state == "diff",    "S1 diff    (all 4 donors)")

# 找 diff 最差样本
bad <- which(state == "diff")[which.min(ls[state == "diff"])]
w("")
w(sprintf("  [S2 剔除 diff 最低深度样本 %s (%s reads)]", samp[bad], format(ls[bad], big.mark=",")))
locF(state == "primary", "S2 primary (unchanged)")
locF(state == "diff" & seq_along(state) != bad, "S2 diff    (drop 1 low-depth)")

# 剔除该样本的整个供体（保持配对平衡）
bd <- donor[bad]
w(""); w(sprintf("  [S3 剔除供体 %s（两组各剩 3 供体 9 样本，功效对等）]", bd))
locF(state == "primary" & donor != bd, "S3 primary (3 donors)")
locF(state == "diff"    & donor != bd, "S3 diff    (3 donors)")

# 严格阈值：只要 |log2FC|>1 也更松？做 1 和 2 两档
w(""); w("===== 部位效应量压缩比（用于量化'抹平'）=====")
r1 <- locF(state == "primary", "primary")
r2 <- locF(state == "diff",    "diff")
w(sprintf("  F 中位数  : %.3f -> %.3f  (%.1f%%)", r1$medF, r2$medF, 100*r2$medF/r1$medF))
w(sprintf("  F 90 分位 : %.3f -> %.3f  (%.1f%%)", r1$q90F, r2$q90F, 100*r2$q90F/r1$q90F))
w(sprintf("  F 最大值  : %.2f -> %.2f  (%.1f%%)", r1$maxF, r2$maxF, 100*r2$maxF/r1$maxF))

close(con); cat("\n[log]", logf, "\n")
