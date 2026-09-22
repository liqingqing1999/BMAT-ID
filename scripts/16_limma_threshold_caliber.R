
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

logf <- file.path(BASE, "logs", "caliber_out.txt")
con <- file(logf, "w")
w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con) }

p <- file.path(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
d <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
w("counts dim: ", paste(dim(d), collapse=" x "))
samp <- colnames(d)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
w("states: ", paste(table(state), collapse=" / "))
w("locs  : ", paste(table(loc), collapse=" / "))
w("donors: ", paste(table(donor), collapse=" / "))
w("")

# ---- 原文口径的低表达过滤：每组(GSM 全 24)中 >=90% 样本 cpm >= 50 太严，先做温和过滤 ----
cpmM <- t(t(as.matrix(d)) / colSums(as.matrix(d)) * 1e6)
keep <- rowSums(cpmM >= 1) >= 4      # 至少在 4 个样本中 cpm>=1
w("genes total = ", nrow(d), " ; after cpm>=1 in >=4 samples: ", sum(keep))
d <- d[keep, ]

run_state <- function(st, tag) {
  idx <- which(state == st)
  m <- as.matrix(d[, idx, drop=FALSE]); m[is.na(m)] <- 0
  dn <- factor(donor[idx])
  lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  colnames(des) <- make.names(colnames(des))
  w(""); w("======== ", tag, "  n=", length(idx), " ========")
  w("design cols: ", paste(colnames(des), collapse=", "))
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  # 找 loc 系数名
  ce <- grep("^lc", colnames(des), value=TRUE)
  w("loc coefs: ", paste(ce, collapse=", "))

  # localization 总效应 F 检验
  f <- topTable(fit, coef=ce, number=Inf, sort.by="F")
  nF05 <- sum(f$adj.P.Val < 0.05); nF01 <- sum(f$adj.P.Val < 0.01)
  w("LOCALIZATION F-test : FDR<0.05 = ", nF05, " ; FDR<0.01 = ", nF01)
  if (nrow(f) > 0) w("   top1: ", rownames(f)[1], " F=", round(f$F[1], 1), " FDR=", signif(f$adj.P.Val[1], 3))

  # pairwise
  cm <- makeContrasts(EPIvMETA = lcepi,
                      EPIvSCAT = lcepi - lcscat,
                      METAvSCAT = lcscat,
                      levels = des)
  fit2 <- eBayes(contrasts.fit(lmFit(v, des), cm))
  allDEG <- list()
  for (cc in colnames(cm)) {
    tt <- topTable(fit2, coef=cc, number=Inf, sort.by="none")
    a05 <- sum(tt$adj.P.Val < 0.05)
    a05fc <- sum(tt$adj.P.Val < 0.05 & abs(tt$logFC) > 2)
    w(sprintf("  %-9s : FDR<0.05 = %4d ; FDR<0.05 & |log2FC|>2 = %4d", cc, a05, a05fc))
    allDEG[[cc]] <- rownames(tt)[tt$adj.P.Val < 0.05]
  }
  # 近似 upsert exclusive：只在一个 pairwise 中显著的基因
  u <- unique(unlist(allDEG))
  cnt <- table(unlist(allDEG))
  w("  union of pairwise DEG (FDR<0.05): ", length(u))
  w("  genes sig in exactly 1 contrast (approx 'exclusive'): ", sum(cnt == 1))
  w("     by which: ", paste(names(cnt)[cnt == 1], collapse=""), " -> ")
  for (cc in colnames(cm)) w("        exclusive-to-", cc, ": ", length(setdiff(allDEG[[cc]], unlist(allDEG[setdiff(colnames(cm), cc)]))))
  w("  shared by all 3 contrasts: ", sum(cnt == 3))
  invisible(list(f = f, deg = allDEG))
}

r_prim <- run_state("primary", "PRIMARY (in vivo)")
r_diff <- run_state("diff",    "DIFF (in vitro)")

w(""); w("======== 汇总对比 ========")
w("                      in vivo(primary)   in vitro(diff)")
w("localization F FDR<0.05 :  %4d            %4d", sum(r_prim$f$adj.P.Val<0.05), sum(r_diff$f$adj.P.Val<0.05))
w("")
w("NOTE: 原文口径 = DESeq2 ~donor+localization, FDR<0.05 且 |log2FC|>2")
w("NOTE: 原文报的是 upsert 'specifically expressed' : in vitro EPI=2 META=11 SCAT=29 ; in vivo EPI=234 META=918 SCAT=1290")
close(con)
cat("\n[log]", logf, "\n")
