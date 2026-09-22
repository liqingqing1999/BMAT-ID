
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
outdir <- file.path(BASE, "results", "tables")
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)

p <- file.path(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)

cpm0 <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep <- rowSums(cpm0 >= 1) >= 4
dd <- d0[keep, ]
Lc <- log2(t(t(as.matrix(dd))/colSums(as.matrix(dd))*1e6) + 1)

res <- list()
for (st in c("primary","diff")) {
  idx <- which(state == st)
  m <- as.matrix(dd[, idx]); m[is.na(m)] <- 0
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  ce <- grep("^lc", colnames(des), value=TRUE)
  tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
  cm <- makeContrasts(EM = lcepi, ES = lcepi - lcscat, MS = lcscat, levels=des)
  f2 <- eBayes(contrasts.fit(lmFit(v, des), cm))
  lf <- sapply(colnames(cm), function(cc) topTable(f2, coef=cc, number=Inf, sort.by="none")$logFC)
  res[[st]] <- data.frame(
    ensg = rownames(tt),
    F = tt$F, FDR = tt$adj.P.Val, P = tt$P.Value,
    lfc_EM = lf[,"EM"], lfc_ES = lf[,"ES"], lfc_MS = lf[,"MS"],
    maxlfc = apply(abs(lf), 1, max),
    stringsAsFactors = FALSE)
  # site means (log2CPM)
  for (g in c("epi","meta","scat")) {
    gi <- idx[loc[idx] == g]
    res[[st]][[paste0("mean_", toupper(g))]] <- rowMeans(Lc[, gi, drop=FALSE])
  }
}

out <- data.frame(
  ensg = res$primary$ensg,
  F_prim = res$primary$F, FDR_prim = res$primary$FDR,
  F_diff = res$diff$F,    FDR_diff = res$diff$FDR,
  maxlfc_prim = res$primary$maxlfc, maxlfc_diff = res$diff$maxlfc,
  EPI_p = res$primary$mean_EPI, META_p = res$primary$mean_META, SCAT_p = res$primary$mean_SCAT,
  EPI_d = res$diff$mean_EPI,    META_d = res$diff$mean_META,    SCAT_d = res$diff$mean_SCAT,
  stringsAsFactors = FALSE)

f <- file.path(outdir, "BMATID_limma_allgenes.csv")
write.csv(out, f, row.names=FALSE)
cat("saved:", f, " rows:", nrow(out), "\n")
cat("FDR_prim<0.05:", sum(out$FDR_prim<0.05), " FDR_diff<0.05:", sum(out$FDR_diff<0.05), "\n")
