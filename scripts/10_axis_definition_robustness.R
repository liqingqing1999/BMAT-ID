
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
logf <- file.path(BASE, "logs/axis_robust_out.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }

p <- file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$","\\1",samp); loc <- sub("^H[0-9]+_([a-z]+)_.*$","\\1",samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$","\\1",samp)
cpm0 <- t(t(as.matrix(d0))/colSums(as.matrix(d0))*1e6); keep <- rowSums(cpm0>=1)>=4
dd <- d0[keep,]; Lc <- log2(t(t(as.matrix(dd))/colSums(as.matrix(dd))*1e6)+1); pi <- which(state=="primary")
ens <- sub("\\..*$","",rownames(Lc)); idxmap <- tapply(seq_along(ens), ens, function(x) x)

mk <- read.delim(file.path(BASE,"results/tables/BMATID_marker_ensg.tsv"), stringsAsFactors=FALSE)
SYMENS <- setNames(mk$ensg, mk$symbol)

score <- function(syms) {
  rows <- unique(unlist(lapply(syms, function(s) { e <- SYMENS[[s]]; if (!is.null(e) && e %in% names(idxmap)) idxmap[[e]] else NULL })))
  Z <- t(scale(t(Lc[rows, pi, drop=FALSE]))); colMeans(Z, na.rm=TRUE)
}
fitlimma <- function(gsel) {
  m <- as.matrix(dd[gsel, pi]); m[is.na(m)] <- 0
  dn <- factor(donor[pi]); lc <- factor(loc[pi], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des)); ce <- grep("^lc", colnames(des), value=TRUE)
  sum(topTable(fit, coef=ce, number=Inf, sort.by="none")$adj.P.Val < 0.05)
}
AX <- list(
  H_all  = c("PTPRC","CD14","LYZ","CD3E","MS4A1","MZB1","JCHAIN","CSF1R","CD19","CD79A"),
  H_noIg = c("PTPRC","CD14","LYZ","CD3E","MS4A1","CD19","CD79A","CSF1R"),        # 剔免疫球蛋白/浆细胞
  B_all  = c("BGLAP","ALPL","IBSP","SPP1","MEPE","DMP1","RUNX2","PHEX","COL1A1","TMEM119","ENPP1"),
  B_noCol= c("BGLAP","ALPL","IBSP","SPP1","MEPE","DMP1","RUNX2","PHEX","TMEM119","ENPP1"),  # 剔胶原
  B_small= c("MEPE","DMP1","PHEX"))                                              # 仅骨细胞三基因
S <- sapply(AX, score); rownames(S) <- colnames(Lc)[pi]
rc <- sapply(colnames(S), function(a) suppressWarnings(as.numeric(cor(t(Lc[,pi]), S[,a]))))
colnames(rc) <- colnames(S)

w("=== 各轴分数之间的相关 ===")
print(round(cor(S), 3))
w(sprintf("基准（全基因）: siteDEG = %d", fitlimma(rep(TRUE, nrow(dd)))))
w("")
w("=== 不同轴定义下：剔除 |r|>0.7 后的部位效应 ===")
for (a in c("H_all","H_noIg","B_all","B_noCol","B_small")) {
  sel <- !(is.na(rc[,a]) | abs(rc[,a]) > 0.7)
  w(sprintf("  %-28s removed=%5d  siteDEG=%4d", a, sum(!sel), fitlimma(sel)))
}

w("")
w("=== 组合：只用非 Ig 的造血轴 + 剔胶原的骨轴 ===")
bad <- rep(FALSE, nrow(dd))
for (a in c("H_noIg","B_noCol")) bad <- bad | (!is.na(rc[,a]) & abs(rc[,a])>0.7)
w(sprintf("  removed=%d  siteDEG=%d", sum(bad), fitlimma(!bad)))
bad2 <- bad | (!is.na(rc[,"B_small"]) & abs(rc[,"B_small"])>0.7)
w(sprintf("  再并 B_small: removed=%d  siteDEG=%d", sum(bad2), fitlimma(!bad2)))
close(con); cat("ok\n")
