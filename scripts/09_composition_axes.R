
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
logf <- file.path(BASE, "logs/axis_out.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }

p <- file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
cpm0 <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep <- rowSums(cpm0 >= 1) >= 4
dd <- d0[keep, ]
Lc <- log2(t(t(as.matrix(dd))/colSums(as.matrix(dd))*1e6) + 1)
pi <- which(state == "primary")

# 三轴标记基因：symbol -> Ensembl，从 BMATID_marker_ensg.tsv 权威读入（勿手写）
mk <- read.delim(file.path(BASE, "results/tables/BMATID_marker_ensg.tsv"), stringsAsFactors=FALSE)
mk <- mk[mk$ensg != "NA", ]
AXES <- split(mk$ensg, mk$axis)
names(AXES) <- c("Hemato", "Bone", "Adipo")[match(names(AXES), c("Hemato","Bone","Adipo"))]
AXES <- AXES[c("Hemato", "Bone", "Adipo")]
SYM  <- split(mk$symbol, mk$axis)

ens <- sub("\\..*$", "", rownames(Lc))
idxmap <- tapply(seq_along(ens), ens, function(x) x)
score <- function(ids) {
  rows <- unique(unlist(lapply(ids, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
  rows <- rows[!is.na(rows)]
  if (!length(rows)) return(rep(NA, length(pi)))
  Z <- t(scale(t(Lc[rows, pi, drop=FALSE])))
  colMeans(Z, na.rm=TRUE)
}
w("=== 三轴标记基因在体内三部位的表达（log2CPM 均值）===")
for (ax in names(AXES)) {
  w(sprintf("--- %s axis (%s)", ax, paste(SYM[[ax]], collapse=", ")))
  for (nm in SYM[[ax]]) {
    id <- mk$ensg[mk$symbol == nm & mk$axis == ax]
    if (!length(id) || !(id[1] %in% names(idxmap))) { w(sprintf("   %-8s (not in matrix)", nm)); next }
    i <- idxmap[[id[1]]][1]
    v <- sapply(c("epi","meta","scat"), function(s) mean(Lc[i, pi[loc[pi]==s]]))
    rng <- max(v) - min(v)
    flag <- if (rng > 1.0) "  <<<" else ""
    w(sprintf("   %-8s EPI=%5.2f META=%5.2f SCAT=%5.2f  range=%.2f%s", nm, v[1], v[2], v[3], rng, flag))
  }
}

S <- sapply(names(AXES), function(a) score(AXES[[a]]))
rownames(S) <- colnames(Lc)[pi]
w("")
w("=== 轴分数之间的相关（体内 12 样本）===")
writeLines(paste(round(cor(S), 3)[1, ], collapse="  "), con)
print(round(cor(S), 3))

# 每个基因与三轴的相关
rc <- sapply(colnames(S), function(a) suppressWarnings(as.numeric(cor(t(Lc[, pi]), S[, a]))))
colnames(rc) <- colnames(S); rownames(rc) <- rownames(Lc)
w("")
w("=== 与三轴的 |r| 分布 ===")
for (a in colnames(rc)) {
  v <- abs(rc[, a]); v <- v[!is.na(v)]
  w(sprintf("  %-7s  |r|>0.5: %5d  >0.6: %5d  >0.7: %5d  >0.8: %4d   (%d genes)",
            a, sum(v>0.5), sum(v>0.6), sum(v>0.7), sum(v>0.8), length(v)))
}

fitlimma <- function(gsel) {
  dd2 <- dd[gsel, ]
  idx <- pi
  m <- as.matrix(dd2[, idx]); m[is.na(m)] <- 0
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  ce <- grep("^lc", colnames(des), value=TRUE)
  sum(topTable(fit, coef=ce, number=Inf, sort.by="none")$adj.P.Val < 0.05)
}

w("")
w("=== 逐轴剔除后，部位效应 FDR<0.05 还剩多少（先过滤后拟合口径）===")
w(sprintf("  %-44s genes=%5d  siteDEG=%4d  (%.0f%% of 2712)", "基准：全部基因", nrow(dd), fitlimma(rep(TRUE, nrow(dd))), 100))
sel <- function(...) {
  bad <- rep(FALSE, nrow(dd))
  for (v in list(...)) bad <- bad | (!is.na(v) & abs(v) > 0.7)
  !bad
}
combos <- list(
  "剔除 Hemato 轴 (|r|>0.7)" = sel(rc[,"Hemato"]),
  "剔除 Hemato + Bone 轴"    = sel(rc[,"Hemato"], rc[,"Bone"]),
  "剔除 Bone 轴 单轴"        = sel(rc[,"Bone"]),
  "剔除 Hemato+Bone+Adipo"   = sel(rc[,"Hemato"], rc[,"Bone"], rc[,"Adipo"]))
for (tag in names(combos)) {
  g <- combos[[tag]]
  n <- fitlimma(g)
  w(sprintf("  %-44s genes=%5d  siteDEG=%4d  (%.0f%%)", tag, sum(g), n, 100*n/2712))
}

# 残差 494 与三轴的关联：这些基因本身像哪一种细胞
rec <- read.delim(file.path(BASE, "results/tables/BMATID_residual494_schemeA.tsv"), check.names=FALSE)
rec <- rec[rec$ensg %in% rownames(rc), ]
w("")
w("=== 残差 494 基因与三轴的相关（中位数 / |r|>0.5 占比）===")
for (a in colnames(rc)) {
  v <- rc[rec$ensg, a]
  w(sprintf("  %-7s median r = %+.3f   |r|>0.5: %d/%d (%.0f%%)", a,
            median(v, na.rm=TRUE), sum(abs(v)>0.5, na.rm=TRUE), length(v),
            100*mean(abs(v)>0.5, na.rm=TRUE)))
}

# ---- 随机轴对照：剔除同样数量、但与任何细胞组成无关的基因 ----
w("")
w("=== 直接归因：2712 个部位 DEG 里有多少与组成轴共变 ===")
sig <- read.csv(file.path(BASE, "results/tables/BMATID_primary_cd45_annotated.csv"), check.names=FALSE)
sig <- sig[sig$FDR_prim < 0.05, ]
gg <- intersect(sig$ensg, rownames(rc))
h <- abs(rc[gg, "Hemato"]) > 0.7; b <- abs(rc[gg, "Bone"]) > 0.7
w(sprintf("  2712 部位 DEG 中：|r|>0.7 vs Hemato = %d (%.0f%%)", sum(h, na.rm=TRUE), 100*mean(h, na.rm=TRUE)))
w(sprintf("                      |r|>0.7 vs Bone   = %d (%.0f%%)", sum(b, na.rm=TRUE), 100*mean(b, na.rm=TRUE)))
w(sprintf("                      任一轴             = %d (%.0f%%)", sum(h|b, na.rm=TRUE), 100*mean(h|b, na.rm=TRUE)))
w(sprintf("                      与三轴都不相关(|r|max<=0.5) = %d (%.0f%%)",
          sum(apply(abs(rc[gg, ]), 1, max, na.rm=TRUE) <= 0.5),
          100*mean(apply(abs(rc[gg, ]), 1, max, na.rm=TRUE) <= 0.5)))

w("")
w("=== 等数量随机剔除对照（与真实剔除同数量，各 8 次）===")
set.seed(2026)
run_nmatch <- function(n_drop, tag) {
  v <- c()
  for (k in 1:8) {
    g <- rep(TRUE, nrow(dd)); g[sample(nrow(dd), min(n_drop, nrow(dd) - 100))] <- FALSE
    v <- c(v, fitlimma(g))
  }
  w(sprintf("  %-40s siteDEG = %s  (median %.0f)", tag, paste(v, collapse=","), median(v)))
  invisible(v)
}
run_nmatch(sum(abs(rc[, "Hemato"]) > 0.7, na.rm = TRUE), "随机剔 2459（= 造血轴数量）")
run_nmatch(sum(apply(abs(rc[, c("Hemato","Bone")]), 1, max, na.rm = TRUE) > 0.7), "随机剔 2995（= 造血+骨轴数量）")
run_nmatch(0.7 * nrow(dd), "随机剔 70% 全部基因（极端对照）")

w("")
w("=== 剔除 Hemato+Bone 后唯一幸存的部位基因 ===")
bad <- rep(FALSE, nrow(dd))
for (a in c("Hemato", "Bone")) bad <- bad | (!is.na(rc[, a]) & abs(rc[, a]) > 0.7)
dd2 <- dd[!bad, ]
idx <- pi
m <- as.matrix(dd2[, idx]); m[is.na(m)] <- 0
dn <- factor(donor[idx]); lc <- factor(loc[idx], levels = c("meta","epi","scat"))
des <- model.matrix(~ dn + lc)
v <- voom(m, des, lib.size = colSums(m), normalize.method = "quantile")
fit <- eBayes(lmFit(v, des))
ce <- grep("^lc", colnames(des), value = TRUE)
tt <- topTable(fit, coef = ce, number = Inf, sort.by = "none")
sv <- tt[tt$adj.P.Val < 0.05, , drop = FALSE]
w(sprintf("  幸存基因数 = %d", nrow(sv)))
if (nrow(sv)) writeLines(paste(rownames(sv), sprintf("F=%.2f", sv$F), sprintf("FDR=%.3f", sv$adj.P.Val)), con)


saveRDS(rc, file.path(BASE, "results/tables/BMATID_axis_correlations.rds"))
saveRDS(S, file.path(BASE, "results/tables/BMATID_axis_scores.rds"))
close(con)
cat("ok\n")
