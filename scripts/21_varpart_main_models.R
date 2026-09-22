## Variance decomposition (variancePartition) of the anatomical-site signal --
## replaces the earlier "correlation filter + random control" procedure with a
## threshold-free variance share.
## 目标：把"部位效应有多少是细胞组成造成的"从"计数 + 阈值"改成"方差占比 + 阈值无关"
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

suppressMessages({
  library(limma); library(variancePartition); library(BiocParallel)
})
args  <- commandArgs(trailingOnly = TRUE)
SMOKE <- length(args) > 0 && args[1] == "smoke"
NG    <- if (SMOKE && length(args) > 1) as.integer(args[2]) else NA
logf  <- file.path(BASE, if (SMOKE) "logs/varpart_smoke.txt" else "logs/varpart_out.txt")
con   <- file(logf, "w")
w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }

w("=========== variancePartition 方差分解 ===========")
w("R ", R.version.string, "  |  variancePartition ", as.character(packageVersion("variancePartition")),
  "  |  lme4 ", as.character(packageVersion("lme4")))
w("模式：", if (SMOKE) paste0("SMOKE (n=", NG, " genes)") else "FULL")
w("")

## ---------- 数据 ----------
p <- file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header = TRUE, check.names = FALSE, row.names = 1)
samp  <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
cpm0  <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep  <- rowSums(cpm0 >= 1) >= 4
dd    <- d0[keep, ]
Lc    <- log2(t(t(as.matrix(dd)) / colSums(as.matrix(dd)) * 1e6) + 1)
if (!is.na(NG)) Lc <- Lc[seq_len(min(NG, nrow(Lc))), , drop = FALSE]
w(sprintf("基因数 %d ；样本 %d（体内 %d / 体外 %d），供体 %d，部位 %s",
          nrow(Lc), length(samp), sum(state == "primary"), sum(state == "diff"),
          length(unique(donor)), paste(sort(unique(loc)), collapse = "/")))
w(sprintf("样本测序深度：min %.2f M，median %.1f M，max %.1f M",
          min(colSums(dd)/1e6), median(colSums(dd)/1e6), max(colSums(dd)/1e6)))
low <- names(which(colSums(dd) < 5e6))
w("低深度样本(<5M)：", if (length(low)) paste(low, collapse = ", ") else "无")
w("")

## ---------- 三轴分数（体内） ----------
mk <- read.delim(file.path(BASE, "results/tables/BMATID_marker_ensg.tsv"), stringsAsFactors = FALSE)
mk <- mk[mk$ensg != "NA", ]
AXES <- split(mk$ensg, mk$axis)[c("Hemato", "Bone", "Adipo")]
ens  <- sub("\\..*$", "", rownames(Lc))
idxmap <- tapply(seq_along(ens), ens, function(x) x)
score <- function(ids, cols) {
  rows <- unique(unlist(lapply(ids, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
  rows <- rows[!is.na(rows)]
  M <- Lc[rows, cols, drop = FALSE]
  sdv <- apply(M, 1, sd)
  drop <- names(sdv)[sdv == 0]
  M <- M[sdv > 0, , drop = FALSE]
  if (!nrow(M)) return(list(s = rep(NA_real_, length(cols)), drop = drop, n = 0))
  Z <- t(scale(t(M)))
  list(s = colMeans(Z, na.rm = TRUE), drop = drop, n = nrow(M))
}
pi <- which(state == "primary"); di <- which(state == "diff")
axP <- lapply(AXES, score, cols = pi)
axD <- lapply(AXES, score, cols = di)
S   <- sapply(axP, function(x) x$s)              # 12 x 3 体内
w("=== 三轴标记：体内可用数 / 体外被丢弃（零方差）===")
for (a in names(AXES)) {
  w(sprintf("  %-7s 体内 %2d / %2d 个可用；体外丢弃：%s", a, axP[[a]]$n, length(AXES[[a]]),
            if (length(axD[[a]]$drop)) paste(axD[[a]]$drop, collapse = ",") else "无"))
}
w("")
w("=== 轴分数相关（体内 12 样本）===")
w(paste(sprintf("%-8s", c("", colnames(S))), collapse = " "))
cm <- round(cor(S), 3)
for (i in seq_len(nrow(cm))) w(sprintf("%-8s %s", rownames(cm)[i], paste(sprintf("%6.3f", cm[i, ]), collapse = " ")))
w("")

## 每基因与三轴的 r（体内）
rc <- sapply(colnames(S), function(a) suppressWarnings(as.numeric(cor(t(Lc[, pi]), S[, a]))))
colnames(rc) <- colnames(S); rownames(rc) <- rownames(Lc)

## ---------- 建模元数据 ----------
meta <- data.frame(
  donor = factor(donor), localization = factor(loc, levels = c("epi","meta","scat")),
  state = factor(state), stringsAsFactors = FALSE)
rownames(meta) <- samp
meta$Hemato <- NA_real_; meta$Bone <- NA_real_; meta$Adipo <- NA_real_
meta$Hemato[pi] <- S[, "Hemato"]; meta$Bone[pi] <- S[, "Bone"]; meta$Adipo[pi] <- S[, "Adipo"]

## ---------- 后端 ----------
PBser <- SerialParam()
PBpar <- tryCatch(if (.Platform$OS.type == "windows") SnowParam(6, type = "SOCK", progressbar = FALSE)
                  else MulticoreParam(6), error = function(e) NULL)

## ---------- 主函数 ----------
runVP <- function(ex, md, form, tag, PB = PBser) {
  md <- droplevels(md)
  ## 协变量若全 NA 则跳过（否则 lmer 报 "data has 0"）
  covs <- setdiff(all.vars(form), c("donor", "localization"))
  if (length(covs)) {
    okrow <- Reduce(`&`, lapply(covs, function(v) is.finite(md[[v]])))
    if (!any(okrow)) { w("  !! ", tag, "：协变量 ", paste(covs, collapse = ","), " 全为 NA，跳过"); return(NULL) }
    if (!all(okrow)) { md <- md[okrow, , drop = FALSE]; ex <- ex[, okrow, drop = FALSE] }
  }
  t0 <- Sys.time()
  res <- tryCatch(
    suppressWarnings(fitExtractVarPartModel(ex, form, md, BPPARAM = PB, quiet = TRUE)),
    error = function(e) { w("  !! ", tag, " 报错: ", conditionMessage(e)); NULL })
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (is.null(res)) return(NULL)
  vp  <- as.data.frame(res)
  nres <- if ("Residuals" %in% colnames(vp)) vp[["Residuals"]] else rep(NA_real_, nrow(vp))
  vp  <- vp[, colnames(vp) != "Residuals", drop = FALSE]
  lcol <- grep("localization", colnames(vp), value = TRUE)[1]
  fin <- is.finite(vp[[lcol]])
  w(sprintf("  [%s] %s  用时 %.1fs  非有限 %d 基因  残差中位 %.1f%%",
            tag, paste(deparse(form), collapse=" "), dt, sum(!fin), 100 * median(nres, na.rm = TRUE)))
  attr(vp, "resid") <- nres
  vp
}

## ---------- M1：体内，非零方差基因（= 之前的 20073 基准集） ----------
w("=== M1 体内 primary：expr ~ (1|donor) + (1|localization) ===")
vpP <- runVP(Lc[, pi, drop = FALSE], meta[pi, , drop = FALSE],
             ~ (1|donor) + (1|localization), "in vivo")
if (!is.null(vpP) && SMOKE && !is.null(PBpar)) {
  w("  -- 并行后端烟测 --")
  vpPp <- runVP(Lc[, pi, drop = FALSE], meta[pi, , drop = FALSE],
                ~ (1|donor) + (1|localization), "in vivo/PAR", PB = PBpar)
  if (!is.null(vpPp)) {
    d <- abs(vpPp[[2]] - vpP[[2]])
    w(sprintf("     串行 vs 并行 最大差异 = %.2e", max(d, na.rm = TRUE)))
  }
}

summVP <- function(vp, tag) {
  if (is.null(vp)) return(invisible(NULL))
  lc <- grep("localization", colnames(vp), value = TRUE)[1]
  dn <- grep("donor", colnames(vp), value = TRUE)[1]
  v  <- vp[[lc]]; v <- v[is.finite(v)]
  d  <- vp[[dn]]; d <- d[is.finite(d)]
  w(sprintf("  %-26s n=%5d  %%loc 中位 %.1f%%（均值 %.1f%%，90分位 %.1f%%）  %%loc>20%%: %d (%.1f%%)  >50%%: %d  >80%%: %d",
            tag, length(v), 100*median(v), 100*mean(v), 100*quantile(v, .9, names = FALSE),
            sum(v > .2), 100*mean(v > .2), sum(v > .5), sum(v > .8)))
  w(sprintf("  %-26s      %%donor 中位 %.1f%%（均值 %.1f%%，90分位 %.1f%%）  %%donor>20%%: %d",
            "", 100*median(d), 100*mean(d), 100*quantile(d, .9, names = FALSE), sum(d > .2)))
  invisible(list(loc = v, donor = d))
}
w("")
w("=== 汇总：部位 / 供体 解释的方差占比 ===")
smP <- summVP(vpP, "体内 in vivo")

## ---------- M2：体外 ----------
w("")
w("=== M2 体外 diff：expr ~ (1|donor) + (1|localization) ===")
vpD <- runVP(Lc[, di, drop = FALSE], meta[di, , drop = FALSE],
             ~ (1|donor) + (1|localization), "in vitro")
smD <- summVP(vpD, "体外 in vitro")

## ---------- M3：体内加组成协变量 ----------
w("")
w("=== M3 体内 + 组成协变量（固定效应）：expr ~ (1|donor) + (1|localization) + Hemato + Bone ===")
vpP2 <- runVP(Lc[, pi, drop = FALSE], meta[pi, , drop = FALSE],
              ~ (1|donor) + (1|localization) + Hemato + Bone, "in vivo + composition")
smP2 <- summVP(vpP2, "体内 + 组成协变量")

## ---------- 共线性 ----------
w("")
w("=== canCorPairs：localization 与组成轴的共线性（体内 12 样本）===")
cc <- tryCatch(suppressWarnings(canCorPairs(~ localization + Hemato + Bone, meta[pi, , drop = FALSE])),
               error = function(e) { w("  canCorPairs 失败: ", conditionMessage(e)); NULL })
if (!is.null(cc)) { print(cc); writeLines(capture.output(print(round(as.matrix(cc), 3))), con) }

## ---------- 阈值无关的方差质量归因 ----------
w("")
w("=== 方差质量归因（阈值无关）：把每个基因的 %loc 当权重，看组成轴基因扛了多少 ===")
if (!is.null(vpP)) {
  lcv <- vpP[[grep("localization", colnames(vpP), value = TRUE)[1]]]
  names(lcv) <- rownames(vpP)
  tot <- sum(lcv, na.rm = TRUE)
  for (nm in list(c("Hemato", "|r(造血轴)|>0.7"), c("Bone", "|r(骨轴)|>0.7"))) {
    m <- abs(rc[names(lcv), nm[1]]) > 0.7
    w(sprintf("  %-18s 基因 %5d  方差质量占部位总质量 %.1f%%", nm[2], sum(m, na.rm = TRUE),
              100 * sum(lcv[m], na.rm = TRUE) / tot))
  }
  mboth <- abs(rc[names(lcv), "Hemato"]) > 0.7 | abs(rc[names(lcv), "Bone"]) > 0.7
  w(sprintf("  %-18s 基因 %5d  方差质量占比 %.1f%%", "任一组成轴>0.7", sum(mboth, na.rm = TRUE),
            100 * sum(lcv[mboth], na.rm = TRUE) / tot))
  w(sprintf("  %-18s 基因 %5d  方差质量占比 %.1f%%", "三轴都不相关(<=0.5)", sum(!mboth & apply(abs(rc[names(lcv), ]), 1, max, na.rm=TRUE) <= .5, na.rm=TRUE),
            100 * sum(lcv[!mboth & apply(abs(rc[names(lcv), ]), 1, max, na.rm=TRUE) <= .5], na.rm = TRUE) / tot))
  # 基因级关联
  w("")
  w("=== 基因级：%loc（体内）与 |r(轴)| 的 Spearman ===")
  for (a in colnames(rc)) {
    y <- abs(rc[names(lcv), a]); g <- is.finite(lcv) & is.finite(y)
    if (sum(g) < 30) { w(sprintf("  |r(%-7s)| : 有效基因仅 %d，跳过", a, sum(g))); next }
    r <- suppressWarnings(cor(lcv[g], y[g], method = "spearman"))
    w(sprintf("  |r(%-7s)| : rho = %+.3f  (n=%d)", a, r, sum(g)))
  }
  ## 按 %loc 分箱看组成轴关联系数
  brk <- c(0, .05, .1, .2, .4, 1)
  w("")
  w("=== 按 %loc 分箱：箱内 |r(Hemato)| / |r(Bone)| 中位数 ===")
  b <- cut(lcv, brk, include.lowest = TRUE, right = FALSE)
  for (lv in levels(b)) {
    idx <- which(b == lv)
    if (!length(idx)) next
    w(sprintf("  %%loc [%s)  基因 %5d   |r(Hemato)| 中位 %.3f   |r(Bone)| 中位 %.3f",
              lv, length(idx), median(abs(rc[names(lcv)[idx], "Hemato"]), na.rm = TRUE),
              median(abs(rc[names(lcv)[idx], "Bone"]), na.rm = TRUE)))
  }
}

## ---------- 敏感性：VST ----------
w("")
w("=== M4 敏感性：改用 DESeq2 VST（体内 / 体外）===")
ok <- suppressWarnings(requireNamespace("DESeq2", quietly = TRUE))
if (ok) {
  suppressMessages(library(DESeq2))
  cds <- DESeqDataSetFromMatrix(countData = round(as.matrix(dd)), colData = meta, design = ~ 1)
  vsd <- tryCatch(vst(cds, blind = TRUE), error = function(e) { w("  vst 失败：", conditionMessage(e)); NULL })
  if (!is.null(vsd)) {
    V <- assay(vsd); rownames(V) <- rownames(Lc)
    if (!is.na(NG)) V <- V[rownames(Lc), , drop = FALSE]
    vpPv <- runVP(V[, pi, drop = FALSE], meta[pi, , drop = FALSE], ~ (1|donor) + (1|localization), "VST in vivo")
    summVP(vpPv, "VST 体内")
    vpDv <- runVP(V[, di, drop = FALSE], meta[di, , drop = FALSE], ~ (1|donor) + (1|localization), "VST in vitro")
    summVP(vpDv, "VST 体外")
  }
}

## ---------- 敏感性：剔除低深度样本 ----------
w("")
w("=== M5 敏感性：剔除低深度样本后（体内）===")
if (length(low) && any(pi %in% which(samp %in% low))) {
  kp <- setdiff(pi, which(samp %in% low))
  vpPk <- runVP(Lc[, kp, drop = FALSE], meta[kp, , drop = FALSE], ~ (1|donor) + (1|localization), "in vivo drop low")
  summVP(vpPk, "体内剔低深度")
}

## ---------- M6 置换零分布（阈值无关的特异性对照，取代"等数量随机剔除"） ----------
w("")
w("=== M6 置换零分布：供体内打乱 localization 标签后重跑（阈值无关的对照）===")
set.seed(2026)
GSUB <- if (is.na(NG)) sample(nrow(Lc), 2000) else seq_len(nrow(Lc))
vpobs <- runVP(Lc[GSUB, pi, drop = FALSE], meta[pi, , drop = FALSE],
               ~ (1|donor) + (1|localization), "obs(subset)")
if (!is.null(vpobs)) {
  v <- vpobs[[grep("localization", colnames(vpobs), value = TRUE)[1]]]; v <- v[is.finite(v)]
  w(sprintf("     观测（同子集 n=%d）  %%loc 中位 %.1f%%  均值 %.1f%%  >20%%: %d (%.1f%%)",
            length(GSUB), 100*median(v), 100*mean(v), sum(v > .2), 100*mean(v > .2)))
}
for (k in 1:3) {
  md <- meta[pi, , drop = FALSE]
  md$localization <- unsplit(lapply(split(md$localization, md$donor), sample), md$donor)
  vpk <- runVP(Lc[GSUB, pi, drop = FALSE], md, ~ (1|donor) + (1|localization), sprintf("perm-%d", k))
  if (!is.null(vpk)) {
    v <- vpk[[grep("localization", colnames(vpk), value = TRUE)[1]]]; v <- v[is.finite(v)]
    w(sprintf("     perm-%d            %%loc 中位 %.1f%%  均值 %.1f%%  >20%%: %d (%.1f%%)",
              k, 100*median(v), 100*mean(v), sum(v > .2), 100*mean(v > .2)))
  }
}

## ---------- 导出 ----------
if (!SMOKE) {
  ex <- list()
  if (!is.null(vpP)) ex$invivo  <- vpP
  if (!is.null(vpD)) ex$invitro <- vpD
  if (!is.null(vpP2)) ex$invivo_comp <- vpP2
  for (nm in names(ex)) {
    df <- as.data.frame(ex[[nm]]); df$ensg <- rownames(df)
    write.csv(df, file.path(BASE, "results/tables", paste0("BMATID_varpart_", nm, ".csv")), row.names = FALSE)
  }
  saveRDS(list(rc = rc, S = S, meta = meta, low = low, axes_drop = lapply(axD, function(x) x$drop)),
          file.path(BASE, "results/tables/BMATID_varpart_inputs.rds"))
  w(""); w("导出：BMATID_varpart_{invivo,invitro,invivo_comp}.csv + BMATID_varpart_inputs.rds")
}
close(con)
cat("done\n")
