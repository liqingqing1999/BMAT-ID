## M6b：2x2 对照 —— 有/无组成协变量 x 真/打乱部位标签（同一 n=2000 子集）
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

suppressMessages({ library(limma); library(variancePartition); library(BiocParallel) })
logf <- file.path(BASE, "logs/varpart_m6b_out.txt")
con  <- file(logf, "w")
w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }

w("=== M6b：2x2 对照（有/无组成协变量 × 真/打乱部位标签）===")
w("目的：严格判定 M3（组成校正）后的残余部位方差是否已等于随机标签")
w("")

p  <- file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header = TRUE, check.names = FALSE, row.names = 1)
samp  <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
cpm0  <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep  <- rowSums(cpm0 >= 1) >= 4
dd    <- d0[keep, ]
Lc    <- log2(t(t(as.matrix(dd)) / colSums(as.matrix(dd)) * 1e6) + 1)

## 三轴分数（体内）
mk <- read.delim(file.path(BASE, "results/tables/BMATID_marker_ensg.tsv"), stringsAsFactors = FALSE)
mk <- mk[mk$ensg != "NA", ]
AXES <- split(mk$ensg, mk$axis)[c("Hemato", "Bone", "Adipo")]
ens <- sub("\\..*$", "", rownames(Lc))
idxmap <- tapply(seq_along(ens), ens, function(x) x)
pi <- which(state == "primary")
score <- function(ids, cols) {
  rows <- unique(unlist(lapply(ids, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
  M <- Lc[rows, cols, drop = FALSE]
  M <- M[apply(M, 1, sd) > 0, , drop = FALSE]
  colMeans(t(scale(t(M))))
}
axP <- lapply(AXES, score, cols = pi)

meta <- data.frame(donor = factor(donor), localization = factor(loc, levels = c("epi","meta","scat")),
                   state = factor(state), stringsAsFactors = FALSE)
rownames(meta) <- samp
mp <- meta[pi, , drop = FALSE]; mp$Hemato <- axP$Hemato; mp$Bone <- axP$Bone

set.seed(2026)
GSUB <- sample(nrow(Lc), 2000)
ex <- Lc[GSUB, pi, drop = FALSE]

getloc <- function(vp) { v <- as.data.frame(vp)[[grep("localization", colnames(as.data.frame(vp)), value = TRUE)[1]]]; v[is.finite(v)] }

fit <- function(form, md, tag) {
  t0 <- Sys.time()
  r <- tryCatch(suppressWarnings(fitExtractVarPartModel(ex, form, md, BPPARAM = SerialParam(), quiet = TRUE)),
                error = function(e) { w("  !! ", tag, " ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  v <- getloc(r)
  w(sprintf("  %-34s 用时 %4.0fs  部位方差：均值 %.1f%%  中位 %.1f%%  >20%%: %d (%.1f%%)",
            tag, as.numeric(difftime(Sys.time(), t0, units = "secs")),
            100*mean(v), 100*median(v), sum(v > .2), 100*mean(v > .2)))
  invisible(v)
}

w("--- 无组成协变量：expr ~ (1|donor) + (1|localization) ---")
fit(~ (1|donor) + (1|localization), mp, "真标签 (obs)")
for (k in 1:3) {
  md <- mp; md$localization <- unsplit(lapply(split(md$localization, md$donor), sample), md$donor)
  fit(~ (1|donor) + (1|localization), md, sprintf("打乱标签 (perm-%d)", k))
}

w("")
w("--- 有组成协变量：+ Hemato + Bone ---")
fit(~ (1|donor) + (1|localization) + Hemato + Bone, mp, "真标签 (obs + comp)")
set.seed(2027)
for (k in 1:3) {
  md <- mp; md$localization <- unsplit(lapply(split(md$localization, md$donor), sample), md$donor)
  fit(~ (1|donor) + (1|localization) + Hemato + Bone, md, sprintf("打乱标签 (perm-%d + comp)", k))
}

w("")
w("读法：若 'obs + comp' ≈ 'perm + comp'，则组成校正后部位项已等于噪声。")
close(con); cat("done\n")
