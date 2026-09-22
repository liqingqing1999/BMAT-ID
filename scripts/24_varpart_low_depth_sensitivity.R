## Low-depth sensitivity: refit the in vitro model with the single low-depth
## library removed, and report the site share.
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
logf <- file.path(BASE, "logs/varpart_m5b_out.txt")
con  <- file(logf, "w")
w <- function(...) { s <- paste0(...); writeLines(s, con); flush(con) }

w("=== M5b：体外剔低深度样本 H410_epi_diff_adip 后的方差分解 ===")
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

meta <- data.frame(donor = factor(donor), localization = factor(loc, levels = c("epi","meta","scat")),
                   state = factor(state), stringsAsFactors = FALSE)
rownames(meta) <- samp
di  <- which(state == "diff")
low <- names(which(colSums(dd) < 5e6))
w("低深度样本：", paste(low, collapse = ", "))

summVP <- function(vp, tag) {
  lc <- grep("localization", colnames(vp), value = TRUE)[1]
  dn <- grep("donor", colnames(vp), value = TRUE)[1]
  v  <- vp[[lc]]; v <- v[is.finite(v)]; d <- vp[[dn]]; d <- d[is.finite(d)]
  w(sprintf("  %-24s n=%5d  部位方差 中位 %.1f%%（均值 %.1f%%，90分位 %.1f%%）  >20%%: %d (%.1f%%)  >50%%: %d  >80%%: %d",
            tag, length(v), 100*median(v), 100*mean(v), 100*quantile(v, .9, names = FALSE),
            sum(v > .2), 100*mean(v > .2), sum(v > .5), sum(v > .8)))
  w(sprintf("  %-24s        供体方差 中位 %.1f%%（均值 %.1f%%）", "", 100*median(d), 100*mean(d)))
}

run1 <- function(cols, tag) {
  md <- droplevels(meta[cols, , drop = FALSE])
  t0 <- Sys.time()
  res <- tryCatch(suppressWarnings(fitExtractVarPartModel(Lc[, cols, drop = FALSE],
              ~ (1|donor) + (1|localization), md, BPPARAM = SerialParam(), quiet = TRUE)),
           error = function(e) { w("  !! ", tag, " 报错: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  w(sprintf("  [%s] 用时 %.1fs  n_sample=%d", tag,
            as.numeric(difftime(Sys.time(), t0, units = "secs")), length(cols)))
  summVP(as.data.frame(res), tag)
}

w("")
w("--- 体外：剔低深度后 11 样本 ---")
run1(setdiff(di, which(samp %in% low)), "体外 剔低深度 11")
w("")
w("完成。")
close(con); cat("done\n")
