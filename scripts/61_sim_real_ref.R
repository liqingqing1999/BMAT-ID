## =====================================================================
## real-data reference on the simulation's own gene pools
##
## Why this exists.  60_sim_validation.R compares the simulated site shares
## against the real ones, but the real side of that comparison was, for the
## adjusted model, only available over the full 20,041-gene universe.  The
## simulation draws a 400-gene pool per replicate that is deliberately enriched
## for the marker genes, so its expected naive site share is not the all-gene
## 0.171 but something higher.  This script removes that mismatch: it redraws
## the SAME pool for each replicate (the pool is the first draw after
## set.seed(), so the stream is reproducible exactly) and fits both real models
## on exactly those genes.
##
## Outputs results/tables/BMATID_sim_validation_realref.csv and a log.
## Read-only with respect to the manuscript.
## =====================================================================
suppressMessages({
  library(variancePartition); library(BiocParallel); library(lme4); library(parallel)
})

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

NG    <- 400L
RREP  <- 30L
NC    <- 6L
SEED  <- 20260929L

dir.create(file.path(BASE, "logs"), showWarnings = FALSE)
dir.create(file.path(BASE, "results/tables"), showWarnings = FALSE)
con <- file(file.path(BASE, "logs/sim_real_ref.txt"), "w")
w <- function(...) { writeLines(paste0(...), con); flush(con) }

w("=== real-data reference on the simulation's own gene pools ===")
w("R ", R.version.string)

## ------------------------------------------------------- substrate (same as sim)
d0 <- read.delim(gzfile(file.path(BASE, "data/raw/GSE291355_counts.tsv.gz")),
                 header = TRUE, check.names = FALSE, row.names = 1)
samp   <- colnames(d0)
donor  <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc    <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state  <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
pi     <- which(state == "primary")
C      <- as.matrix(d0[, pi, drop = FALSE]); storage.mode(C) <- "double"
lib    <- colSums(C)
donP   <- donor[pi]
locChr <- loc[pi]
locP   <- factor(locChr, levels = c("epi", "meta", "scat"))
nS     <- length(pi)
dons   <- sort(unique(donP))

cpmR <- t(t(C) / lib * 1e6)
LRr  <- log2(cpmR + 1)
cpmAll <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep   <- rowSums(cpmAll >= 1) >= 4
LR     <- LRr[keep, , drop = FALSE]
ng     <- nrow(LR)

mk     <- read.delim(file.path(BASE, "results/tables/BMATID_marker_ensg.tsv"), stringsAsFactors = FALSE)
mk     <- mk[mk$ensg != "NA", ]
ens    <- sub("\\..*$", "", rownames(LR))
idxmap <- tapply(seq_along(ens), ens, function(x) x)
axisScore <- function(ids) {
  rows <- unique(unlist(lapply(ids, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
  M <- LR[rows, , drop = FALSE]
  s <- apply(M, 1, sd); M <- M[s > 0, , drop = FALSE]
  colMeans(t(scale(t(M))))
}
hemIds <- mk$ensg[mk$axis == "Hemato"]; bonIds <- mk$ensg[mk$axis == "Bone"]
cHv <- as.numeric(scale(axisScore(hemIds)))
cBv <- as.numeric(scale(axisScore(bonIds)))
mkIdx <- unique(unlist(lapply(c(hemIds, bonIds),
                              function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))

w(sprintf("genes %d | markers %d | libraries %d", ng, length(mkIdx), nS))

## ------------------------------------------------- full-universe reference (cache, not refit)
## The all-gene decomposition already exists: the main pipeline wrote it to
## results/tables/BMATID_varpart_invivo.csv (naive) and BMATID_varpart_invivo_comp.csv (adjusted), and
## the manuscript publishes those values (17.1% -> 4.2%).  Refitting 20,073 genes here cost
## ~20 min of wall time for a number that was already known, so only the per-pool part below,
## which IS new, is computed.
mdFull <- data.frame(donor = factor(donP),
                     localization = factor(locChr, levels = c("epi", "meta", "scat")),
                     Hemato = cHv, Bone = cBv, stringsAsFactors = FALSE)
rownames(mdFull) <- colnames(LR)
fitShare <- function(M, form, md) {
  vp <- suppressWarnings(fitExtractVarPartModel(M, form, md, BPPARAM = SerialParam(), quiet = TRUE))
  v <- as.data.frame(vp)
  v[[grep("localization", colnames(v), value = TRUE)[1]]]
}
cN <- read.csv(file.path(BASE, "results/tables/BMATID_varpart_invivo.csv"))
cA <- read.csv(file.path(BASE, "results/tables/BMATID_varpart_invivo_comp.csv"))
stopifnot("localization" %in% colnames(cN), "localization" %in% colnames(cA))
w(sprintf("full universe (cached, %d genes): naive %.4f | adjusted %.4f | %%zero %.1f%% -> %.1f%%",
          nrow(cN), mean(cN$localization, na.rm = TRUE), mean(cA$localization, na.rm = TRUE),
          100 * mean(cN$localization == 0, na.rm = TRUE),
          100 * mean(cA$localization == 0, na.rm = TRUE)))

## ------------------------------------------------- per-replicate matched pools
REF <- function(r) {
  set.seed(SEED + 100000 + 1000 * r)
  pool <- if (length(mkIdx) >= NG) mkIdx[seq_len(NG)] else
          c(mkIdx, sample(setdiff(seq_len(ng), mkIdx), NG - length(mkIdx)))
  M  <- LR[pool, , drop = FALSE]
  MD <- mdFull
  shN <- fitShare(M, ~ (1|donor) + (1|localization), MD)
  shA <- fitShare(M, ~ (1|donor) + (1|localization) + Hemato + Bone, MD)
  data.frame(rep = r,
             real_naive_mean = mean(shN, na.rm = TRUE),
             real_adj_mean   = mean(shA, na.rm = TRUE),
             real_naive_med  = median(shN, na.rm = TRUE),
             real_adj_med    = median(shA, na.rm = TRUE),
             real_naive_zero = mean(shN == 0, na.rm = TRUE),
             real_adj_zero   = mean(shA == 0, na.rm = TRUE))
}

cl <- makeCluster(NC)
clusterEvalQ(cl, suppressMessages({ library(variancePartition); library(BiocParallel); library(lme4) }))
clusterExport(cl, c("LR", "mdFull", "mkIdx", "ng", "NG", "SEED", "fitShare"),
              envir = environment())

res <- do.call(rbind, parLapply(cl, seq_len(RREP), REF))
stopCluster(cl)

w("")
w("--- real data, the simulation's own 400-gene pools ---")
w(sprintf("naive    mean over replicates %.4f | median of replicate means %.4f",
          mean(res$real_naive_mean), median(res$real_naive_mean)))
w(sprintf("adjusted mean over replicates %.4f | median of replicate means %.4f",
          mean(res$real_adj_mean), median(res$real_adj_mean)))
w(sprintf("range across replicates: naive %.4f-%.4f | adjusted %.4f-%.4f",
          min(res$real_naive_mean), max(res$real_naive_mean),
          min(res$real_adj_mean), max(res$real_adj_mean)))

out <- file.path(BASE, "results/tables/BMATID_sim_validation_realref.csv")
write.csv(res, out, row.names = FALSE)
w("")
w("exported: results/tables/BMATID_sim_validation_realref.csv")
w("DONE")
close(con)
cat("real-ref DONE\n")
