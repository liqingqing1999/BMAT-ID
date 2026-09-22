## M7: variance decomposition with REFERENCE-ANCHORED composition covariates
## (immune / skeletal signatures derived from the independent GSE169396 atlas)
## Compare with the model that uses the self-built marker axes instead.
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

suppressMessages({ library(variancePartition); library(BiocParallel) })
logf <- file.path(BASE, "logs/deconv_varpart_out.txt")
con  <- file(logf, "w"); w <- function(...) { writeLines(paste0(...), con); flush(con) }

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
pi    <- which(state == "primary")
w("genes ", nrow(Lc), " | primary samples ", length(pi))

L1 <- read.csv(file.path(BASE, "results/tables/BMATID_deconv_L1_scores.csv"),
               row.names = 1, check.names = FALSE)
w("L1 cols: ", paste(colnames(L1), collapse = ","))
L1 <- L1[samp, , drop = FALSE]
stopifnot(identical(rownames(L1), samp))

meta <- data.frame(donor = factor(donor),
                   localization = factor(loc, levels = c("epi", "meta", "scat")),
                   Immune_ref = L1$immune, Skeletal_ref = L1$skeletal)
rownames(meta) <- samp
w("reference covariates (12 primary):")
w(paste(sprintf("%-26s immune=%7.2f skeletal=%7.2f", samp[pi], meta$Immune_ref[pi], meta$Skeletal_ref[pi]),
        collapse = "\n"))
w("")
w("cor(Immune_ref, Skeletal_ref) = ", round(cor(meta$Immune_ref[pi], meta$Skeletal_ref[pi]), 3))

run <- function(form, tag) {
  t0 <- Sys.time()
  res <- tryCatch(suppressWarnings(fitExtractVarPartModel(Lc[, pi, drop = FALSE], form,
                                                          meta[pi, , drop = FALSE],
                                                          BPPARAM = SerialParam(), quiet = TRUE)),
                  error = function(e) { w("  !! ", tag, " err: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(NULL)
  vp <- as.data.frame(res)
  nres <- if ("Residuals" %in% colnames(vp)) vp[["Residuals"]] else rep(NA_real_, nrow(vp))
  vp2 <- vp[, colnames(vp) != "Residuals", drop = FALSE]
  lc <- grep("localization", colnames(vp2), value = TRUE)[1]
  x  <- vp2[[lc]]
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  w(sprintf("\n== %s : %s   (%.0fs)", tag, paste(deparse(form), collapse = ""), dt))
  w(sprintf("   localization: mean %.2f%%  median %.2f%%  >20%%: %d (%.1f%%)  >50%%: %d  ==0: %d (%.1f%%)",
            100 * mean(x, na.rm = TRUE), 100 * median(x, na.rm = TRUE),
            sum(x > .2, na.rm = TRUE), 100 * sum(x > .2, na.rm = TRUE) / sum(is.finite(x)),
            sum(x > .5, na.rm = TRUE),
            sum(x == 0, na.rm = TRUE), 100 * sum(x == 0, na.rm = TRUE) / sum(is.finite(x))))
  for (cn in setdiff(colnames(vp2), lc))
    w(sprintf("   %-16s mean %.2f%%  median %.2f%%", cn,
              100 * mean(vp2[[cn]], na.rm = TRUE), 100 * median(vp2[[cn]], na.rm = TRUE)))
  w(sprintf("   residual median %.1f%%", 100 * median(nres, na.rm = TRUE)))
  vp2
}

## reference: M1 (no covariates) re-run for a matched comparison
run(~ (1 | donor) + (1 | localization), "M1ref (no covariates)")
## M7: reference-anchored composition covariates
vp7 <- run(~ (1 | donor) + (1 | localization) + Immune_ref + Skeletal_ref, "M7 (reference-anchored covariates)")

if (!is.null(vp7)) {
  out <- data.frame(ensg = rownames(vp7), vp7, check.names = FALSE)
  write.csv(out, file.path(BASE, "results/tables/BMATID_varpart_M7_reference.csv"), row.names = FALSE)
  w("\nwrote BMATID_varpart_M7_reference.csv")
}
w("\nDONE")
close(con)
cat("done\n")
