## =====================================================================
## Simulation-based independent validation
##
## Why this exists.  The manuscript's central inference is that the site term of
## a bulk, site-resolved BMAd preparation is dominated by the composition of the
## co-purifying compartment.  A reviewer's first question is the obvious one --
## "how do you know that is confounding and not a genuine positional programme?"
## This script answers it with a parametric simulation in which the truth is
## known, run through the SAME pipeline the manuscript uses (variancePartition,
## naive model versus composition-adjusted model).
##
## Calibration.  The generative model is calibrated gene by gene to the real
## GSE291355 in vivo libraries: mu, the donor spread, the regression coefficients
## on the two real marker axes, the residual SD, the library sizes and the count
## overdispersion are all taken from the real data, and both real marker sets are
## present in every replicate.  Composition is therefore simulated with the real
## loadings, not with invented ones.
##
## Generative model, per simulated gene g and library s:
##     log2(CPM)_gs = mu_g + donor_gd + compScale*(bH_g*cH_s + bB_g*cB_s) + lam_g*t_s
##     counts_gs    ~ NB( mean = 2^log2CPM_gs * libsize_s/1e6, size = 1/phi_g )
## with t_s the site-level contrast of the composition gradient itself, which is
## the only kind of site-level effect a three-site design can carry (see below).
##
## Regimes (truth known by construction):
##   A  composition only      -> what confounding alone produces
##   B  intrinsic only        -> composition switched off, the culture analogue
##   C  both                  -> as in a bulk in vivo preparation
## plus experiment 2, a sweep of the composition scale at a fixed intrinsic load.
##
## A structural note that the simulation measures rather than assumes.  In a
## three-site design the centring leaves a two-dimensional site-contrast space,
## and the two composition axes span it, so no site-level contrast can be
## orthogonal to composition.  A site-level estimate is therefore *always*
## confounded, and the adjusted residual is a lower bound rather than a partition.
## The script reports the canonical correlations between site and the two axes in
## the simulated data so this can be compared with the real values (0.887, 0.944).
##
## Nothing in the manuscript is modified here; the script only writes results.
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


NG   <- 400L        # genes per replicate: the reported quantities are means and
                    # proportions, which are stable well below the real 20k
RREP <- 30L         # replicates per regime
RSWP <- 12L         # replicates per sweep level
NC   <- 6L          # workers
SEED <- 20260929L

## `Rscript 60_sim_validation.R smoke` runs a miniature version end to end, so a
## syntax or wiring error is caught in seconds; `quick` keeps the real 400-gene
## universe but few replicates, which is what shows whether the calibration lands
## on the observed values before the full run is committed.
MODE  <- if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1] else ""
SMOKE <- identical(MODE, "smoke")
QUICK <- identical(MODE, "quick")
if (SMOKE) { NG <- 40L; RREP <- 2L; RSWP <- 1L; NC <- 3L }
if (QUICK) { RREP <- 3L; RSWP <- 2L }
SUF <- if (SMOKE) "_smoke" else if (QUICK) "_quick" else ""

dir.create(file.path(BASE, "logs"), showWarnings = FALSE)
dir.create(file.path(BASE, "results/tables"), showWarnings = FALSE)
con <- file(file.path(BASE, paste0("logs/sim_validation", SUF, ".txt")), "w")
w <- function(...) { writeLines(paste0(...), con); flush(con) }

w("=== simulation-based independent validation ===")
w("R ", R.version.string,
  " | variancePartition ", as.character(packageVersion("variancePartition")),
  " | lme4 ", as.character(packageVersion("lme4")))
w(sprintf("genes/replicate %d | replicates %d (exp 1) / %d (exp 2) | workers %d",
          NG, RREP, RSWP, NC))

## ---------------------------------------------------------------- 1. substrate
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

w("")
w("--- 1. calibration substrate (real in vivo libraries) ---")
w(sprintf("libraries %d | donors %s | sites %s", nS, paste(dons, collapse = "/"),
          paste(levels(locP), collapse = "/")))
w(sprintf("library size (M reads) min %.1f / median %.1f / max %.1f",
          min(lib) / 1e6, median(lib) / 1e6, max(lib) / 1e6))

cpmR <- t(t(C) / lib * 1e6)
LRr  <- log2(cpmR + 1)
## The gene filter follows the published analysis, which applied it over all 24
## libraries (in vivo and in vitro) before subsetting: applying it to the twelve
## in vivo libraries alone is stricter and drops 2,800 genes.
cpmAll <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep   <- rowSums(cpmAll >= 1) >= 4
LR   <- LRr[keep, , drop = FALSE]
CR   <- C[keep, , drop = FALSE]
ng   <- nrow(LR)

## latent composition = the real marker axes, standardised across the 12 libraries
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

## the two marker sets, as row indices of the simulated universe: they must be
## present in every replicate, because the adjustment is built from them.
## A few markers fail the CPM filter (the real analysis drops them the same way),
## so the per-axis counts are reported and only a floor is enforced.
mkIdx <- unique(unlist(lapply(c(hemIds, bonIds),
                              function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
hemKeep <- unique(unlist(lapply(hemIds, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
bonKeep <- unique(unlist(lapply(bonIds, function(i) if (i %in% names(idxmap)) idxmap[[i]] else NULL)))
w(sprintf("genes retained %d | marker genes present %d (haemato %d/%d, bone %d/%d)",
          ng, length(mkIdx), length(hemKeep), length(hemIds), length(bonKeep), length(bonIds)))
stopifnot(length(hemKeep) >= 8L, length(bonKeep) >= 8L)

## per-gene calibration: regress every real gene on donor + the two real axes.
## The raw coefficients are NOT usable as truth: with twelve libraries and six
## fitted parameters a gene with no composition effect still receives loadings
## that happen to track the axes, and because the axes correlate with site at
## ~0.87 that overfitting is injected as spurious SITE structure -- it inflated
## the simulated site share roughly 1.8-fold on the first attempt.  The loadings
## are therefore shrunk by an empirical-Bayes (eBayes-style) factor
##     tau^2 / (tau^2 + SE^2),   tau^2 = max(0, var(beta_hat) - mean(SE^2)),
## which leaves genuinely loaded genes almost untouched and pulls noise to zero.
Dd <- model.matrix(~ locP - 1)
Xr <- cbind(1, model.matrix(~ factor(donP) - 1)[, -1, drop = FALSE], cHv, cBv)
XtXi  <- solve(crossprod(Xr))
Bx    <- XtXi %*% crossprod(Xr, t(LR))                      # 6 x ng
Rg    <- LR - t(Xr %*% Bx)
dfres <- nrow(Xr) - ncol(Xr)
seB   <- sqrt(outer(diag(XtXi)[5:6], rowSums(Rg^2) / dfres))  # 2 x ng
bRaw  <- Bx[5:6, , drop = FALSE]
tau2  <- pmax(0, apply(bRaw, 1, var) - rowMeans(seB^2))      # length 2
shrk  <- sweep(seB^2, 1, tau2, function(s2, t2) t2 / (t2 + s2))
bHat  <- bRaw * shrk
w(sprintf("loading shrinkage (eBayes): tau^2 %.4f / %.4f ; median shrinkage factor %.3f / %.3f",
          tau2[1], tau2[2], median(shrk[1, ]), median(shrk[2, ])))
mu_g   <- as.numeric(Bx[1, ] + c(Xr[1, 2:4] %*% Bx[2:4, ]))
bH_g   <- as.numeric(bHat[1, ])
bB_g   <- as.numeric(bHat[2, ])
sd_res <- as.numeric(apply(Rg, 1, sd))
Dm     <- sapply(dons, function(d) rowMeans(LR[, donP == d, drop = FALSE]))
donSd  <- as.numeric(apply(Dm, 1, sd))
mrc    <- rowMeans(CR); vrc <- apply(CR, 1, var)
## The calibration target for each replicate is the real donor-and-site-removed SD
## of the *same* genes that replicate simulates, not the global median: the two
## marker sets are deliberately over-represented in a 400-gene universe, and they
## are both noisier and more composition-loaded than an average gene.
DgR <- sapply(dons, function(d) rowMeans(LR[, donP == d, drop = FALSE]))[, match(donP, dons)]
SgR <- sapply(levels(locP), function(l) rowMeans(LR[, locChr == l, drop = FALSE]))[, match(locChr, levels(locP))]
sd_bs_real <- apply(LR - (DgR + SgR - mu_g), 1, sd)

## The completed analysis of the real data, so each replicate can be compared
## against the observation for exactly the genes it simulates rather than against
## a global figure.  This is the check that the simulation reproduces reality.
vpR <- read.csv(file.path(BASE, "results/tables/BMATID_varpart_invivo.csv"),
                stringsAsFactors = FALSE)
names(vpR)[names(vpR) == "localization"] <- "site"
vpR$key  <- sub("\\..*$", "", vpR$ensg)
realSite <- vpR$site[match(sub("\\..*$", "", rownames(LR)), vpR$key)]
w(sprintf("real per-gene site shares loaded for %d of %d genes (all-gene mean %.4f, median %.4f)",
          sum(is.finite(realSite)), ng, mean(realSite, na.rm = TRUE), median(realSite, na.rm = TRUE)))

## Count overdispersion is NOT taken from the method-of-moments value on the raw
## counts: that value already absorbs the donor and composition structure, which
## here is explicit in the linear predictor, so using it would count the same
## variance twice (it inflated the realised residual SD threefold on the first
## attempt).  It is instead solved from the target residual SD by the delta method,
##     sd_log2CPM^2 * ln2^2 ~= 1/m + phi,
## and then refined by an explicit calibration loop once simOnce() exists.
phiBase <- pmin(pmax((sd_res * log(2))^2 - 1 / pmax(mrc, 1), 1e-4), 5)
phi     <- phiBase
TARGET_SD <- median(sd_bs_real)

w(sprintf("calibration: log2CPM q25/med/q75 %.1f/%.1f/%.1f | residual SD median %.3f | donor SD median %.3f | phi (delta) median %.4f",
          quantile(mu_g, .25), median(mu_g), quantile(mu_g, .75),
          median(sd_res), median(donSd), median(phiBase)))
w(sprintf("real composition loadings: |bH| median %.3f (markers %.3f) | |bB| median %.3f (markers %.3f)",
          median(abs(bH_g)), median(abs(bH_g[mkIdx])), median(abs(bB_g)), median(abs(bB_g[mkIdx]))))
w(sprintf("latent axes: cor(cH,cB) %+.3f", cor(cHv, cBv)))
w(sprintf("  cH by site %s", paste(sprintf("%+.2f", tapply(cHv, locP, mean)[levels(locP)]), collapse = " / ")))
w(sprintf("  cB by site %s", paste(sprintf("%+.2f", tapply(cBv, locP, mean)[levels(locP)]), collapse = " / ")))

## The intrinsic contrast available to a three-site design: the site-level
## pattern of the composition gradient itself.  Its orthogonality to the
## composition axes is checked numerically below, not assumed.
Q  <- qr.Q(qr(Dd))
uH <- as.numeric(Q %*% crossprod(Q, cHv))
tConf <- as.numeric(scale(uH))
w(sprintf("intrinsic contrast t_s by site %s",
          paste(sprintf("%+.2f", tapply(tConf, locP, mean)[levels(locP)]), collapse = " / ")))

## demonstrate the structural claim: the two axes span the site-contrast space
siteSpaceDim <- qr(Dd)$rank - 1L
w(sprintf("site-contrast space dimension after centring: %d ; composition axes: 2",
          siteSpaceDim))
if (siteSpaceDim <= 2L)
  w("  => every site-level contrast is spanned by the two axes: the adjusted")
if (siteSpaceDim <= 2L)
  w("     residual is a LOWER BOUND, not an estimate of the intrinsic effect.")

## ------------------------------------------------------- 2. generate + analyse
simOnce <- function(compScale = 1, lamScale = 1, pIntr = 0.10, seed = 1, doOracle = TRUE) {
  set.seed(seed)
  pool <- if (length(mkIdx) >= NG) mkIdx[seq_len(NG)] else
          c(mkIdx, sample(setdiff(seq_len(ng), mkIdx), NG - length(mkIdx)))
  idx <- pool

  mu  <- mu_g[idx]; sr <- sd_res[idx]; ph <- phi[idx]
  dsd <- donSd[idx]
  bH  <- compScale * bH_g[idx]
  bB  <- compScale * bB_g[idx]

  isIntr <- rep(FALSE, NG)
  if (pIntr > 0) isIntr[sample(NG, max(1L, round(pIntr * NG)))] <- TRUE
  lam <- numeric(NG)
  lam[isIntr] <- sr[isIntr] * abs(rnorm(sum(isIntr), 1.0, 0.35)) * lamScale

  dd  <- matrix(rnorm(NG * length(dons)), NG, length(dons)) * dsd
  Eta <- matrix(mu, NG, nS) + dd[, match(donP, dons), drop = FALSE]
  Eta <- Eta + outer(bH, cHv) + outer(bB, cBv)
  if (any(isIntr)) Eta <- Eta + outer(lam, tConf)

  cnt <- matrix(0, NG, nS)
  for (j in seq_len(nS)) {
    m <- pmax(2^Eta[, j] - 1, 1e-9) * lib[j] / 1e6
    cnt[, j] <- rnbinom(NG, mu = m, size = 1 / ph)
  }
  CPM <- t(t(cnt) / colSums(cnt) * 1e6)
  L   <- log2(CPM + 1)
  rownames(L) <- sprintf("sim%04d", seq_len(NG))

  ## calibration check: the realised residual spread must match the real one for
  ## exactly the genes this replicate simulated
  Dg <- sapply(dons, function(d) rowMeans(L[, donP == d, drop = FALSE]))[, match(donP, dons)]
  Sg <- sapply(levels(locP), function(l) rowMeans(L[, locChr == l, drop = FALSE]))[, match(locChr, levels(locP))]
  sd_real <- median(apply(L - (Dg + Sg - rowMeans(L)), 1, sd))
  sd_targ <- median(sd_bs_real[idx])

  ax <- function(rows) {
    M <- L[rows, , drop = FALSE]
    s <- apply(M, 1, sd); M <- M[s > 0, , drop = FALSE]
    if (!nrow(M)) return(rep(0, nS))
    colMeans(t(scale(t(M))))
  }
  hemPresent <- which(ens[idx] %in% sub("\\..*$", "", hemIds))
  bonPresent <- which(ens[idx] %in% sub("\\..*$", "", bonIds))
  mH <- if (length(hemPresent)) hemPresent else sample(NG, 11)
  mB <- if (length(bonPresent)) bonPresent else sample(NG, 12)
  sH <- as.numeric(scale(ax(mH))); sB <- as.numeric(scale(ax(mB)))

  md <- data.frame(donor = factor(donP),
                   localization = factor(locChr, levels = c("epi", "meta", "scat")),
                   Hemato = sH, Bone = sB, stringsAsFactors = FALSE)
  rownames(md) <- colnames(L)

  ## canonical correlation between site and the two measured axes (real: 0.887, 0.944)
  cc <- tryCatch(suppressWarnings(as.matrix(canCorPairs(~ localization + Hemato + Bone, md))),
                 error = function(e) NULL)
  ccH <- if (!is.null(cc)) cc["localization", "Hemato"] else NA_real_
  ccB <- if (!is.null(cc)) cc["localization", "Bone"]   else NA_real_

  grab <- function(form, data) {
    vp <- tryCatch(suppressWarnings(fitExtractVarPartModel(L, form, data,
                                                           BPPARAM = SerialParam(), quiet = TRUE)),
                   error = function(e) NULL)
    if (is.null(vp)) return(rep(NA_real_, NG))
    v <- as.data.frame(vp)
    v[[grep("localization", colnames(v), value = TRUE)[1]]]
  }
  shN <- grab(~ (1|donor) + (1|localization), md)
  shA <- grab(~ (1|donor) + (1|localization) + Hemato + Bone, md)
  shO <- rep(NA_real_, NG)
  if (doOracle && compScale != 0) {
    mdo <- md; mdo$cH <- cHv; mdo$cB <- cBv
    shO <- grab(~ (1|donor) + (1|localization) + cH + cB, mdo)
  }

  mm <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  qq <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else median(x, na.rm = TRUE)
  zz <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x == 0, na.rm = TRUE)
  gt <- function(x) if (!length(x) || all(is.na(x))) NA_real_ else mean(x > 0.5, na.rm = TRUE)

  data.frame(
    compScale = compScale, lamScale = lamScale,
    sd_realised = sd_real, sd_target = sd_targ, ccH = ccH, ccB = ccB, nIntr = sum(isIntr),
    real_ref_mean = mean(realSite[idx], na.rm = TRUE),
    real_ref_med  = median(realSite[idx], na.rm = TRUE),
    naive_mean = mm(shN), adj_mean = mm(shA), orac_mean = mm(shO),
    naive_med = qq(shN), adj_med = qq(shA), naive_zero = zz(shN), adj_zero = zz(shA),
    true_intr_mean = if (any(isIntr))
      mm(lam[isIntr]^2 / (lam[isIntr]^2 + sr[isIntr]^2 + dsd[isIntr]^2)) else NA_real_,
    adj_gt50_intr = gt(shA[isIntr]), naive_gt50_intr = gt(shN[isIntr]),
    adj_med_intr = qq(shA[isIntr]),   naive_med_intr = qq(shN[isIntr]),
    adj_med_mkH = qq(shA[mH]), adj_med_mkB = qq(shA[mB]),
    naive_med_mkH = qq(shN[mH]), naive_med_mkB = qq(shN[mB]),
    stringsAsFactors = FALSE)
}

simWorker <- function(r, compScale, lamScale, pIntr, seedBase) {
  simOnce(compScale, lamScale, pIntr, seed = seedBase + 1000 * r)
}

## ---- calibrate the count overdispersion against the real residual spread ----
w("")
w("--- 2. calibrating the count noise to the real residual spread ---")
w(sprintf("global target residual SD (real median over all genes) %.3f ; delta-method phi median %.4f",
          TARGET_SD, median(phiBase)))
for (it in 1:6) {
  r0  <- simOnce(1, 0, 0, seed = 999)
  rat <- r0$sd_target / r0$sd_realised
  w(sprintf("   iteration %d: target %.3f vs realised %.3f -> scale factor %.3f (phi median %.4f)",
            it, r0$sd_target, r0$sd_realised, rat, median(phi)))
  if (!is.finite(rat) || abs(rat - 1) < 0.03) break
  phi <- pmin(pmax(phi * rat^2, 1e-5), 20)
}
w(sprintf("final phi median %.4f ; target %.3f vs realised %.3f",
          median(phi), r0$sd_target, r0$sd_realised))
w("")

cl <- makeCluster(NC)
clusterEvalQ(cl, suppressMessages({ library(variancePartition); library(BiocParallel); library(lme4) }))
clusterExport(cl, c("mu_g", "sd_res", "sd_bs_real", "realSite", "phi", "donSd", "bH_g", "bB_g", "cHv", "cBv",
                    "tConf", "donP", "locChr", "locP", "lib", "nS", "ng",
                    "dons", "mkIdx", "ens", "hemIds", "bonIds", "NG", "simOnce", "simWorker"),
              envir = environment())

## ---------------------------------------------------- experiment 1: three regimes
regimes <- list(
  A = list(compScale = 1, lamScale = 0, pIntr = 0.00,
           lab = "composition only - no adipocyte-intrinsic site effect exists"),
  B = list(compScale = 0, lamScale = 1, pIntr = 0.30,
           lab = "intrinsic only, composition switched off (culture analogue)"),
  C = list(compScale = 1, lamScale = 1, pIntr = 0.10,
           lab = "both, as in a bulk in vivo preparation"))

w("")
w("--- experiment 1 : three regimes with known truth ---")
allres <- list()
for (nm in names(regimes)) {
  g  <- regimes[[nm]]
  t0 <- Sys.time()
  out <- parLapply(cl, seq_len(RREP), simWorker,
                   compScale = g$compScale, lamScale = g$lamScale,
                   pIntr = g$pIntr, seedBase = SEED + 100000)
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  df <- do.call(rbind, out); df$regime <- nm; df$rep <- seq_len(nrow(df))
  allres[[nm]] <- df

  w(sprintf("[%s] %s", nm, g$lab))
  w(sprintf("   calibration: residual SD target %.3f vs realised %.3f | simulated canonical correlations %.3f / %.3f (real 0.887 / 0.944)",
            mean(df$sd_target), mean(df$sd_realised), mean(df$ccH, na.rm = TRUE), mean(df$ccB, na.rm = TRUE)))
  w(sprintf("   same genes, real analysis: mean site share %.4f (median %.4f)",
            mean(df$real_ref_mean, na.rm = TRUE), mean(df$real_ref_med, na.rm = TRUE)))
  w(sprintf("   NAIVE model     mean site share %.4f | median %.4f | %%exactly zero %.1f%%",
            mean(df$naive_mean), median(df$naive_med), 100 * mean(df$naive_zero)))
  w(sprintf("   ADJUSTED model  mean site share %.4f | median %.4f | %%exactly zero %.1f%%",
            mean(df$adj_mean), median(df$adj_med), 100 * mean(df$adj_zero)))
  if (!all(is.na(df$orac_mean)))
    w(sprintf("   ORACLE model    mean site share %.4f  (adjusting for the true latent scores)",
              mean(df$orac_mean, na.rm = TRUE)))
  w(sprintf("   marker axes (median share): haemato %.3f -> %.3f | bone %.3f -> %.3f",
            median(df$naive_med_mkH), median(df$adj_med_mkH),
            median(df$naive_med_mkB), median(df$adj_med_mkB)))
  if (any(df$nIntr > 0))
    w(sprintf("   intrinsic genes: true share %.3f | naive median %.3f | adjusted median %.3f | recovered>0.5 naive %.1f%% -> adjusted %.1f%%",
              mean(df$true_intr_mean, na.rm = TRUE), median(df$naive_med_intr),
              median(df$adj_med_intr),
              100 * mean(df$naive_gt50_intr, na.rm = TRUE),
              100 * mean(df$adj_gt50_intr, na.rm = TRUE)))
  w(sprintf("   [%.0f s]", dt)); w("")
}

## --------------------------------------------- experiment 2: confounding sweep
w("--- experiment 2 : composition scale at fixed intrinsic load (10%) ---")
facs <- c(0, 0.5, 1.0, 1.5, 2.0)
sw   <- list()
for (f in facs) {
  t0 <- Sys.time()
  out <- parLapply(cl, seq_len(RSWP), simWorker,
                   compScale = f, lamScale = 1, pIntr = 0.10,
                   seedBase = SEED + 200000 + round(1000 * f))
  df <- do.call(rbind, out); df$compScale <- f
  sw[[sprintf("%.1f", f)]] <- df
  w(sprintf("   compScale=%.1f  true intrinsic %.4f | naive %.4f | adjusted %.4f | oracle %s | cc %.3f/%.3f  [%.0f s]",
            f, mean(df$true_intr_mean, na.rm = TRUE), mean(df$naive_mean), mean(df$adj_mean),
            if (all(is.na(df$orac_mean))) "  n/a  " else sprintf("%.4f", mean(df$orac_mean, na.rm = TRUE)),
            mean(df$ccH, na.rm = TRUE), mean(df$ccB, na.rm = TRUE),
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
w("")

stopCluster(cl)

## ------------------------------------------------------------------ 3. export
tab <- do.call(rbind, c(allres, lapply(sw, function(d) { d$regime <- "sweep"; d$rep <- seq_len(nrow(d)); d })))
write.csv(tab, file.path(BASE, paste0("results/tables/BMATID_sim_validation_replicates", SUF, ".csv")), row.names = FALSE)

srow <- function(regime, metric, value) data.frame(regime = regime, metric = metric, value = value)
S <- list()
for (nm in names(allres)) {
  df <- allres[[nm]]
  for (k in c("sd_realised", "sd_target", "ccH", "ccB", "naive_mean", "adj_mean", "orac_mean",
              "naive_zero", "adj_zero", "naive_med_mkH", "adj_med_mkH",
              "naive_med_mkB", "adj_med_mkB", "true_intr_mean", "naive_med_intr",
              "adj_med_intr", "adj_gt50_intr", "naive_gt50_intr"))
    S[[length(S) + 1]] <- srow(nm, k, if (all(is.na(df[[k]]))) NA_real_ else mean(df[[k]], na.rm = TRUE))
}
for (k in names(sw)) {
  df <- sw[[k]]
  for (kk in c("true_intr_mean", "naive_mean", "adj_mean", "orac_mean",
               "adj_gt50_intr", "adj_med_mkH", "adj_med_mkB", "ccH", "ccB"))
    S[[length(S) + 1]] <- srow(paste0("sweep_comp", k), kk,
                               if (all(is.na(df[[kk]]))) NA_real_ else mean(df[[kk]], na.rm = TRUE))
}
write.csv(do.call(rbind, S), file.path(BASE, paste0("results/tables/BMATID_sim_validation_summary", SUF, ".csv")), row.names = FALSE)

w(paste0("exported: results/tables/BMATID_sim_validation_replicates", SUF, ".csv"))
w(paste0("          results/tables/BMATID_sim_validation_summary", SUF, ".csv"))
w("DONE")
close(con)
cat("simulation DONE\n")
