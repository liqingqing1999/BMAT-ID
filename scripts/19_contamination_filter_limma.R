
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
logf <- file.path(BASE, "logs", "clean_out2.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con); flush(con) }

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

pi <- which(state == "primary")
# PTPRC = ENSG00000081237
gi <- grep("^ENSG00000081237", rownames(Lc))
w("PTPRC rows found: ", length(gi))
cd45 <- as.numeric(Lc[gi[1], pi])
w("PTPRC(CD45) primary log2CPM: ", paste(round(cd45, 2), collapse=", "))
w("PTPRC range: ", round(min(cd45),2), " - ", round(max(cd45),2))

Lp <- Lc[, pi]
rc <- as.numeric(cor(t(Lp), cd45))
names(rc) <- rownames(Lp)
w("")
w("与 CD45 的相关分布：")
w("  NA 数（全零基因导致）: ", sum(is.na(rc)))
print(round(quantile(abs(rc), c(0,.25,.5,.75,.9,.95,.99,1), na.rm=TRUE), 3))
for (t in c(0.5, 0.6, 0.7, 0.8)) w(sprintf("  |r| > %.1f : %d genes (%.1f%%)", t, sum(abs(rc)>t, na.rm=TRUE), 100*mean(abs(rc)>t, na.rm=TRUE)))

# Every run() call is also recorded, so that the machine-readable ladder written at
# the end is exactly what this script computed -- 52_fig2_attribution_ladder.py reads
# that file instead of carrying the numbers in its own source.
LADDER <- list()
run <- function(gsel, tag) {
  dd2 <- dd[gsel, ]
  res <- c()
  for (st in c("primary","diff")) {
    idx <- which(state == st)
    m <- as.matrix(dd2[, idx]); m[is.na(m)] <- 0
    dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
    des <- model.matrix(~ dn + lc)
    v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
    fit <- eBayes(lmFit(v, des))
    ce <- grep("^lc", colnames(des), value=TRUE)
    tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
    res <- c(res, sum(tt$adj.P.Val < 0.05))
  }
  w(sprintf("  %-46s genes=%5d | FDR<0.05 primary=%5d  diff=%d", tag, sum(gsel), res[1], res[2]))
  LADDER[[length(LADDER) + 1]] <<- list(tag = tag, genes = sum(gsel), primary = res[1])
  invisible(res)
}

w("")
w("===== 剔除与 CD45 共变基因后，部位效应还剩多少（limma 同口径）=====")
run(rep(TRUE, nrow(dd)), "基准：全部基因")
for (t in c(0.8, 0.7, 0.6, 0.5)) run(!is.na(rc) & abs(rc) <= t, sprintf("剔除 |r(CD45)|>%.1f", t))

w("")
w("===== 对照：随机剔除等量基因（排除'剔除本身降低检出力'的假象）=====")
set.seed(42)          # fixed seed: the equal-number random control is reproducible
RND <- list()
for (t in c(0.7, 0.5)) {
  n_drop <- sum(!is.na(rc) & abs(rc) > t)
  rand <- sample(nrow(dd), n_drop)
  gsel <- rep(TRUE, nrow(dd)); gsel[rand] <- FALSE
  r <- run(gsel, sprintf("随机剔除 %d 个基因（= |r|>%.1f 的数量）", n_drop, t))
  # NB: plain <- here, not <<- : a top-level for() body already runs in the global
  # environment, and `RND[[i]] <<- v` would have to *read* RND from the parent first.
  RND[[length(RND) + 1]] <- list(removed = n_drop, primary = r[1])
}

w("")
w("===== 体外 diff 同样做（体外 CD45 归零，无法算相关，直接给全基因结果）=====")
run(rep(TRUE, nrow(dd)), "diff 全基因")

# ---- machine-readable ladder (Figure 2 reads this; relative to log2 CPM >= 1 filter) ----
# the two x/y vectors are paired for the figure: x = genes removed, y = genes with a
# site effect still present, both starting at the no-removal baseline.
cd45_removed <- c(0, vapply(c(0.8, 0.7, 0.6, 0.5),
                            function(t) sum(!is.na(rc) & abs(rc) > t), 0))
cd45_kept    <- c(LADDER[[1]]$primary, vapply(LADDER[2:5], function(x) x$primary, 0))
stopifnot(length(cd45_removed) == length(cd45_kept))
out <- file.path(BASE, "results", "tables", "BMATID_ladder_data.json")
writeLines(c("{",
  sprintf(' "baseline_genes": %d,', nrow(dd)),
  sprintf(' "baseline_site_dependent": %d,', LADDER[[1]]$primary),
  sprintf(' "removed_cd45": [%s],', paste(cd45_removed, collapse = ", ")),
  sprintf(' "site_dependent_after": [%s],', paste(cd45_kept, collapse = ", ")),
  sprintf(' "random_removed": [%s],', paste(vapply(RND, function(x) x$removed, 0), collapse = ", ")),
  sprintf(' "random_site_dependent": [%s]', paste(vapply(RND, function(x) x$primary, 0), collapse = ", ")),
  "}"), out)
w("")
w("saved: ", out)

close(con); cat("\n[log]", logf, "\n")
