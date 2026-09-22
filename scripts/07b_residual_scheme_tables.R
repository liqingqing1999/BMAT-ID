# ---- residual-site gene schemes: two filtering calibers, one partial-F fit --------
# Writes the two tables that `09_*`, `41_*` and `42_*` read:
#   results/tables/BMATID_residual494_schemeA.tsv    (filter first, then fit)
#   results/tables/BMATID_residual1248_schemeB.tsv   (fit all genes, then filter)
# Input: data/raw/GSE291355_counts.tsv.gz and nothing else.
#
# The two calibers differ in where the CD45 co-variation filter is applied:
#   A  the limma prior is re-estimated on the non-contaminant genes only
#      (filter, then fit) -> the conservative set used in the manuscript;
#   B  all genes are fitted, then the non-contaminant ones are kept
#      (fit, then filter) -> the permissive sensitivity set.
# The intersection of A and B is reported in the log; in this dataset A is a subset
# of B, so the intersection has the same members as A and no third table is written.
#
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
logf <- file.path(BASE, "logs", "residual_schemes_out.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con); flush(con) }

p <- file.path(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp  <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)
cpm0  <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
keep  <- rowSums(cpm0 >= 1) >= 4
dd    <- d0[keep, ]
w("genes retained (CPM >= 1 in >= 4 samples): ", nrow(dd))
Lc <- log2(t(t(as.matrix(dd)) / colSums(as.matrix(dd)) * 1e6) + 1)

pi   <- which(state == "primary")
gi   <- grep("^ENSG00000081237", rownames(Lc))          # PTPRC / CD45
cd45 <- as.numeric(Lc[gi[1], pi])
rc   <- suppressWarnings(as.numeric(cor(t(Lc[, pi]), cd45))); names(rc) <- rownames(Lc)
w("is.na(r) (all-zero genes): ", sum(is.na(rc)))

fitlimma <- function(gsel) {
  dd2 <- dd[gsel, ]
  idx <- pi
  m <- as.matrix(dd2[, idx]); m[is.na(m)] <- 0
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  ce <- grep("^lc", colnames(des), value=TRUE)
  tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
  cm <- makeContrasts(EM = lcepi, ES = lcepi - lcscat, MS = lcscat, levels = des)
  f2 <- eBayes(contrasts.fit(lmFit(v, des), cm))
  lf <- sapply(colnames(cm), function(cc) topTable(f2, coef=cc, number=Inf, sort.by="none")$logFC)
  list(ensg = rownames(tt), F = tt$F, FDR = tt$adj.P.Val, P = tt$P.Value,
       lfc_EM = lf[,"EM"], lfc_ES = lf[,"ES"], lfc_MS = lf[,"MS"])
}

noncon <- (is.na(rc) | abs(rc) <= 0.7)
w("non-contaminant genes (fit first): ", sum(noncon), " / ", length(noncon))

# caliber A: filter first, then fit (the limma prior is re-estimated on the clean set)
A <- fitlimma(noncon)
sigA <- A$FDR < 0.05
w("caliber A  filter-then-fit : FDR<0.05 = ", sum(sigA), " / ", length(sigA))

# caliber B: fit all genes, then filter
B <- fitlimma(rep(TRUE, nrow(dd)))
sigB <- (B$FDR < 0.05) & noncon
w("caliber B  fit-then-filter: FDR<0.05 & |r|<=0.7 = ", sum(sigB))

# intersection (reported only; see the header note)
inter <- sigA & sigB[match(A$ensg, B$ensg)]
w("A n B = ", sum(inter), "   (A is a subset of B: ", all(inter == sigA), ")")

Lcs <- Lc
site_mean <- function(gsel, s) rowMeans(Lcs[gsel, pi[loc[pi] == s], drop=FALSE])

mk <- function(o, sig, tag) {
  d <- data.frame(ensg = o$ensg, r_cd45 = rc[o$ensg], F_prim = o$F, FDR_prim = o$FDR,
                  lfc_EM = o$lfc_EM, lfc_ES = o$lfc_ES, lfc_MS = o$lfc_MS)
  d$EPI  <- site_mean(o$ensg, "epi")
  d$META <- site_mean(o$ensg, "meta")
  d$SCAT <- site_mean(o$ensg, "scat")
  d <- d[sig, ]
  mx <- pmax(d$EPI, d$META, d$SCAT); mn <- pmin(d$EPI, d$META, d$SCAT)
  d$range <- mx - mn
  d$top_site <- c("EPI","META","SCAT")[apply(cbind(d$EPI, d$META, d$SCAT), 1, which.max)]
  d <- d[order(-d$F_prim), ]
  f <- file.path(BASE, sprintf("results/tables/%s.tsv", tag))
  write.table(d, f, row.names=FALSE, sep="\t", quote=FALSE)
  w("saved: ", f, "  rows=", nrow(d))
  w("    top_site distribution: ", paste(names(table(d$top_site)), table(d$top_site), sep="=", collapse="  "))
  invisible(d)
}

dA <- mk(A, sigA, "BMATID_residual494_schemeA")
dB <- mk(B, sigB, "BMATID_residual1248_schemeB")
close(con)
cat("\n[log]", logf, "\n")
