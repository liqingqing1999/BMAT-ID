
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
logf <- file.path(BASE, "logs/residuals_out.txt")
con <- file(logf, "w"); w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con); flush(con) }

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
gi <- grep("^ENSG00000081237", rownames(Lc))
cd45 <- as.numeric(Lc[gi[1], pi])
Lp <- Lc[, pi]
rc <- suppressWarnings(as.numeric(cor(t(Lp), cd45)))
names(rc) <- rownames(Lp)

# limma on primary, same model as export
idx <- pi
m <- as.matrix(dd[, idx]); m[is.na(m)] <- 0
dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
des <- model.matrix(~ dn + lc)
v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
fit <- eBayes(lmFit(v, des))
ce <- grep("^lc", colnames(des), value=TRUE)
tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
cm <- makeContrasts(EM = lcepi, ES = lcepi - lcscat, MS = lcscat, levels = des)
f2 <- eBayes(contrasts.fit(lmFit(v, des), cm))
lf <- sapply(colnames(cm), function(cc) topTable(f2, coef=cc, number=Inf, sort.by="none")$logFC)

res <- data.frame(
  ensg     = rownames(tt),
  r_cd45   = rc[rownames(tt)],
  F_prim   = tt$F,
  FDR_prim = tt$adj.P.Val,
  P_prim   = tt$P.Value,
  lfc_EM   = lf[, "EM"], lfc_ES = lf[, "ES"], lfc_MS = lf[, "MS"],
  EPI_p    = rowMeans(Lc[, idx[loc[idx] == "epi"],  drop = FALSE]),
  META_p   = rowMeans(Lc[, idx[loc[idx] == "meta"], drop = FALSE]),
  SCAT_p   = rowMeans(Lc[, idx[loc[idx] == "scat"], drop = FALSE]),
  stringsAsFactors = FALSE)

outdir <- file.path(BASE, "results/tables")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
f <- file.path(outdir, "BMATID_primary_cd45_annotated.csv")
write.csv(res, f, row.names = FALSE)

sig   <- res$FDR_prim < 0.05
noncon <- (is.na(res$r_cd45) | abs(res$r_cd45) <= 0.7)
w("total kept genes                : ", nrow(res))
w("primary FDR<0.05                : ", sum(sig))
w("primary FDR<0.05 & |r|<=0.7     : ", sum(sig & noncon))
w("primary FDR<0.01 & |r|<=0.7     : ", sum(res$FDR_prim < 0.01 & noncon, na.rm=TRUE))
w("primary FDR<0.05 & |r|<=0.5     : ", sum(sig & (is.na(res$r_cd45) | abs(res$r_cd45) <= 0.5)))
w("saved: ", f)

# also export just the residual significant gene list for enrichment/downstream
sink(file.path(BASE, "results/tables/BMATID_residual_site_genes.tsv"))
r <- res[sig & noncon, c("ensg","r_cd45","F_prim","FDR_prim","lfc_EM","lfc_ES","lfc_MS","EPI_p","META_p","SCAT_p")]
r <- r[order(-r$F_prim), ]
write.table(r, row.names = FALSE, sep = "\t", quote = FALSE)
sink()
w("residual gene list rows: ", nrow(r))
close(con)
