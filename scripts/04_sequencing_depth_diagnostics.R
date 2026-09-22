
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
logf <- file.path(BASE, "logs", "diag_out.txt")
con <- file(logf, "w")
w <- function(...) { s <- paste0(...); cat(s, "\n"); writeLines(s, con) }

p <- file.path(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
d0 <- read.delim(gzfile(p), header=TRUE, check.names=FALSE, row.names=1)
samp <- colnames(d0)
donor <- sub("^(H[0-9]+)_.*$", "\\1", samp)
loc   <- sub("^H[0-9]+_([a-z]+)_.*$", "\\1", samp)
state <- sub("^H[0-9]+_[a-z]+_([a-z]+)_adip$", "\\1", samp)

w("===== 1) library size（测序深度是否够）=====")
for (st in c("primary","diff")) {
  idx <- which(state == st)
  ls <- colSums(d0[, idx])
  w(sprintf("  %-8s median=%s  min=%s  max=%s", st,
            format(median(ls), big.mark=","), format(min(ls), big.mark=","), format(max(ls), big.mark=",")))
}

cpm <- t(t(as.matrix(d0)) / colSums(as.matrix(d0)) * 1e6)
L <- log2(cpm + 1)

w(""); w("===== 2) 组内平均 Pearson 相关（log2CPM，全基因）=====")
for (st in c("primary","diff")) {
  idx <- which(state == st)
  C <- cor(L[, idx]); diag(C) <- NA
  w(sprintf("  %-8s overall within-group mean r = %.4f", st, mean(C, na.rm=TRUE)))
  for (g in c("epi","meta","scat")) {
    ii <- which(state == st & loc == g)
    if (length(ii) > 1) {
      cc <- cor(L[, ii]); diag(cc) <- NA
      w(sprintf("        %-5s within mean r = %.4f", g, mean(cc, na.rm=TRUE)))
    }
  }
  # 跨部位
  for (gg in combn(c("epi","meta","scat"), 2, simplify=FALSE)) {
    a <- which(state == st & loc == gg[1]); b <- which(state == st & loc == gg[2])
    w(sprintf("        %s vs %s  mean r = %.4f", gg[1], gg[2], mean(cor(L[, a], L[, b]))))
  }
}

w(""); w("===== 3) p 值分布诊断（localization F 检验，看是否只是功效不足）=====")
keep <- rowSums(cpm >= 1) >= 4
dd <- d0[keep, ]
for (st in c("primary","diff")) {
  idx <- which(state == st)
  m <- as.matrix(dd[, idx]); m[is.na(m)] <- 0
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  v <- voom(m, des, lib.size=colSums(m), normalize.method="quantile")
  fit <- eBayes(lmFit(v, des))
  ce <- grep("^lc", colnames(des), value=TRUE)
  tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
  pv <- tt$P.Value
  w(sprintf("  %-8s  p<0.05 未校正 = %5d (%.1f%%)  | p<0.01 = %5d | FDR<0.05 = %5d",
            st, sum(pv < 0.05), 100*mean(pv < 0.05), sum(pv < 0.01), sum(tt$adj.P.Val < 0.05)))
  # p 值十分位均匀性（如果全部 p>0.9 说明模型退化）
  qs <- quantile(pv, probs=seq(0,1,0.1))
  w("           p 值分位: "); w(paste(sprintf("            %.2f", qs), collapse="\n"))
  # moderated F
  w(sprintf("           moderated F: median=%.2f  max=%.2f", median(tt$F), max(tt$F)))
}

w(""); w("===== 4) 反推原文口径：'specifically expressed' 的 exclusive 计数 =====")
w("  原文 Fig4B(in vitro)= EPI 2 / META 11 / SCAT 29 ; Fig4D(in vivo)= EPI 234 / META 918 / SCAT 1290")
for (T in c(0, 0.5, 1, 5, 10)) {
  w(sprintf("  ---- 阈值 cpm >= %g ----", T))
  for (st in c("primary","diff")) {
    idx <- which(state == st); grp <- loc[idx]
    out <- c()
    for (g in c("epi","meta","scat")) {
      a <- idx[grp == g]; b <- idx[grp != g]
      na_ <- rowSums(cpm[, a, drop=FALSE] >= T); nb <- rowSums(cpm[, b, drop=FALSE] >= T)
      out <- c(out, sum(na_ == length(a) & nb == 0))
    }
    w(sprintf("    %-8s EPI=%5d  META=%5d  SCAT=%5d", st, out[1], out[2], out[3]))
  }
}

w(""); w("===== 5) 稳健性：改用 limma-trend（log2CPM + trend）复核 diff 是否仍为 0 =====")
for (st in c("primary","diff")) {
  idx <- which(state == st)
  m <- as.matrix(dd[, idx]); m[is.na(m)] <- 0
  l2 <- normalizeBetweenArrays(log2(t(t(m)/colSums(m))*1e6 + 1), method="quantile")
  dn <- factor(donor[idx]); lc <- factor(loc[idx], levels=c("meta","epi","scat"))
  des <- model.matrix(~ dn + lc)
  fit <- eBayes(lmFit(l2, des), trend=TRUE)
  ce <- grep("^lc", colnames(des), value=TRUE)
  tt <- topTable(fit, coef=ce, number=Inf, sort.by="none")
  w(sprintf("  %-8s limma-trend: FDR<0.05 = %d ; p<0.05 = %d", st, sum(tt$adj.P.Val<0.05), sum(tt$P.Value<0.05)))
}

close(con); cat("\n[log]", logf, "\n")
