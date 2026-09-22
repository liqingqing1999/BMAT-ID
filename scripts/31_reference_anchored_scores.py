## Reference-anchored composition estimate for GSE291355 (human BMAd bulk)
## Reference = GSE169396 full atlas (human femoral head, 18,205 cells, 26 clusters).
## Goal: replace self-built marker axes with signatures derived from an INDEPENDENT
## single-cell dataset.
##
## NOTE (found while annotating): the femoral-head atlas is overwhelmingly immune;
## osteochondral cells form one cluster (178 cells) and no pure MSC/stromal cluster exists.
## The reference can therefore anchor the HAEMATOPOIETIC fraction robustly, but not the
## bone fraction. This limitation is itself reported.
# ---- portable project root -------------------------------------------------
# Nothing to edit: the root is the parent of this script's own directory.
# Override with the environment variable BMAT_ID_DIR if you keep the scripts
# somewhere else.  See README.md, section "Running the pipeline".
import os as _os
os = _os
BASE = _os.environ.get("BMAT_ID_DIR") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
for _d in ("logs", "results/tables", "results/figures"):
    _os.makedirs(_os.path.join(BASE, _d), exist_ok=True)
# ---------------------------------------------------------------------------

import io, os, gzip, json
import numpy as np
import pandas as pd
import scipy.optimize as opt

T   = os.path.join(BASE, "results", "tables") + "/"
RAW = os.path.join(BASE, "data", "raw") + "/"
L   = os.path.join(BASE, "logs", "deconv_out.txt")
log = io.open(L, "w", encoding="utf-8", buffering=1)
def say(*a):
    s = " ".join(str(x) for x in a); log.write(s + "\n"); print(s)

## ---------- 1. gene-id bridge ----------
ft = os.environ.get("BMAT_REF_FEATURES", "")
if not ft or not os.path.exists(ft):
    raise SystemExit("Set BMAT_REF_FEATURES to GSM5201883_S1_features.tsv.gz "
                     "(GSE169396 supplementary file); see README.md, section "
                     "\"External inputs\".")
ens2sym = {}
with gzip.open(ft, "rt", encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        if len(p) < 2: continue
        ens2sym[p[0].split(".")[0]] = p[1]
say("bridge:", len(ens2sym), "ensembl -> symbol")

## ---------- 2. reference pseudo-bulk ----------
pb = pd.read_csv(T + "BMATID_ref_pseudobulk_clusters.csv")
clusters = [c for c in pb.columns if c != "ensg"]
mat = pb.set_index("ensg")[clusters].astype(float)
say("pseudo-bulk:", mat.shape)

## ---------- 3. cluster -> super-type (assigned from the specificity table) ----------
lab = {
 "0":"immune","1":"immune","2":"immune","3":"immune","4":"immune","5":"immune","6":"immune",
 "7":"immune","8":"immune","9":"immune","12":"immune","14":"immune","17":"immune","18":"immune",
 "19":"immune","21":"immune","22":"immune","24":"immune","25":"immune",
 "10":"erythroid","13":"erythroid","15":"erythroid",
 "11":"skeletal",                    # osteo-chondral cluster (178 cells)
 "16":"vasculature","20":"vasculature",
 "23":"dropped",                     # no marker specificity
}
miss = [c for c in clusters if c not in lab]
say("clusters without a label:", miss)
lab = {c: v for c, v in lab.items() if c in clusters}

cby = pd.read_csv(T + "BMATID_ref_cluster_by_group.csv", index_col=0)
cby.index = cby.index.astype(str)
sup = {}
for c, l in lab.items():
    if l == "dropped": continue
    sup.setdefault(l, []).append(c)
say("\n== super-types ==")
for k, v in sorted(sup.items()):
    cells = int(cby.loc[[c for c in v if c in cby.index]].sum().sum())
    genes = int((mat[v] > 0).any(axis=1).sum())
    say("  %-12s %2d clusters, %6d cells, %5d genes detected" % (k, len(v), cells, genes))

agg = pd.DataFrame({k: mat[v].sum(axis=1) for k, v in sup.items()})
agg_cpm = agg / agg.sum(axis=0) * 1e6
say("signature matrix:", agg_cpm.shape)

## ---------- 4. bulk ----------
bulk = pd.read_csv(RAW + "GSE291355_counts.tsv.gz", sep="\t", index_col=0)
say("\nbulk:", bulk.shape)
say("bulk columns:", list(bulk.columns))
bulk.index = [str(i).split(".")[0] for i in bulk.index]
Bcpm = bulk.astype(float)
Bcpm = Bcpm / Bcpm.sum(axis=0) * 1e6

bmap = {}
for e in Bcpm.index:
    s = ens2sym.get(e)
    if s and s not in bmap: bmap[s] = e
say("bulk rows mapped to symbol:", len(bmap), "/", len(Bcpm))

common = sorted(set(agg_cpm.index) & set(bmap))
say("common symbols:", len(common))
Bl = pd.DataFrame({s: Bcpm.loc[bmap[s]] for s in common}).T
say("bulk subset:", Bl.shape)

state = {c: ("in vitro" if ("diff" in str(c).lower() or "vitro" in str(c).lower()) else "in vivo")
         for c in Bl.columns}
say("\nsample state:", json.dumps(state, ensure_ascii=False))
vivo = [c for c, v in state.items() if v == "in vivo"]
vitro = [c for c, v in state.items() if v == "in vitro"]
say("in vivo n=%d  in vitro n=%d" % (len(vivo), len(vitro)))

## ---------- 5. data-driven markers from the reference ----------
say("\n===== data-driven signatures from the reference =====")
lgR = np.log2(agg_cpm + 1)
spec = pd.DataFrame(index=agg_cpm.index)
for k in agg_cpm.columns:
    others = lgR[[c for c in agg_cpm.columns if c != k]].mean(axis=1)
    spec[k] = lgR[k] - others
sig_genes = {}
say("\n== top data-driven markers per super-type ==")
for k in agg_cpm.columns:
    cand = spec[k][(agg_cpm[k] > 20) & (spec[k] > 1)]
    cand = cand.sort_values(ascending=False)
    sig_genes[k] = list(cand.index[:200])
    say("  %-12s %4d genes;  e.g. %s" % (k, len(cand), ", ".join(sig_genes[k][:12])))
spec.to_csv(T + "BMATID_ref_supertype_specificity.csv", encoding="utf-8")

## ---------- 6. L1: reference-anchored scores on bulk ----------
say("\n===== L1: reference-anchored cell-type scores on bulk =====")
zB = np.log2(Bl + 1)
zB = zB.sub(zB.mean(axis=1), axis=0).div(zB.std(axis=1).replace(0, np.nan), axis=0)
L1 = pd.DataFrame(index=Bl.columns)
for k, gs in sig_genes.items():
    gs = [g for g in gs if g in zB.index]
    if len(gs) < 5: continue
    L1[k] = zB.loc[gs].mean(axis=0) * np.sqrt(len(gs))
say(L1.round(2).to_string())
L1.to_csv(T + "BMATID_deconv_L1_scores.csv", encoding="utf-8")

ptprc = bmap.get("PTPRC")
if ptprc:
    cd45 = np.log2(Bcpm.loc[ptprc] + 1)
    say("\n== correlation of each reference-anchored score with CD45 (PTPRC) log2CPM ==")
    for k in L1.columns:
        r = np.corrcoef(cd45.values.astype(float), L1[k].values.astype(float))[0, 1]
        say("   %-12s r = %+.3f" % (k, r))

say("\n== site-wise means of the reference-anchored immune score ==")
import collections
site = collections.defaultdict(list)
for c in Bl.columns:
    l = str(c).lower()
    s = "META" if "meta" in l else ("EPI" if "epi" in l else ("SCAT" if "scat" in l else "?"))
    site[(s, state[c])].append(c)
for k in sorted(site):
    say("   %-12s n=%d  immune score %.2f  (samples: %s)" % (
        str(k), len(site[k]), L1.loc[site[k], "immune"].mean(), ",".join(map(str, site[k]))))

## ---------- 7. L2: NNLS deconvolution ----------
say("\n===== L2: NNLS deconvolution =====")
sel = agg_cpm[(agg_cpm.max(axis=1) > 20)]
X = sel.loc[sel.index.intersection(Bl.index)]
Y = Bl.loc[X.index]
say("matrix:", X.shape, " (types: %s)" % ",".join(X.columns))
props = {}
resid = {}
for s in Y.columns:
    x, r = opt.nnls(X.values.astype(float), Y[s].values.astype(float))
    t = x.sum()
    props[s] = x / t if t > 0 else x
    resid[s] = r / (np.linalg.norm(Y[s].values.astype(float)) + 1e-9)
P = pd.DataFrame(props, index=X.columns)
say("\n== estimated composition (%) ==")
say((P * 100).round(1).to_string())
say("\nrelative residual:", {str(k): round(v, 3) for k, v in resid.items()})
P.to_csv(T + "BMATID_deconv_proportions.csv", encoding="utf-8")

say("\n== immune fraction (immune + erythroid) by site ==")
Pmix = P.copy()
Pmix.loc["myeloid_immune"] = P.loc[["immune", "erythroid"]].sum()
for k in sorted(site):
    say("   %-12s  immune %.3f  skeletal %.3f  vasculature %.3f  residual %.3f" % (
        str(k), Pmix.loc["myeloid_immune", site[k]].mean(),
        Pmix.loc["skeletal", site[k]].mean(), Pmix.loc["vasculature", site[k]].mean(),
        np.mean([resid[s] for s in site[k]])))

say("\nDONE")
log.close()
