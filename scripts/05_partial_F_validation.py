# -*- coding: utf-8 -*-
"""诊断 partial F 实现"""
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

import os, re, gzip, csv
import numpy as np
from scipy import stats

LOG = os.path.join(BASE, 'logs/diag2_out.txt')
buf = []
def w(s=''):
    buf.append(str(s)); print(s)

lines = gzip.open(os.path.join(BASE,'data/raw/GSE291355_counts.tsv.gz'),'rt',
                  encoding='utf-8', errors='replace').read().split('\n')
hdr = lines[0].split('\t')
samp = [h.strip() for h in hdr]
donor = [re.match(r'(H\d+)_', s).group(1) for s in samp]
loc   = [re.match(r'H\d+_([a-z]+)_', s).group(1) for s in samp]
state = [re.match(r'H\d+_[a-z]+_([a-z]+)_adip', s).group(1) for s in samp]
w('samples: %d ; state counts: %s' % (len(samp), {s: state.count(s) for s in set(state)}))

rows = list(csv.DictReader(open(os.path.join(BASE,'results/tables/BMATID_limma_allgenes.csv'), encoding='utf-8')))
ensg = [r['ensg'].split('.')[0] for r in rows]
kecp = set(ensg)
G = []; X = []
for _ln in lines[1:]:
    if not _ln.strip(): continue
    p = _ln.split('\t')
    g = p[0].split('.')[0]
    if g not in kecp: continue
    G.append(g); X.append([float(v) if v not in ('','NA') else 0.0 for v in p[1:]])
X = np.array(X)
w('X shape: %s' % (X.shape,))
w('order identical to limma? %s' % (G == ensg))

L = np.log2(X / X.sum(axis=0) * 1e6 + 1)
pi = [i for i, s in enumerate(state) if s == 'primary']
w('pi (primary idx): %s' % pi)
w('primary sample names: %s' % [samp[i] for i in pi])
w('primary donors: %s' % [donor[i] for i in pi])
w('primary locs: %s' % [loc[i] for i in pi])

ds = sorted(set(donor[i] for i in pi)); ls = ['meta','epi','scat']
w('donor levels: %s ; loc levels: %s' % (ds, ls))
Xf = []; Xr = []
for i in pi:
    b = [1.0] + [1.0 if donor[i]==d else 0.0 for d in ds[1:]]
    Xr.append(b)
    Xf.append(b + [1.0 if loc[i]==l else 0.0 for l in ls[1:]])
Xf = np.array(Xf); Xr = np.array(Xr)
w('Xf shape %s rank %d ; Xr shape %s rank %d' % (Xf.shape, np.linalg.matrix_rank(Xf), Xr.shape, np.linalg.matrix_rank(Xr)))

Yp = L[:, pi].T  # samples x genes
w('Yp shape: %s' % (Yp.shape,))

def rss(Xm, Y):
    B = np.linalg.lstsq(Xm, Y, rcond=None)[0]
    R = Y - Xm @ B
    return (R**2).sum(axis=0)

rf = rss(Xf, Yp); rr = rss(Xr, Yp)
w('rf: min %.4g med %.4g max %.4g ; #rf==0: %d' % (rf.min(), np.median(rf), rf.max(), int((rf<=0).sum())))
w('rr: min %.4g med %.4g max %.4g' % (rr.min(), np.median(rr), rr.max()))
df_f = Yp.shape[0] - np.linalg.matrix_rank(Xf)
df_r = Yp.shape[0] - np.linalg.matrix_rank(Xr)
F = ((rr - rf) / (df_r - df_f)) / (rf / df_f)
w('F: nan %d ; min %.3f med %.3f max %.3f' % (int(np.isnan(F).sum()), np.nanmin(F), np.nanmedian(F), np.nanmax(F)))
p = stats.f.sf(F, df_r-df_f, df_f)
w('p: nan %d ; min %.3g ; #p<0.05 = %d' % (int(np.isnan(p).sum()), np.nanmin(p), int((p<0.05).sum())))

# 与 limma 的 F 对比
lF = np.array([float(r['F_prim']) for r in rows])
lFDR = np.array([float(r['FDR_prim']) for r in rows])
m = ~np.isnan(F)
w('')
w('corr(my F, limma F) = %.4f  (n=%d)' % (np.corrcoef(F[m], lF[m])[0,1], m.sum()))
w('my F at limma FDR<0.05 genes: med %.2f ; limma F at same: med %.2f' % (np.median(F[m][lFDR<0.05]), np.median(lF[lFDR<0.05])))
w('my FDR<0.05 v.s. limma FDR<0.05 overlap: %d' % int(((p<0.05) & (lFDR<0.05)).sum()))

with open(LOG, 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(buf))
print('\n[log]', LOG)
