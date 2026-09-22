# -*- coding: utf-8 -*-
"""Functional-module decomposition of the site effect, plus the haematopoietic-
   contamination correction. Writes results/tables/BMATID_module_retention.csv."""
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

import os, re, gzip, json, ssl, urllib.request, urllib.parse
import numpy as np

LOG = os.path.join(BASE, 'logs/modules_out.txt')
os.makedirs(os.path.dirname(LOG), exist_ok=True)
buf = []
def w(s=''):
    buf.append(str(s)); print(s)

ctx = ssl.create_default_context(); ctx.check_hostname=False; ctx.verify_mode=ssl.CERT_NONE

# ---------- 1) 读 limma 结果 ----------
import csv
f = os.path.join(BASE, 'results/tables/BMATID_limma_allgenes.csv')
rows = list(csv.DictReader(open(f, encoding='utf-8')))
w('limma rows: %d' % len(rows))
ensg = [r['ensg'].split('.')[0] for r in rows]
for r in rows:
    for k in ['F_prim','FDR_prim','F_diff','FDR_diff','maxlfc_prim','maxlfc_diff','EPI_p','META_p','SCAT_p']:
        r[k] = float(r[k])
idx_by_ensg = {e: i for i, e in enumerate(ensg)}

# ---------- 2) 模块基因 ----------
MODULES = {
 'A_骨/矿化':      ['SPP1','RUNX2','COL1A1','COL1A2','ALPL','SP7','BGLAP','IBSP','MGP','MMP13','MMP14',
                    'DLX5','BARX1','SPARC','COL2A1','SOX9','ACAN','BMP2','BMP4','SOST','PHEX','ENPP1','ANKH'],
 'B_造血龛':       ['CXCL12','VCAM1','KITLG','LEPR','ANGPT1','JAG1','OSM','IL7','TPO'],
 'C_造血/浆细胞(污染指纹)': ['PTPRC','CD14','CD68','LYZ','C1QA','C1QB','C1QC','IGHG1','IGHG2','IGHG3',
                    'IGKC','IGLC1','IGLC2','JCHAIN','CD3E','CD3D','CD19','MS4A1','CD79A','MZB1','IGLL5','SDC1'],
 'D_成脂核心':     ['PPARG','CEBPA','CEBPD','FABP4','PLIN1','PLIN4','ADIPOQ','LEP','LPL','CFD','GPAM',
                    'SCD','FASN','DGAT1','DGAT2','CIDEC','PNPLA2','FABP5'],
 'E_内皮':         ['VWF','PECAM1','CDH5','CLDN5'],
 'F_MSC/干性':     ['ENG','THY1','NT5E','MCAM','ALCAM','PDPN','LEPR'],
}

# ---------- 3) symbol -> ENSG ----------
allsym = sorted({g for v in MODULES.values() for g in v})
sym2ensg = {}
try:
    u = 'https://mygene.info/v3/query'
    data = urllib.parse.urlencode({'q': ','.join(allsym), 'scopes': 'symbol',
                                   'fields': 'ensembl.gene,symbol', 'species': 'human'}).encode()
    req = urllib.request.Request(u, data=data, headers={'Content-Type':'application/x-www-form-urlencoded'})
    res = json.loads(urllib.request.urlopen(req, context=ctx, timeout=90).read())
    hits = res if isinstance(res, list) else res.get('hits', [])
    for h in hits:
        q = h.get('query')
        e = h.get('ensembl')
        if not q or h.get('notfound') or not e:
            continue
        if isinstance(e, list):
            e = e[0]
        gid = e.get('gene') if isinstance(e, dict) else None
        if gid:
            sym2ensg[q.upper()] = gid
except Exception as ex:
    w('mygene ERROR: %s' % ex)
w('mapped symbols: %d / %d' % (len(sym2ensg), len(allsym)))

# ---------- 4) 逐模块汇总 ----------
w('')
w('=' * 118)
w('模块层面的"部位效应"：体内 vs 体外')
w('=' * 118)
w('%-26s %5s %9s %9s %10s %10s %8s %8s' % ('module','n','medF_prim','medF_diff','medF_ratio','medlfc_p','medlfc_d','ratio'))
w('-' * 118)
module_genes = {}
MODRET = []          # machine-readable retention table -> results/tables/BMATID_module_retention.csv
for mod, syms in MODULES.items():
    ii = []
    for s in syms:
        e = sym2ensg.get(s.upper())
        if e and e in idx_by_ensg:
            ii.append((s, idx_by_ensg[e]))
    if not ii:
        w('%-26s %5d  (none found)' % (mod, 0)); continue
    Fp = np.array([rows[i]['F_prim'] for _, i in ii])
    Fd = np.array([rows[i]['F_diff'] for _, i in ii])
    Lp = np.array([rows[i]['maxlfc_prim'] for _, i in ii])
    Ld = np.array([rows[i]['maxlfc_diff'] for _, i in ii])
    sp = np.sum([rows[i]['FDR_prim'] < 0.05 for _, i in ii])
    sd = np.sum([rows[i]['FDR_diff'] < 0.05 for _, i in ii])
    module_genes[mod] = ii
    MODRET.append({'module': mod, 'n_genes': len(ii),
                   'median_F_in_vivo': np.median(Fp), 'median_F_in_vitro': np.median(Fd),
                   'retention_pct': 100 * np.median(Fd) / max(np.median(Fp), 1e-9)})
    w('%-26s %5d %9.2f %9.2f %9.1f%% %10.2f %8.2f %7.1f%%   (sig prim=%d diff=%d)' % (
        mod, len(ii), np.median(Fp), np.median(Fd),
        100*np.median(Fd)/max(np.median(Fp),1e-9),
        np.median(Lp), np.median(Ld), 100*np.median(Ld)/max(np.median(Lp),1e-9),
        sp, sd))

w('')
w('--- 各模块基因明细（按体内 F 降序）---')
for mod, ii in module_genes.items():
    w('')
    w('### ' + mod)
    w('%-10s %9s %9s %9s %9s %9s %9s %9s %8s' % ('sym','EPI_p','META_p','SCAT_p','F_prim','F_diff','lfc_p','lfc_d','FDR_p'))
    for s, i in sorted(ii, key=lambda x: -rows[x[1]]['F_prim']):
        r = rows[i]
        w('%-10s %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %9.2f %8.2g' % (
            s, r['EPI_p'], r['META_p'], r['SCAT_p'], r['F_prim'], r['F_diff'],
            r['maxlfc_prim'], r['maxlfc_diff'], r['FDR_prim']))

# ---------- 4b) machine-readable retention table (Figure 3 reads this) ----------
csvret = os.path.join(BASE, 'results/tables/BMATID_module_retention.csv')
with open(csvret, 'w', newline='', encoding='utf-8') as fh:
    cw = csv.writer(fh)
    cw.writerow(['module', 'n_genes', 'median_F_in_vivo', 'median_F_in_vitro', 'retention_pct'])
    for r in MODRET:
        cw.writerow([r['module'], r['n_genes'], '%.6f' % r['median_F_in_vivo'],
                     '%.6f' % r['median_F_in_vitro'], '%.6f' % r['retention_pct']])
w('')
w('wrote %s (%d modules)' % (csvret, len(MODRET)))

# ---------- 5) 污染校正 ----------
w('')
w('=' * 118)
w('造血污染校正：剔除与 PTPRC(CD45) 共变的基因后，部位效应还剩多少？')
w('=' * 118)

_lines = gzip.open(os.path.join(BASE,'data/raw/GSE291355_counts.tsv.gz'),'rt',
                   encoding='utf-8', errors='replace').read().split('\n')
hdr = _lines[0].split('\t')
samp = [h.strip() for h in hdr]
donor = [re.match(r'(H\d+)_', s).group(1) for s in samp]
loc   = [re.match(r'H\d+_([a-z]+)_', s).group(1) for s in samp]
state = [re.match(r'H\d+_[a-z]+_([a-z]+)_adip', s).group(1) for s in samp]
G = []; X = []
kecp = set(ensg)
for _ln in _lines[1:]:
    if not _ln.strip(): continue
    p = _ln.split('\t')
    g = p[0].split('.')[0]
    if g not in kecp: continue
    G.append(g); X.append([float(v) if v not in ('','NA') else 0.0 for v in p[1:]])
X = np.array(X); L = np.log2(X / X.sum(axis=0) * 1e6 + 1)
w('matrix: %d genes x %d samples' % X.shape)
gpos = {g: i for i, g in enumerate(G)}

pi = [i for i, s in enumerate(state) if s == 'primary']
w('PTPRC(CD45) 体内表达(log2CPM): ' + ', '.join('%.2f' % L[gpos['ENSG00000081237'], i] for i in pi))
cd45 = L[gpos['ENSG00000081237'], pi]
w('PTPRC 体内 range: %.2f - %.2f' % (cd45.min(), cd45.max()))

# 与 CD45 的相关（仅在 primary 12 样本上）
Yp = L[:, pi]
Yc = Yp - Yp.mean(axis=1, keepdims=True)
cc = (Yc * (cd45 - cd45.mean())).sum(axis=1) / np.sqrt((Yc**2).sum(axis=1) * ((cd45-cd45.mean())**2).sum())
n_hi = int((np.abs(cc) > 0.7).sum()); n_hi8 = int((np.abs(cc) > 0.8).sum())
w('|r(CD45)| > 0.7 : %d genes (%.1f%%) ; > 0.8 : %d' % (n_hi, 100*n_hi/len(G), n_hi8))

# partial F 检验： y ~ donor + loc  vs  y ~ donor
def design_primary():
    ds = sorted(set(donor[i] for i in pi)); ls = ['meta','epi','scat']
    Xf = []; Xr = []
    for i in pi:
        b = [1.0] + [1.0 if donor[i]==d else 0.0 for d in ds[1:]]
        Xr.append(b)
        Xf.append(b + [1.0 if loc[i]==l else 0.0 for l in ls[1:]])
    return np.array(Xf), np.array(Xr)

Xf, Xr = design_primary()
def pf(Y):
    df_f = Y.shape[0] - np.linalg.matrix_rank(Xf)
    df_r = Y.shape[0] - np.linalg.matrix_rank(Xr)
    def rss(X):
        B = np.linalg.lstsq(X, Y, rcond=None)[0]
        R = Y - X @ B
        return (R**2).sum(axis=0)
    rf, rr = rss(Xf), rss(Xr)
    F = ((rr - rf) / (df_f - df_r)) / (rf / df_f)
    return F, df_f - df_r, df_f

def bh(p):
    p = np.asarray(p, float); n = len(p)
    o = np.argsort(p); rank = np.empty(n); rank[o] = np.arange(1, n+1)
    q = p * n / rank
    q = np.minimum.accumulate(q[o][::-1])[::-1]
    out = np.empty(n); out[o] = np.clip(q, 0, 1)
    return out

from scipy import stats
# Yp 是 genes x samples；partial F 需要 samples x genes
F, d1, d2 = pf(Yp.T)
p = stats.f.sf(F, d1, d2); q = bh(p)
w('')
w('【全基因】部位效应 FDR<0.05 = %d' % int((q < 0.05).sum()))
mask = np.abs(cc) <= 0.7
q_clean = bh(p[mask])
w('【剔除 |r(CD45)|>0.7 的 %d 个基因后】剩 %d 基因，FDR<0.05 = %d' % (n_hi, mask.sum(), int((q_clean < 0.05).sum())))
mask2 = np.abs(cc) <= 0.5
w('【更严：剔除 |r(CD45)|>0.5 的 %d 个基因后】剩 %d 基因，FDR<0.05 = %d' % (
    int((np.abs(cc) > 0.5).sum()), mask2.sum(), int((bh(p[mask2]) < 0.05).sum())))

# 体外同样做（体外 CD45 全 0 → 相关不可算，直接看体外 FDR）
w('')
w('注：体外 diff 全部样本 CD45 表达 = %.2f（培养清除造血细胞）' % L[gpos['ENSG00000081237'], [i for i,s in enumerate(state) if s=="diff"]].max())

with open(LOG, 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(buf))
print('\n[log]', LOG)
