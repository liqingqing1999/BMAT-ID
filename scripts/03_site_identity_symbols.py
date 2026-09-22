# -*- coding: utf-8 -*-
"""把 top 部位身份基因映射为 symbol，并按部位分类（EPI/META/SCAT-high）"""
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

import os, re, gzip, json, numpy as np, urllib.request, ssl

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
log = open(os.path.join(BASE, 'logs/symbols_out.txt'), 'w', encoding='utf-8')
def w(s=''):
    print(s); log.write(str(s) + '\n')

rows = gzip.open(os.path.join(BASE, 'data/raw/GSE291355_counts.tsv.gz'), 'rt', encoding='utf-8', errors='replace').read().split('\n')
hdr = rows[0].split('\t')
samp = [(i + 1, h.strip()) for i, h in enumerate(hdr) if re.match(r'H\d+_', h.strip())]
names = [n for _, n in samp]; cols = [i for i, _ in samp]
genes, X = [], []
for ln in rows[1:]:
    if not ln.strip(): continue
    p = ln.split('\t')
    if len(p) <= max(cols): continue
    genes.append(p[0].strip())
    X.append([float(p[c]) if p[c] not in ('', 'NA') else 0.0 for c in cols])
X = np.array(X); L = np.log2(X / X.sum(axis=0) * 1e6 + 1)
lab = [re.match(r'(H\d+)_(epi|meta|scat)_(primary|diff)_adip', n.strip(), re.I).groups() for n in names]
donor = np.array([l[0] for l in lab]); site = np.array([l[1].upper() for l in lab]); state = np.array([l[2] for l in lab])
SITES = ['EPI', 'META', 'SCAT']

# F 值
def site_F(mask_state):
    idx = np.where(state == mask_state)[0]; Ls = L[:, idx]; ss = site[idx]
    gs = [np.where(ss == s)[0] for s in SITES]
    gm = Ls.mean(axis=1, keepdims=True)
    b = sum(len(g) * (Ls[:, g].mean(axis=1, keepdims=True) - gm) ** 2 for g in gs) / (len(gs) - 1)
    wv = sum(((Ls[:, g] - Ls[:, g].mean(axis=1, keepdims=True)) ** 2).sum(axis=1, keepdims=True) for g in gs) / (Ls.shape[1] - len(gs))
    return (b / (wv + 1e-12)).flatten()

Fp = site_F('primary'); Fd = site_F('diff')

# 取 top 150 部位身份基因（体内）
top = np.argsort(Fp)[::-1][:150]
top_ids = [genes[i].split('.')[0] for i in top]

# MyGene 批量映射
sym = {}
try:
    data = urllib.parse.urlencode({'q': ','.join(top_ids), 'scopes': 'ensembl.gene', 'fields': 'symbol,name',
                                   'species': 'human'}).encode()
    req = urllib.request.Request('https://mygene.info/v3/query', data=data,
                                 headers={'Content-Type': 'application/x-www-form-urlencoded'})
    res = json.loads(urllib.request.urlopen(req, context=ctx, timeout=60).read())
    hits = res if isinstance(res, list) else res.get('hits', [])
    for h in hits:
        q = h.get('query')
        if q and not h.get('notfound'):
            sym[q] = (h.get('symbol') or h.get('name') or '?')
except Exception as e:
    w('MyGene ERROR: %s' % e)

def site_means(gi):
    out = {}
    for s in SITES:
        idx = np.where((state == 'primary') & (site == s))[0]
        out[s] = float(L[gi, idx].mean())
    return out

w('')
w('=' * 88)
w('体内最强的 60 个"部位身份基因"：符号 / 体内三部位均值 / 体外 F 值')
w('=' * 88)
w('%-4s %-14s %8s %8s %8s %6s %9s %8s' % ('#', 'symbol', 'EPI', 'META', 'SCAT', 'range', 'F_prim', 'F_diff'))
cnt = {}
for k, gi in enumerate(top[:60]):
    eid = genes[gi].split('.')[0]
    s = sym.get(eid, eid)
    m = site_means(gi)
    rng = max(m.values()) - min(m.values())
    hi = max(m, key=m.get)
    cnt[hi] = cnt.get(hi, 0) + 1
    w('%-4d %-14s %8.2f %8.2f %8.2f %6.2f %9.1f %8.2f   [%s-high]' % (
        k + 1, s, m['EPI'], m['META'], m['SCAT'], rng, Fp[gi], Fd[gi], hi))

w('')
w('top150 中最高的部位分布: %s' % cnt)
w('')
w('--- 全部 top150 的 symbol 列表（按 F 降序）---')
w(', '.join('%s(%.0f)' % (sym.get(genes[i].split('.')[0], genes[i].split('.')[0]), Fp[i]) for i in top[:150]))

log.close()
print('DONE')
