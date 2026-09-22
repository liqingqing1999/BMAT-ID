# -*- coding: utf-8 -*-
"""Core quantification: how much of the *site-identity* signal does adipogenic
   culture remove, and which site does the residual pattern drift towards?
   Input: the public GSE291355 counts (24 libraries = 4 donors x 3 sites x 2 states)."""
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

import os, re, gzip, numpy as np

log = open(os.path.join(BASE, 'logs/quant_out.txt'), 'w', encoding='utf-8')
def w(s=''):
    print(s); log.write(str(s) + '\n')

# ---- 读计数表 ----
rows = gzip.open(os.path.join(BASE, 'data/raw/GSE291355_counts.tsv.gz'), 'rt', encoding='utf-8', errors='replace').read().split('\n')
hdr = rows[0].split('\t')
print('header n=%d  first=%s  second=%s' % (len(hdr), hdr[0][:30], hdr[1][:30]))
# 24 个样本列名（数据行首列是基因 ID → 样本列偏移 +1）
samp = [(i + 1, h.strip()) for i, h in enumerate(hdr) if re.match(r'H\d+_', h.strip())]
names = [n for _, n in samp]
cols = [i for i, _ in samp]
gi = 0  # 首列是基因 ID
genes, X = [], []
for ln in rows[1:]:
    if not ln.strip(): continue
    p = ln.split('\t')
    if len(p) <= max(cols): continue
    genes.append(p[gi].strip())
    X.append([float(p[c]) if p[c] not in ('', 'NA') else 0.0 for c in cols])
X = np.array(X, dtype=float)
w('matrix: %d genes x %d samples' % X.shape)

# ---- log2CPM ----
lib = X.sum(axis=0)
cpm = X / lib * 1e6
L = np.log2(cpm + 1)

# ---- 解析标签 ----
lab = []
for n in names:
    m = re.match(r'(H\d+)_(epi|meta|scat)_(primary|diff)_adip', n.strip(), re.I)
    lab.append((m.group(1), m.group(2).upper(), m.group(3)))
donor = np.array([l[0] for l in lab])
site  = np.array([l[1] for l in lab])
state = np.array([l[2] for l in lab])
w('samples: %s' % names)

SITES = ['EPI', 'META', 'SCAT']

# ---- 1. diff 样本与各部位 primary 的平均相关 ----
C = np.corrcoef(L.T)
def sub(mask):
    return np.where(mask)[0]
w('')
w('=' * 70)
w('1) 体外诱导样本(diff) 与 各部位体内原代(primary) 的平均相关')
w('=' * 70)
for s in SITES:
    a = sub(state == 'diff')
    b = sub((state == 'primary') & (site == s))
    v = [C[i, j] for i in a for j in b]
    w('  diff  vs  primary_%-5s : r = %.3f  (n=%d pairs)' % (s, np.mean(v), len(v)))
w('  （参考）primary 内部各部位之间：')
for i, s1 in enumerate(SITES):
    for s2 in SITES[i+1:]:
        a = sub((state == 'primary') & (site == s1)); b = sub((state == 'primary') & (site == s2))
        v = [C[i, j] for i in a for j in b]
        w('     primary_%s vs primary_%s : r = %.3f' % (s1, s2, np.mean(v)))
w('  （参考）diff 内部各部位之间：')
for i, s1 in enumerate(SITES):
    for s2 in SITES[i+1:]:
        a = sub((state == 'diff') & (site == s1)); b = sub((state == 'diff') & (site == s2))
        v = [C[i, j] for i in a for j in b]
        w('     diff_%s vs diff_%s : r = %.3f' % (s1, s2, np.mean(v)))

# ---- 2. 部位效应：每个基因的 F 统计量（状态内）----
def site_effect(mask_state):
    idx = sub(state == mask_state)
    Ls = L[:, idx]
    ss = site[idx]
    groups = [np.where(ss == s)[0] for s in SITES]
    n_tot = Ls.shape[1]
    gm = Ls.mean(axis=1, keepdims=True)
    ss_between = sum(len(g) * (Ls[:, g].mean(axis=1, keepdims=True) - gm) ** 2 for g in groups)
    df_b = len(groups) - 1
    ss_within = sum(((Ls[:, g] - Ls[:, g].mean(axis=1, keepdims=True)) ** 2).sum(axis=1, keepdims=True) for g in groups)
    df_w = n_tot - len(groups)
    ms_b = ss_between / df_b
    ms_w = ss_within / df_w
    F = (ms_b / (ms_w + 1e-12)).flatten()
    return F

F_prim = site_effect('primary')
F_diff = site_effect('diff')

def bh_fdr_p(F, df1, df2):
    from math import lgamma, log, exp
    # F 检验 p 值（用 scipy 若可用）
    try:
        from scipy.stats import f as fdist
        return fdist.sf(F, df1, df2)
    except Exception:
        return None

p_prim = bh_fdr_p(F_prim, 2, 9)
p_diff = bh_fdr_p(F_diff, 2, 9)

def bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

w('')
w('=' * 70)
w('2) "部位身份基因"数量（状态内 site 效应，线性模型 F 检验）')
w('=' * 70)
if p_prim is not None:
    q_prim = bh(p_prim); q_diff = bh(p_diff)
    for alpha in [0.05, 0.01, 0.001]:
        np_ = int((q_prim < alpha).sum()); nd = int((q_diff < alpha).sum())
        w('  FDR<%-6s : 体内 primary 有 %5d 个部位差异基因 | 体外 diff 只剩 %4d 个  (保留 %.1f%%)'
          % (alpha, np_, nd, 100.0 * nd / max(np_, 1)))
    w('')
    w('  F 统计量分布：')
    w('     primary : 中位 %.2f  90%%分位 %.2f  最大 %.1f' % (np.median(F_prim), np.percentile(F_prim, 90), F_prim.max()))
    w('     diff    : 中位 %.2f  90%%分位 %.2f  最大 %.1f' % (np.median(F_diff), np.percentile(F_diff, 90), F_diff.max()))
    # 前 20 个部位身份基因在 diff 中的 F
    top = np.argsort(F_prim)[::-1][:30]
    w('')
    w('  体内部位身份最强的 30 个基因，其 F 值在体内 vs 体外：')
    w('     %-22s %8s %8s %8s' % ('gene', 'F_prim', 'F_diff', 'ratio'))
    for i in top:
        w('     %-22s %8.2f %8.2f %8.3f' % (genes[i], F_prim[i], F_diff[i], F_diff[i] / F_prim[i]))
else:
    w('  (scipy 不可用，跳过 p 值)')

log.close()
print('DONE')
