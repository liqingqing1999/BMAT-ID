# -*- coding: utf-8 -*-
"""决定性敏感性分析：剔除"造血污染指纹基因"后，部位身份还剩多少被抹平？"""
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

log = open(os.path.join(BASE, 'logs/clean_out.txt'), 'w', encoding='utf-8')
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
    genes.append(p[0].strip().split('.')[0])
    X.append([float(p[c]) if p[c] not in ('', 'NA') else 0.0 for c in cols])
X = np.array(X); L = np.log2(X / X.sum(axis=0) * 1e6 + 1)
lab = [re.match(r'(H\d+)_(epi|meta|scat)_(primary|diff)_adip', n.strip(), re.I).groups() for n in names]
donor = np.array([l[0] for l in lab]); site = np.array([l[1].upper() for l in lab]); state = np.array([l[2] for l in lab])
gidx = {g: i for i, g in enumerate(genes)}

prim = np.where(state == 'primary')[0]
CD45 = L[gidx['ENSG00000081237'], prim]

# ---- 计算每个基因与 CD45 的相关（12 个 primary 样本）→ 污染指纹 ----
w('=' * 96)
w('1) 用 CD45(PTPRC) 反推"造血污染指纹"：与 CD45 高度共变的基因')
w('=' * 96)
Lp = L[:, prim]
cd = (CD45 - CD45.mean()) / (CD45.std() + 1e-12)
Xc = (Lp - Lp.mean(axis=1, keepdims=True)) / (Lp.std(axis=1, keepdims=True) + 1e-12)
r = (Xc * cd).mean(axis=1)
# 只在对 CD45 有表达意义的基因上评估（避免低表达基因的伪相关）
expressed = (Lp.max(axis=1) > 3)
w('  可评估基因（primary 中 max log2CPM>3）: %d / %d' % (expressed.sum(), len(genes)))
for thr in [0.5, 0.6, 0.7, 0.8, 0.9]:
    n = int(((r > thr) & expressed).sum())
    w('   与 CD45 相关 r > %.1f 的基因: %5d 个' % (thr, n))

CONTAM = (r > 0.7) & expressed
w('')
w('  污染指纹基因示例（r>0.7，按 r 降序前 25，映射前先看 Ensembl ID）:')
order = np.argsort(r)[::-1]
cnt = 0
for i in order:
    if not CONTAM[i]: continue
    w('     %-20s r=%.3f  primary均值=%.2f' % (genes[i], r[i], Lp[i].mean()))
    cnt += 1
    if cnt >= 25: break

# ---- 剔除污染基因后重算 site 效应 ----
def site_F(mask_state, keep):
    idx = np.where(state == mask_state)[0]
    Ls = L[keep][:, idx]; ss = site[idx]
    gs = [np.where(ss == s)[0] for s in ['EPI', 'META', 'SCAT']]
    gm = Ls.mean(axis=1, keepdims=True)
    b = sum(len(g) * (Ls[:, g].mean(axis=1, keepdims=True) - gm) ** 2 for g in gs) / 2
    wv = sum(((Ls[:, g] - Ls[:, g].mean(axis=1, keepdims=True)) ** 2).sum(axis=1, keepdims=True) for g in gs) / (Ls.shape[1] - 3)
    return (b / (wv + 1e-12)).flatten()

from scipy.stats import f as fdist
def bh(p):
    n = len(p); o = np.argsort(p); q = np.empty(n)
    q[o] = np.minimum.accumulate((p[o] * n / (np.arange(n) + 1))[::-1])[::-1]
    return np.clip(q, 0, 1)

w('')
w('=' * 96)
w('2) 三种基因集下的"部位身份基因"数 —— 关键：剔除污染后还剩多少被抹平')
w('=' * 96)
SETS = {
 '全部基因': np.ones(len(genes), bool),
 '剔除 CD45 指纹(r>0.7)': ~CONTAM,
 '剔除 CD45 指纹(r>0.5)': ~((r > 0.5) & expressed),
}
for nm, keep in SETS.items():
    Fp = site_F('primary', keep); Fd = site_F('diff', keep)
    qp = bh(fdist.sf(Fp, 2, 9)); qd = bh(fdist.sf(Fd, 2, 9))
    w('')
    w('  【%s】保留 %d 个基因' % (nm, keep.sum()))
    for a in [0.05, 0.01]:
        np_ = int((qp < a).sum()); nd = int((qd < a).sum())
        w('     FDR<%-5s 体内 %4d 个 | 体外 %3d 个   (保留 %.1f%%)' % (a, np_, nd, 100.0 * nd / max(np_, 1)))
    w('     F 分布  primary 中位 %.2f/90%%%.2f  |  diff 中位 %.2f/90%%%.2f' % (
        np.median(Fp), np.percentile(Fp, 90), np.median(Fd), np.percentile(Fd, 90)))
    # 存下非污染 top 基因
    if nm.startswith('剔除 CD45 指纹(r>0.7)'):
        ids = np.where(keep)[0]
        topl = ids[np.argsort(Fp)[::-1][:120]]
        w('     非污染 top120 部位身份基因（Ensembl）:')
        w('       ' + ', '.join(genes[i] for i in topl))
        np.save(os.path.join(BASE, 'data/processed/noncontam_top.npy'), topl)

# ---- 成骨/成脂关键基因的部位差异（污染不敏感的基因）----
w('')
w('=' * 96)
w('3) 关键功能基因的"部位落差"在 体内 vs 体外（真实生物学，非污染）')
w('=' * 96)
KEY = {'RUNX2': 'ENSG00000124813', 'SP7': 'ENSG00000170374', 'COL1A1': 'ENSG00000108821',
       'ALPL': 'ENSG00000162551', 'SPP1': 'ENSG00000118785', 'IBSP': 'ENSG00000029559',
       'PPARG': 'ENSG00000132170', 'FABP4': 'ENSG00000170323', 'ADIPOQ': 'ENSG00000181092',
       'PLIN1': 'ENSG00000166819', 'CEBPA': 'ENSG00000245848', 'CXCL12': 'ENSG00000107562',
       'VCAM1': 'ENSG00000162692', 'LEPR': 'ENSG00000116678', 'SFRP1': 'ENSG00000104332'}
w('%-9s %26s %10s | %26s %10s' % ('gene', '体内 primary  E/M/S', '落差', '体外 diff  E/M/S', '落差'))
for g, eid in KEY.items():
    i = gidx.get(eid)
    if i is None: continue
    mp = [L[i, np.where((state == 'primary') & (site == s))[0]].mean() for s in ['EPI', 'META', 'SCAT']]
    md = [L[i, np.where((state == 'diff') & (site == s))[0]].mean() for s in ['EPI', 'META', 'SCAT']]
    rp = max(mp) - min(mp); rd = max(md) - min(md)
    w('%-9s  %7.2f %7.2f %7.2f %9.2f | %7.2f %7.2f %7.2f %9.2f   %s' % (
        g, mp[0], mp[1], mp[2], rp, md[0], md[1], md[2], rd,
        '体外落差↓' if rd < rp * 0.5 else ('体外落差保持' if rd > rp * 0.8 else '部分↓')))

log.close()
print('DONE')
