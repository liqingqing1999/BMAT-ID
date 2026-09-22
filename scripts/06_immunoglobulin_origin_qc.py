# -*- coding: utf-8 -*-
"""质控：部位身份基因里的免疫球蛋白基因，是 BMAd 自身表达还是造血细胞污染？"""
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

log = open(os.path.join(BASE, 'logs/qc_out.txt'), 'w', encoding='utf-8')
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

MARK = {
 'PTPRC(CD45)': 'ENSG00000081237', 'CD3E(T)': 'ENSG00000198851', 'CD19(B)': 'ENSG00000177455',
 'MS4A1(CD20)': 'ENSG00000156738', 'MZB1(浆)': 'ENSG00000170476', 'JCHAIN(浆)': 'ENSG00000132465',
 'SDC1(CD138)': 'ENSG00000115884', 'IGHG2': 'ENSG00000211893', 'IGKC': 'ENSG00000211592',
 'PECAM1(CD31)': 'ENSG00000110799', 'VWF(内皮)': 'ENSG00000110799',
 'ADIPOQ(脂肪)': 'ENSG00000181092', 'PLIN1(脂肪)': 'ENSG00000166819', 'FABP4(脂肪)': 'ENSG00000170323',
 'PPARG(脂肪)': 'ENSG00000132170', 'LEP(脂肪)': 'ENSG00000174697', 'PLIN4': 'ENSG00000167699',
 'RUNX2(骨)': 'ENSG00000124813', 'SP7(骨)': 'ENSG00000170374', 'COL1A1(骨)': 'ENSG00000108821',
 'CXCL12(龛)': 'ENSG00000107562', 'VCAM1(龛)': 'ENSG00000162692', 'KITL(龛)': 'ENSG00000068758',
 'PTPRC-v2': 'ENSG00000081237',
}

w('=' * 100)
w('标志基因 log2CPM —— 按 状态×部位 均值（primary=体内原代, diff=体外诱导）')
w('=' * 100)
hdr2 = '%-16s' % 'gene'
for st in ['primary', 'diff']:
    for s in ['EPI', 'META', 'SCAT']:
        hdr2 += ' %9s' % ('%s_%s' % (st[:4], s))
w(hdr2)
for name, eid in MARK.items():
    if name.endswith('-v2'): continue
    i = gidx.get(eid)
    if i is None:
        w('%-16s   (not in matrix)' % name); continue
    line = '%-16s' % name
    for st in ['primary', 'diff']:
        for s in ['EPI', 'META', 'SCAT']:
            idx = np.where((state == st) & (site == s))[0]
            line += ' %9.2f' % L[i, idx].mean()
    w(line)

# 每个样本的 PTPRC / MZB1 用于看离散
w('')
w('=' * 100)
w('逐样本 PTPRC(CD45) 与 MZB1 —— 若某部位系统性偏高即污染信号')
w('=' * 100)
for st in ['primary', 'diff']:
    for s in ['EPI', 'META', 'SCAT']:
        idx = np.where((state == st) & (site == s))[0]
        vals45 = L[gidx['ENSG00000081237'], idx]
        valsmz = L[gidx.get('ENSG00000170476', 0), idx]
        w('  %-8s %-5s  CD45=[%s]  MZB1=[%s]' % (st, s, ' '.join('%5.2f' % v for v in vals45),
                                                 ' '.join('%5.2f' % v for v in valsmz)))

# 免疫基因（IGH/IGK/IGL 前缀）在 top150 中的占比，按部位
w('')
w('=' * 100)
w('免疫球蛋白基因在"部位身份"中的贡献（体内 primary 的 site F 值 top150）')
w('=' * 100)
allsym = {}
import re as _re
def site_F(mask_state):
    idx = np.where(state == mask_state)[0]; Ls = L[:, idx]; ss = site[idx]
    gs = [np.where(ss == s)[0] for s in ['EPI', 'META', 'SCAT']]
    gm = Ls.mean(axis=1, keepdims=True)
    b = sum(len(g) * (Ls[:, g].mean(axis=1, keepdims=True) - gm) ** 2 for g in gs) / 2
    wv = sum(((Ls[:, g] - Ls[:, g].mean(axis=1, keepdims=True)) ** 2).sum(axis=1, keepdims=True) for g in gs) / (Ls.shape[1] - 3)
    return (b / (wv + 1e-12)).flatten()
Fp = site_F('primary')
top = np.argsort(Fp)[::-1][:150]
w('  （基因名可能是 Ensembl ID，用 ID 前缀判断）')
w('  top150 全部为 Ensembl ID，免疫球蛋白基因需另行映射；此处给出 F 值排名统计：')
w('   top150 中 site 效应最强：#1 F=%.0f, #50 F=%.1f, #150 F=%.1f' % (Fp[top[0]], Fp[top[49]], Fp[top[149]]))
w('   primary 中 FDR<0.05 的部位差异基因数 = 112（见 quant_out.txt）')

log.close()
print('DONE')
