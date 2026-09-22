## 关键论据：体内 vs 体外，部位标签与组成轴的共线性对比
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

import os, gzip, re, csv
import numpy as np

OUT  = os.path.join(BASE, 'logs/varpart_collinearity_out.txt')
con  = open(OUT, 'w', encoding='utf-8')
def w(*a):
    s = ' '.join(str(x) for x in a); print(s); con.write(s + '\n'); con.flush()

# 表达矩阵
rows = gzip.open(os.path.join(BASE, 'data/raw/GSE291355_counts.tsv.gz'), 'rt',
                 encoding='utf-8', errors='replace').read().split('\n')
hdr = rows[0].split('\t')
samp = [h.strip() for h in hdr if re.match(r'H\d+_', h.strip())]
genes, X = [], []
for ln in rows[1:]:
    if not ln.strip():
        continue
    p = ln.split('\t')
    if len(p) < len(samp) + 1:
        continue
    genes.append(p[0])
    X.append([float(v) if v not in ('', 'NA') else 0.0 for v in p[1:1 + len(samp)]])
X = np.array(X, float)
lab = [re.match(r'(H\d+)_(epi|meta|scat)_(primary|diff)_adip', s, re.I).groups() for s in samp]
site = np.array([l[1].lower() for l in lab]); state = np.array([l[2] for l in lab])
L = np.log2(X / X.sum(0) * 1e6 + 1)
gidx = {g.split('.')[0]: i for i, g in enumerate(genes)}   # marker 表不带版本号
pi = np.where(state == 'primary')[0]; di = np.where(state == 'diff')[0]

mk = [l.rstrip('\n').split('\t') for l in open(os.path.join(BASE, 'results/tables/BMATID_marker_ensg.tsv'),
                                               encoding='utf-8') if l.strip()][1:]
AX = {}
for ax, sym, ens in mk:
    AX.setdefault(ax, []).append(ens)

def axis_score(ens_list, idx):
    rr = [gidx[e] for e in ens_list if e in gidx]
    Z = L[np.ix_(rr, idx)]
    keep = Z.std(1) > 0
    Z = Z[keep]
    if Z.shape[0] == 0:
        return np.full(len(idx), np.nan), 0
    Z = (Z - Z.mean(1, keepdims=True)) / (Z.std(1, keepdims=True) + 1e-9)
    return Z.mean(0), int(keep.sum())

w('=' * 96)
w('  部位标签与组成轴的共线性：体内 vs 体外')
w('=' * 96)
w('')
w('  统计量：组成轴分数 ~ 部位标签 的 R²（组间方差占比）+ 轴分数在三个部位的均值')
w('  预期：体内 R² 很高（部位标签等价于组成梯度）；体外 R² 应接近 0（无造血/骨细胞）')
w('')

for tag, idx in [('体内 primary (n=12)', pi), ('体外 diff (n=12)', di)]:
    w('### %s' % tag)
    SS_tot_all = []
    for ax in ['Hemato', 'Bone', 'Adipo']:
        s, used = axis_score(AX[ax], idx)
        g = site[idx]
        # 组间 R²（eta²）
        ss_tot = ((s - s.mean()) ** 2).sum()
        ss_bet = sum(len(s[g == k]) * (s[g == k].mean() - s.mean()) ** 2 for k in ['epi', 'meta', 'scat'] if (g == k).any())
        r2 = ss_bet / ss_tot if ss_tot > 0 else float('nan')
        mus = {k: (s[g == k].mean() if (g == k).any() else float('nan')) for k in ['epi', 'meta', 'scat']}
        w('  %-8s 标记 %2d 个  R²(部位) = %.3f   EPI %+.2f  META %+.2f  SCAT %+.2f   极差 %.2f' %
          (ax, used, r2, mus['epi'], mus['meta'], mus['scat'],
           np.nanmax(list(mus.values())) - np.nanmin(list(mus.values()))))
    w('')

# 直接算：部位 one-hot 与组成轴的多元 R²
w('  【多元】用两个组成轴同时预测部位（多分类，Fisher 判别方向上的解释度）')
for tag, idx in [('体内', pi), ('体外', di)]:
    sH, _ = axis_score(AX['Hemato'], idx); sB, _ = axis_score(AX['Bone'], idx)
    g = site[idx]
    A = np.column_stack([np.ones(len(idx)), sH, sB])
    y = (g == 'meta').astype(float) - (g == 'scat').astype(float)   # meta vs scat 对比
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ beta
    r2 = 1 - ((y - pred) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    w('    %s：预测 META-vs-SCAT 方向  R² = %.3f' % (tag, r2))

w('')
w('  结论：若体内 R² 高、体外 R² ≈ 0 → "部位=组成"的共线只在体内成立，')
w('        因此体外部位效应的塌陷（17.1%→3.6%）可作为"体内部位效应主要是组成"的直接证据。')
con.close()
print('->', OUT)
