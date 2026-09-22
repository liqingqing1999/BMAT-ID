## 诊断 M3 双峰：共线导致的方差分配不稳定？
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

import os, csv, math
T = os.path.join(BASE, "results", "tables")
OUT = os.path.join(BASE, "logs", "varpart_diag_out.txt")
con = open(OUT, 'w', encoding='utf-8')
def w(*a):
    s = ' '.join(str(x) for x in a); print(s); con.write(s + '\n')

def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

vivo  = {r['ensg']: r for r in load('BMATID_varpart_invivo.csv')}
vitro = {r['ensg']: r for r in load('BMATID_varpart_invitro.csv')}
comp  = {r['ensg']: r for r in load('BMATID_varpart_invivo_comp.csv')}

F = lambda d, g, k: float(d[g][k]) if g in d and d[g].get(k) not in (None, '', 'NA') else float('nan')

w('=' * 104)
w('  M3 双峰诊断：部位项的高方差是否与组成项"共享"（共线症状）')
w('=' * 104)
w('')
w('  判据：若 %loc(M3) 高的基因同时 %Hemato+%Bone 也高 → 两者在抢同一批方差（不可分离）')
w('        若 %loc(M3) 高的基因 %Hemato+%Bone 低 → 该部位方差独立于组成（可解读为真信号）')
w('')

bins = [(0, 1e-9, '= 0'), (1e-9, .1, '0–10%'), (.1, .2, '10–20%'), (.2, .5, '20–50%'),
        (.5, .8, '50–80%'), (.8, 1.01, '80–100%')]
w('  %-10s %6s %10s %12s %12s %12s' % ('%loc(M3)', 'n', '%loc M1中位', '%Hemato中位', '%Bone中位', '组成合计中位'))
def med(v):
    v = sorted(x for x in v if not math.isnan(x))
    return v[len(v) // 2] if v else float('nan')

for lo, hi, tag in bins:
    gs = [g for g in comp if lo <= F(comp, g, 'localization') < hi]
    if not gs:
        continue
    m1  = med([F(vivo, g, 'localization') for g in gs])
    mh  = med([F(comp, g, 'Hemato') for g in gs])
    mb  = med([F(comp, g, 'Bone') for g in gs])
    mt  = med([F(comp, g, 'Hemato') + F(comp, g, 'Bone') for g in gs])
    w('  %-10s %6d %10.3f %12.3f %12.3f %12.3f' % (tag, len(gs), m1, mh, mb, mt))

w('')
w('  【解读】若 %loc(M3) 越高的箱，组成合计中位越低 → 两者互补，部位项的残余是"独立"的')
w('')

# 关键：loc(M3)>0.5 的基因，在 M1 里 %loc 是多少？在 M2（体外）是多少？
hi3 = [g for g in comp if F(comp, g, 'localization') > .5]
w('  部位方差(M3) > 50%% 的基因：%d 个' % len(hi3))
w('    其中 部位方差(M1) 也 > 0.5 的：%d (%.1f%%)' % (
    sum(1 for g in hi3 if F(vivo, g, 'localization') > .5),
    100 * sum(1 for g in hi3 if F(vivo, g, 'localization') > .5) / len(hi3)))
w('    其中 体外部位方差(M2) > 0.2 的：%d (%.1f%%)' % (
    sum(1 for g in hi3 if F(vitro, g, 'localization') > .2),
    100 * sum(1 for g in hi3 if F(vitro, g, 'localization') > .2) / len(hi3)))
w('')
w('  部位方差(M3) = 0 的基因：%d 个（占 %.1f%%）' % (
    sum(1 for g in comp if F(comp, g, 'localization') <= 1e-9),
    100 * sum(1 for g in comp if F(comp, g, 'localization') <= 1e-9) / len(comp)))
w('    → 这些基因的部位方差被组成协变量完全吸收')

# M1 -> M3 的分层保留
w('')
w('  【M1 → M3 分层保留】')
for lo, hi, tag in [(.5, 1.01, 'M1 %loc>50%'), (.8, 1.01, 'M1 %loc>80%'), (.2, .5, 'M1 %loc 20-50%'), (0, .2, 'M1 %loc<20%')]:
    a = [g for g in vivo if lo <= F(vivo, g, 'localization') < hi]
    b = [g for g in a if g in comp and F(comp, g, 'localization') > lo]
    if not a:
        continue
    mt = med([F(comp, g, 'Hemato') + F(comp, g, 'Bone') for g in a])
    w('  %-16s n=%5d  →  M3 仍 >%s : %5d (%.1f%%)   该组组成合计中位 %.3f' %
      (tag, len(a), ('50%' if lo == .5 else '80%' if lo == .8 else '20%'), len(b), 100 * len(b) / len(a), mt))

# 残差项对比
w('')
w('  【残差方差（Residuals）中位】M1 = %.3f   M2 = %.3f   M3 = %.3f' % (
    1 - med([F(vivo, g, 'localization') + F(vivo, g, 'donor') for g in vivo]),
    1 - med([F(vitro, g, 'localization') + F(vitro, g, 'donor') for g in vitro]),
    1 - med([F(comp, g, 'localization') + F(comp, g, 'donor') + F(comp, g, 'Hemato') + F(comp, g, 'Bone') for g in comp])))

con.close()
print('->', OUT)
