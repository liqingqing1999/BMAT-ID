# -*- coding: utf-8 -*-
"""Figure 1 (four panels)
   a  variance composition across four mixed models (incl. reference-anchored)
   b  ECDF of per-gene site-attributable variance
   c  POSITIVE CONTROL: axis-patterning TFs keep site variance, composition markers lose it
   d  survival of strong site genes by tier
"""
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

import os, csv, json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG  = os.path.join(BASE, 'results/figures')
T    = os.path.join(BASE, 'results/tables')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 400, 'figure.facecolor': 'white', 'axes.facecolor': 'white'})
C_VIVO, C_VITRO, C_COMP = '#2F4B7C', '#F28E2B', '#8C8C8C'

def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))
def gv(rows, key='localization'):
    out = {}
    for r in rows:
        try: out[r['ensg']] = float(r.get(key, ''))
        except (TypeError, ValueError): pass
    return out

def gvd(d, key='localization'):
    out = {}
    for k, r in d.items():
        try: out[k] = float(r.get(key, ''))
        except (TypeError, ValueError): pass
    return out

d1 = {r['ensg']: r for r in load('BMATID_varpart_invivo.csv')}
d2 = {r['ensg']: r for r in load('BMATID_varpart_invitro.csv')}
d3 = {r['ensg']: r for r in load('BMATID_varpart_invivo_comp.csv')}
d4 = {r['ensg']: r for r in load('BMATID_varpart_M7_reference.csv')}   # reference-anchored model
M1, M2, M3, M4 = gvd(d1), gvd(d2), gvd(d3), gvd(d4)

def mean_of(d, key):
    v = [float(x[key]) for x in d.values() if x.get(key) not in (None, '', 'NA')]
    return 100 * float(np.mean(v))

comp_M1 = dict(loc=mean_of(d1, 'localization'), donor=mean_of(d1, 'donor'), hem=0.0, bone=0.0)
comp_M2 = dict(loc=mean_of(d2, 'localization'), donor=mean_of(d2, 'donor'), hem=0.0, bone=0.0)
comp_M3 = dict(loc=mean_of(d3, 'localization'), donor=mean_of(d3, 'donor'),
               hem=mean_of(d3, 'Hemato'), bone=mean_of(d3, 'Bone'))
# reference-anchored model: hand-built axes replaced by independent-atlas signatures
comp_M7 = dict(loc=mean_of(d4, 'localization'), donor=mean_of(d4, 'donor'),
               hem=mean_of(d4, 'Immune_ref'), bone=mean_of(d4, 'Skeletal_ref'))
for c in (comp_M1, comp_M2, comp_M3, comp_M7):
    c['resid'] = 100 - (c['loc'] + c['donor'] + c['hem'] + c['bone'])
print('M1', {k: round(v, 1) for k, v in comp_M1.items()})
print('M2', {k: round(v, 1) for k, v in comp_M2.items()})
print('M3', {k: round(v, 1) for k, v in comp_M3.items()})
print('M7', {k: round(v, 1) for k, v in comp_M7.items()})
l1, l2, l3 = np.array(list(M1.values())), np.array(list(M2.values())), np.array(list(M3.values()))
l4 = np.array(list(gvd(d4).values()))     # reference-anchored model

# 2 x 2 but not equal-width: panel c (33 gene labels + 4 group columns) needs
# the widest box, so the four panels are placed in figure coordinates directly.
fig = plt.figure(figsize=(7.0, 7.2))
BOX = {
    'a': [0.078, 0.642, 0.315, 0.288],
    'b': [0.527, 0.642, 0.428, 0.288],
    'c': [0.078, 0.078, 0.500, 0.462],
    'd': [0.640, 0.078, 0.315, 0.462],
}

# ---------------- a ----------------
ax = fig.add_axes(BOX['a'])
segs = [('loc', 'Site (localization)', '#2F4B7C'), ('donor', 'Donor', '#9DB4D6'),
        ('hem', 'Composition (haemato.)', '#C44E52'), ('bone', 'Composition (bone)', '#8172B2')]
names = ['In vivo', 'In vitro', 'In vivo\n+ axes', 'In vivo\n+ reference']
data = [comp_M1, comp_M2, comp_M3, comp_M7]
N = len(data)
bottom = np.zeros(N); handles = []
for key, lab, col in segs:
    vals = np.array([d[key] for d in data])
    if vals.sum() <= 0: continue
    b = ax.bar(range(N), vals, 0.6, bottom=bottom, color=col, label=lab,
               edgecolor='white', linewidth=0.5, zorder=3)
    handles.append(b); bottom = bottom + vals
rb = ax.bar(range(N), [d['resid'] for d in data], 0.6, bottom=bottom, color='#E8E8E8',
            label='Residual', edgecolor='white', linewidth=0.5, zorder=3)
for i in range(N):
    ax.text(i, 101.5, '%.1f' % data[i]['loc'], ha='center', va='bottom',
            fontsize=7.2, fontweight='bold', color='#2F4B7C')
ax.set_xticks(range(N)); ax.set_xticklabels(names, fontsize=6.6)
ax.set_ylabel('Variance component\n(mean % of gene variance)', fontsize=7.6)
ax.set_ylim(0, 112); ax.set_xlim(-0.6, N - 0.4); ax.set_yticks([0, 20, 40, 60, 80, 100])
fig.legend(handles=handles + [rb], fontsize=6.2, frameon=False, loc='upper center',
           bbox_to_anchor=(0.5, 0.988), ncol=5, columnspacing=0.9, handlelength=1.0,
           handletextpad=0.4)
ax.set_title('a  Variance composition', fontsize=8.6, loc='left', fontweight='bold')

# ---------------- b ----------------
ax = fig.add_axes(BOX['b'])
def ecdf(v):
    x = np.sort(v); y = np.arange(1, len(x) + 1) / len(x)
    return np.concatenate([[0], x]), np.concatenate([[0], y])
for v, c, lab in [(l1, C_VIVO, 'In vivo'), (l2, C_VITRO, 'In vitro'),
                  (l3, C_COMP, 'In vivo + composition'),
                  (l4, '#C44E52', 'In vivo + reference')]:
    xs, ys = ecdf(100 * v); ax.plot(xs, ys, color=c, lw=1.5, label=lab)
ax.set_xlabel('Site-attributable variance per gene (%)', fontsize=7.6)
ax.set_ylabel('Cumulative fraction of genes', fontsize=7.6)
ax.set_xlim(0, 100); ax.set_ylim(0, 1.02)
# exact-zero fractions, computed from the fitted models (tolerance: component == 0)
# The ECDF is discontinuous at x = 0: values between 0 and ~1e-9 are visually packed
# onto the axis, so the *exact-zero* line sits below the curve's apparent intercept.
# Both tiers are therefore labelled, with exact zero as the claim and 1e-9 as transparency.
A3, A4 = np.asarray(list(M3.values())), np.asarray(list(M4.values()))
FZ3, FZ4 = float(np.mean(A3 == 0)), float(np.mean(A4 == 0))
FZ3b, FZ4b = float(np.mean(A3 <= 1e-9)), float(np.mean(A4 <= 1e-9))
print('exact zero: M3=%.1f%% ref=%.1f%% | below 1e-9: M3=%.1f%% ref=%.1f%%'
      % (100 * FZ3, 100 * FZ4, 100 * FZ3b, 100 * FZ4b))
ax.axhline(FZ3, color='#BBB', lw=0.7, ls=':')
ax.axhline(FZ4, color='#E4B9B9', lw=0.7, ls=':')
ax.text(97, 0.72, '%.1f%% of genes: site term = 0\n(%.1f%% with reference)\n'
        'below 1e-9: %.1f%% / %.1f%%' % (100 * FZ3, 100 * FZ4, 100 * FZ3b, 100 * FZ4b),
        fontsize=6.0, ha='right', va='top', color='#777')
ax.legend(fontsize=6.3, frameon=False, loc='lower right')
ax.set_title('b  Site variance per gene', fontsize=8.6, loc='left', fontweight='bold')

# ---------------- d ----------------
ax = fig.add_axes(BOX['d'])
groups = [('>80%', .8, 1.01), ('50\u201380%', .5, .8), ('20\u201350%', .2, .5)]
xpos = np.arange(3); w = 0.34; n_m1, n_m3 = [], []
for tag, lo, hi in groups:
    g_ = [g for g in M1 if lo <= M1[g] < hi]
    n_m1.append(len(g_)); n_m3.append(sum(1 for g in g_ if g in M3 and M3[g] > lo))
ax.bar(xpos - w / 2, n_m1, w, color='#B8C4DA', label='In vivo alone', zorder=3)
ax.bar(xpos + w / 2, n_m3, w, color=C_COMP, label='After composition', zorder=3)
for i in range(3):
    ax.text(i - w / 2, n_m1[i] * 1.35, '%d' % n_m1[i], ha='center', fontsize=6.4, color='#4A5B78')
    ax.text(i + w / 2, n_m3[i] * 1.35, '%d' % n_m3[i], ha='center', fontsize=6.4, color='#555')
    ax.text(i, max(n_m1[i], n_m3[i]) * 0.30, '%.0f%%' % (100 * n_m3[i] / max(1, n_m1[i])),
            ha='center', fontsize=7.4, fontweight='bold', color='#B23A3A')
ax.set_yscale('log'); ax.set_ylim(8, 60000)
ax.set_xticks(xpos); ax.set_xticklabels([g[0] for g in groups], fontsize=7.4)
ax.set_xlabel('Site variance in vivo (%)', fontsize=7.6)
ax.set_ylabel('Number of genes (log)', fontsize=7.6)
ax.legend(fontsize=6.2, frameon=False, loc='upper right', borderaxespad=0.15)
ax.set_title('d  Only 11% of the strongest\nsite genes survive', fontsize=8.6, loc='left', fontweight='bold')

# ---------------- c (positive control, bottom-left) ----------------
axc = fig.add_axes(BOX['c'])
groups_c = {
 "Axis-patterning TF": dict(color="#C1272D", marker="o", genes=[
    ("TBX5", .99, .97), ("RSPO1", .97, .92), ("DLX5", .95, .51), ("BARX1", .95, .59),
    ("TBX18", .81, .70), ("PYGO1", .78, .70), ("TBX3", .71, .70), ("EN1", .70, .00),
    ("SHOX2", .59, .58), ("MEOX2", .58, .64), ("HOXA6", .54, .65), ("TBX15", .51, .58)]),
 "Haematopoietic": dict(color="#2E5F8A", marker="s", genes=[
    ("LYZ", .69, .00), ("PTPRC", .67, .00), ("CD79A", .59, .00),
    ("CD19", .56, .00), ("MS4A1", .53, .00), ("CD3E", .49, .00)]),
 "Bone / mineral.": dict(color="#1B7B4B", marker="^", genes=[
    ("RUNX2", .93, .35), ("DMP1", .86, .00), ("MEPE", .83, .00), ("SPP1", .81, .00),
    ("SP7", .77, .00), ("ALPL", .77, .00), ("IBSP", .76, .00), ("COL1A1", .63, .00)]),
 "Adipogenic core": dict(color="#E08214", marker="D", genes=[
    ("FABP4", .47, .00), ("CIDEC", .44, .00), ("CFD", .43, .00), ("PLIN4", .41, .00),
    ("PLIN1", .38, .00), ("ADIPOQ", .38, .00), ("PPARG", .32, .00)]),
}
X0, X1 = 0.0, 1.0
axc.axvline(X0, color="#BBBBBB", lw=0.9); axc.axvline(X1, color="#BBBBBB", lw=0.9)
axc.axhspan(-0.014, 0.014, color="#F0F0F0", zorder=0)
span = 0.75
offs = {gn: (i - (len(groups_c) - 1) / 2) / max(len(groups_c) - 1, 1) * span
        for i, gn in enumerate(groups_c)}

def decollide(vals, gap):
    """Order-preserving minimum-gap layout; only ever moves labels up."""
    idx = sorted(range(len(vals)), key=lambda i: vals[i])
    out = list(vals)
    for k in range(1, len(idx)):
        lo, hi = idx[k - 1], idx[k]
        if out[hi] - out[lo] < gap:
            out[hi] = out[lo] + gap
    return out

group_med = {}
for gn, g in groups_c.items():
    o = offs[gn]
    a_vals = [a for _, a, _ in g['genes']]
    lab_y = decollide(a_vals, 0.040)
    for (name, a, b), ly in zip(g['genes'], lab_y):
        axc.plot([X0 + o, X1 + o], [a, b], '-', color=g['color'], lw=1.0, alpha=0.75, zorder=3)
        axc.plot([X0 + o], [a], g['marker'], color=g['color'], ms=3.8, zorder=4)
        axc.plot([X1 + o], [b], g['marker'], color=g['color'], ms=3.8, zorder=4)
        if abs(ly - a) > 0.004:
            axc.plot([X0 + o - 0.055, X0 + o], [ly, a], '-', color=g['color'], lw=0.4,
                     alpha=0.45, zorder=2)
        axc.text(X0 + o - 0.075, ly, name, fontsize=5.4, color=g['color'], ha='right', va='center')
    mm1 = float(np.median(a_vals)); mm3 = float(np.median([x[2] for x in g['genes']]))
    group_med[gn] = (o, mm1, mm3)
    axc.plot([X0 + o, X1 + o], [mm1, mm3], '-', color=g['color'], lw=3.2, alpha=0.95,
             solid_capstyle='round', zorder=5)

# right-hand group labels: the three groups that collapse to zero share y = 0,
# so they are de-collided (with leader lines) to remain legible
order  = sorted(group_med, key=lambda gn: group_med[gn][2])
lab_ys = decollide([group_med[gn][2] for gn in order], 0.085)
for gn, ly in zip(order, lab_ys):
    o, mm1, mm3 = group_med[gn]
    col = groups_c[gn]['color']
    if abs(ly - mm3) > 0.004:
        axc.plot([X1 + o, X1 + o + 0.05], [mm3, ly], '-', color=col, lw=0.4, alpha=0.45, zorder=2)
    axc.text(X1 + o + 0.06, ly, '%s\n%.2f \u2192 %.2f' % (gn, mm1, mm3), fontsize=6.2,
             color=col, va='center', fontweight='bold')
axc.text(-0.88, 0.026, 'variance \u2192 0', fontsize=5.8, color='#888888', va='bottom')
axc.set_xlim(-0.98, 2.08); axc.set_ylim(-0.06, 1.14)
axc.set_xticks([X0, X1])
axc.set_xticklabels(['site + donor\nonly', '+ composition\naxes'], fontsize=7.4)
axc.set_ylabel('Site-attributable variance (per gene)', fontsize=7.6)
axc.set_title('c  Positive control: known positional-identity genes\n'
              '     survive, composition markers do not', fontsize=8.6, loc='left', fontweight='bold')
axc.tick_params(axis='y', labelsize=7.2)

fig.savefig(os.path.join(FIG, 'BMATID_Fig1_site_composition.png'), bbox_inches='tight')
fig.savefig(os.path.join(FIG, 'BMATID_Fig1_site_composition.pdf'), bbox_inches='tight')
plt.close(fig)
print('wrote BMATID_Fig1_site_composition.png/.pdf')
print('panel c group medians:', {gn: (round(np.median([x[1] for x in g['genes']]), 3),
                                     round(np.median([x[2] for x in g['genes']]), 3))
                                 for gn, g in groups_c.items()})
