# -*- coding: utf-8 -*-
"""Figure 1 -- how much of the transcriptome is a property of the anatomical site (2 panels)

   a  variance composition across four mixed models (incl. reference-anchored)
   b  ECDF of per-gene site-attributable variance

   HISTORY -- WHY THIS FILE IS NO LONGER CALLED "site_composition"
   Until 2026-09-22 this script drew a four-panel Figure 1.  Its two extra panels
   answered a different question -- "what is left after correction?" -- which the
   manuscript argues last, and figure numbering has to follow the order of first
   citation.  The canvas was therefore split:

       Fig. 1 = a + b          -> Results 3.1 / 3.2   (this file)
       Fig. 4 = the old c + d  -> Results 3.5         (54_fig4_positive_control.py)

   The old Figures 2 and 3 were exchanged in the same pass: the module panel is
   now Figure 2 (52_fig2_module_retention.py) and the removal ladder Figure 3
   (53_fig3_removal_ladder.py).  That swap is a pure file rename -- only the
   panel letters are printed on the canvas, never the figure number.

Input : results/tables/BMATID_varpart_invivo.csv, BMATID_varpart_invitro.csv,
        BMATID_varpart_invivo_comp.csv  (21_varpart_main_models.R)
        BMATID_varpart_M7_reference.csv (32_varpart_reference_anchored.R)
Output: results/figures/BMATID_Fig1_site_variance.{png,pdf}
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

fig = plt.figure(figsize=(7.4, 3.6))
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.36], wspace=0.26,
                      left=0.085, right=0.985, top=0.795, bottom=0.155)

# ---------------- a ----------------
ax = fig.add_subplot(gs[0, 0])
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
ax.set_xticks(range(N)); ax.set_xticklabels(names, fontsize=6.8)
ax.set_ylabel('Variance component\n(mean % of gene variance)', fontsize=7.6)
ax.set_ylim(0, 112); ax.set_xlim(-0.6, N - 0.4); ax.set_yticks([0, 20, 40, 60, 80, 100])
fig.legend(handles=handles + [rb], fontsize=6.2, frameon=False, loc='upper center',
           bbox_to_anchor=(0.5, 0.988), ncol=5, columnspacing=0.9, handlelength=1.0,
           handletextpad=0.4)
ax.set_title('a  Variance composition', fontsize=8.6, loc='left', fontweight='bold')

# ---------------- b ----------------
ax = fig.add_subplot(gs[0, 1])
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

fig.savefig(os.path.join(FIG, 'BMATID_Fig1_site_variance.png'), bbox_inches='tight')
fig.savefig(os.path.join(FIG, 'BMATID_Fig1_site_variance.pdf'), bbox_inches='tight')
plt.close(fig)
print('wrote BMATID_Fig1_site_variance.png/.pdf')
