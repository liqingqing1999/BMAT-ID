# -*- coding: utf-8 -*-
"""Figure 5 -- simulation-based independent validation.

Panel a: three regimes with known truth (A composition only, B intrinsic only,
         C both) x three estimators (naive / composition-adjusted / oracle),
         against the real data on the same gene pools.
Panel b: composition-strength sweep at a fixed intrinsic load -- the naive
         estimate tracks the confounder, the adjusted estimate does not.

Input : results/tables/BMATID_sim_validation_summary.csv
        results/tables/BMATID_sim_validation_realref.csv
        (produced by 60_sim_validation.R and 61_sim_real_ref.R)
Output: results/figures/BMATID_Fig5_simulation_validation.{png,pdf}
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

import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG = os.path.join(BASE, "results", "figures")
T = os.path.join(BASE, "results", "tables")
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'xtick.major.width': 0.8, 'ytick.major.width': 0.8, 'savefig.dpi': 400,
    'figure.facecolor': 'white', 'axes.facecolor': 'white', 'pdf.fonttype': 42,
})
C_NAIVE, C_ADJ, C_ORAC, C_TRUE = '#C44E52', '#2F4B7C', '#8C8C8C', '#55A868'


def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))


S = {}
for r in load('BMATID_sim_validation_summary.csv'):
    S[(r['regime'], r['metric'])] = float(r['value']) if r['value'] not in ('', 'NA') else np.nan


def g(regime, metric):
    return S.get((regime, metric), np.nan)


REAL = load('BMATID_sim_validation_realref.csv')
real_naive = float(np.mean([float(r['real_naive_mean']) for r in REAL]))
real_adj = float(np.mean([float(r['real_adj_mean']) for r in REAL]))
print('real on the simulation pools: naive %.4f  adjusted %.4f' % (real_naive, real_adj))

REG = [('A', 'Composition only\n(no intrinsic effect)'),
       ('B', 'Intrinsic only\n(composition off)'),
       ('C', 'Both\n(bulk in vivo analogue)')]
naive = [g(r, 'naive_mean') for r, _ in REG]
adj = [g(r, 'adj_mean') for r, _ in REG]
orac = [g(r, 'orac_mean') for r, _ in REG]
print('regimes  naive', ['%.4f' % v for v in naive])
print('regimes  adj  ', ['%.4f' % v for v in adj])
print('regimes  orac ', ['%.4f' % v for v in orac])

SW = ['0.0', '0.5', '1.0', '1.5', '2.0']
xs = np.array([float(v) for v in SW])
sw_naive = np.array([g('sweep_comp' + v, 'naive_mean') for v in SW])
sw_adj = np.array([g('sweep_comp' + v, 'adj_mean') for v in SW])
sw_orac = np.array([g('sweep_comp' + v, 'orac_mean') for v in SW])
sw_true = np.array([g('sweep_comp' + v, 'true_intr_mean') for v in SW])
print('sweep    naive', ['%.4f' % v for v in sw_naive])
print('sweep    adj  ', ['%.4f' % v for v in sw_adj])

fig = plt.figure(figsize=(7.1, 2.85))
gs = fig.add_gridspec(1, 2, width_ratios=[1.16, 1.0], wspace=0.30)

# ------------------------------------------------------------------ a: regimes
ax = fig.add_subplot(gs[0, 0])
x = np.arange(len(REG))
w = 0.26
bars = [('Naive\n(site only)', naive, C_NAIVE),
        ('Adjusted\n(+ composition)', adj, C_ADJ),
        ('Oracle\n(true axes)', orac, C_ORAC)]
for k, (lab, vals, col) in enumerate(bars):
    ax.bar(x + (k - 1) * w, vals, w, color=col, label=lab, edgecolor='white', linewidth=0.5, zorder=3)
ax.axhline(real_naive, color=C_NAIVE, ls=(0, (4, 2)), lw=1.0, zorder=4)
ax.axhline(real_adj, color=C_ADJ, ls=(0, (4, 2)), lw=1.0, zorder=4)
ax.text(2.46, real_naive, 'real, naive', color=C_NAIVE, fontsize=7.0, va='bottom', ha='right')
ax.text(2.46, real_adj, 'real, adjusted', color=C_ADJ, fontsize=7.0, va='bottom', ha='right')
ax.set_xticks(x)
ax.set_xticklabels([lab for _, lab in REG], fontsize=7.4)
ax.set_ylabel('Site-attributable variance\n(mean over genes, %)', fontsize=7.8)
ax.set_ylim(0, 0.30)
ax.set_yticks(np.arange(0, 0.31, 0.05))
ax.set_yticklabels(['%d' % round(v * 100) for v in np.arange(0, 0.31, 0.05)], fontsize=7.4)
ax.legend(frameon=False, fontsize=6.8, loc='upper left', ncol=1, handlelength=1.1,
          handletextpad=0.5, labelspacing=0.35)
ax.set_title('Known truth: only composition produces a site term', fontsize=8.0, pad=5)
for k, (lab, vals, col) in enumerate(bars):
    for i, v in enumerate(vals):
        if np.isfinite(v):
            ax.text(x[i] + (k - 1) * w, v + 0.004, '%.1f' % (v * 100), ha='center',
                    va='bottom', fontsize=6.2, color='#333333')
ax.text(-0.42, 0.288, 'a', fontsize=10, fontweight='bold', va='top')

# -------------------------------------------------------------------- b: sweep
ax = fig.add_subplot(gs[0, 1])
ax.plot(xs, sw_true, color=C_TRUE, ls=(0, (1, 1.2)), lw=1.3, marker='o', ms=3.2,
        label='True intrinsic load', zorder=3)
ax.plot(xs, sw_naive, color=C_NAIVE, lw=1.5, marker='s', ms=3.2, label='Naive (site only)', zorder=4)
ax.plot(xs, sw_adj, color=C_ADJ, lw=1.5, marker='^', ms=3.4, label='Adjusted (+ composition)', zorder=5)
ax.plot(xs, sw_orac, color=C_ORAC, lw=1.0, ls=':', marker='d', ms=2.8,
        label='Oracle (true axes)', zorder=3)
ax.set_xlabel('Composition strength (scale factor on the real loadings)', fontsize=7.8)
ax.set_ylabel('Site-attributable variance\n(mean over genes, %)', fontsize=7.8)
ax.set_xticks(xs)
ax.set_xticklabels(['%g' % v for v in xs], fontsize=7.4)
ax.set_ylim(0, 0.42)
ax.set_yticks(np.arange(0, 0.41, 0.10))
ax.set_yticklabels(['%d' % round(v * 100) for v in np.arange(0, 0.41, 0.10)], fontsize=7.4)
ax.legend(frameon=False, fontsize=6.8, loc='upper left', handlelength=1.4,
          handletextpad=0.5, labelspacing=0.35)
ax.set_title('The naive estimate tracks the confounder', fontsize=8.0, pad=5)
ax.annotate('', xy=(2.0, sw_naive[-1]), xytext=(2.0, sw_adj[-1]),
            arrowprops=dict(arrowstyle='<->', lw=0.7, color='#555555'))
ax.text(1.93, (sw_naive[-1] + sw_adj[-1]) / 2, '%.1f-fold' % (sw_naive[-1] / sw_adj[-1]),
        fontsize=6.6, color='#555555', ha='right', va='center', rotation=90)
ax.text(-0.22, 0.405, 'b', fontsize=10, fontweight='bold', va='top')

fig.savefig(os.path.join(FIG, 'BMATID_Fig5_simulation_validation.png'),
            bbox_inches='tight', facecolor='white')
fig.savefig(os.path.join(FIG, 'BMATID_Fig5_simulation_validation.pdf'),
            bbox_inches='tight', facecolor='white')
print('[F5] simulation figure written')
