# -*- coding: utf-8 -*-
"""Supplementary Figure S1 -- per-library reference-anchored composition scores.

Left : the four reference signatures (immune, erythroid, skeletal, vasculature)
       across the 24 libraries, in vivo (n = 12) versus in vitro (n = 12).
Right: donor-matched in vivo versus in vitro score for every signature.

Input : results/tables/BMATID_deconv_L1_scores.csv
        (produced by 31_reference_anchored_scores.py)
Output: results/figures/BMATID_FigS1_reference_scores.{png,pdf}
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

import os, csv, io
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TB  = os.path.join(BASE, "results", "tables")
FIG = os.path.join(BASE, "results", "figures")

rows = list(csv.DictReader(io.open(os.path.join(TB, 'BMATID_deconv_L1_scores.csv'), encoding='utf-8')))
assert len(rows) == 24, len(rows)
sig = ['immune', 'erythroid', 'skeletal', 'vasculature']
col = {'immune': '#2F4B7C', 'erythroid': '#C44E52', 'skeletal': '#55A868', 'vasculature': '#8172B2'}
order = sorted(r[''] for r in rows)
# 排序：先体内（primary）按部位 EPI/META/SCAT，再体外（diff）
site_rank = {'epi': 0, 'meta': 1, 'scat': 2}
prim = sorted([r for r in rows if r[''].endswith('primary_adip')],
              key=lambda r: (site_rank[r[''].split('_')[1]], r[''].split('_')[0]))
vit = sorted([r for r in rows if r[''].endswith('diff_adip')],
             key=lambda r: (site_rank[r[''].split('_')[1]], r[''].split('_')[0]))
seq = prim + vit
lbl = [r[''].replace('_primary_adip', ' · in vivo').replace('_diff_adip', ' · in vitro')
       .replace('_', ' ').replace('epi', 'EPI').replace('meta', 'META').replace('scat', 'SCAT')
       for r in seq]

fig, axes = plt.subplots(1, 2, figsize=(13.2, 4.4), gridspec_kw=dict(width_ratios=[1.55, 1.0], wspace=0.28))
ax = axes[0]
x = np.arange(len(seq)); w = 0.2
for k, s in enumerate(sig):
    v = [float(r[s]) for r in seq]
    ax.bar(x + (k - 1.5) * w, v, w, label=s, color=col[s], edgecolor='white', linewidth=0.4, zorder=3)
ax.axvline(11.5, color='#333', lw=1.0, ls='--')
ax.text(5.5, ax.get_ylim()[1] * 0.94, 'in vivo (n = 12)', ha='center', fontsize=8.4, color='#333')
ax.text(17.5, ax.get_ylim()[1] * 0.94, 'in vitro (n = 12)', ha='center', fontsize=8.4, color='#333')
ax.set_xticks(x); ax.set_xticklabels(lbl, rotation=90, fontsize=5.6)
ax.set_ylabel('Reference-anchored score\n(mean z x sqrt(n))', fontsize=7.8)
ax.legend(fontsize=7, frameon=False, ncol=4, loc='upper center', bbox_to_anchor=(0.5, -0.42))
ax.set_title('All four reference signatures collapse in vitro', fontsize=9, loc='left', fontweight='bold')
ax.axhline(0, color='#999', lw=0.7)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)

ax2 = axes[1]
for s in sig:
    pv = [float(r[s]) for r in prim]; vv = [float(r[s]) for r in vit]
    ax2.scatter(pv, vv, s=26, color=col[s], label=s, edgecolor='white', linewidth=0.5, zorder=3)
lo = min(float(r[s]) for r in rows for s in sig); hi = max(float(r[s]) for r in rows for s in sig)
ax2.plot([lo, hi], [lo, hi], ls=':', color='#999', lw=1.0, zorder=1)
ax2.axhline(0, color='#BBB', lw=0.7); ax2.axvline(0, color='#BBB', lw=0.7)
ax2.set_xlabel('in vivo score', fontsize=7.8); ax2.set_ylabel('in vitro score', fontsize=7.8)
ax2.set_title('Donor-matched in vivo vs in vitro', fontsize=9, loc='left', fontweight='bold')
ax2.legend(fontsize=6.6, frameon=False, loc='lower right')
for sp in ('top', 'right'):
    ax2.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'BMATID_FigS1_reference_scores.png'), dpi=300, bbox_inches='tight')
fig.savefig(os.path.join(FIG, 'BMATID_FigS1_reference_scores.pdf'), bbox_inches='tight')
plt.close(fig)
print('[S1] supplementary figure written')
