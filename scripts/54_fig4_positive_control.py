# -*- coding: utf-8 -*-
"""Figure 4 -- what survives composition correction (2 panels)

   a  POSITIVE CONTROL: known positional-identity genes keep their site variance
      after composition adjustment, while classical haematopoietic,
      bone/mineralisation and adipogenic markers lose it
   b  survival of the strong site genes, by site-variance tier

   WHY THIS FILE EXISTS
   Up to 2026-09-22 these two panels were panels c and d of a four-panel
   Figure 1.  They answer a different question from panels a/b -- "what is left
   after correction?" -- which the manuscript argues last (Results 3.5), and
   figure numbering has to follow the order of first citation.  The canvas was
   therefore split and the two panels re-based c/d -> a/b.  Only the panel
   letters are printed on the canvas, so this split is the single part of the
   renumbering that needed a redraw; the Fig. 2 / Fig. 3 exchange is a rename.

Reads : results/tables/BMATID_varpart_invivo.csv, BMATID_varpart_invivo_comp.csv
        (21_varpart_main_models.R) and BMATID_varpart_allgenes_merged.csv
        (50_fig1_panel_data.py).  The panel-a gene values are read from that
        table, never typed in here.
Output: results/figures/BMATID_Fig4_positive_control.{png,pdf}
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
import os, csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

FIG  = os.path.join(BASE, 'results/figures')
T    = os.path.join(BASE, 'results/tables')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8.5, 'axes.linewidth': 0.8,
                     'axes.spines.top': False, 'axes.spines.right': False,
                     'savefig.dpi': 400, 'figure.facecolor': 'white', 'axes.facecolor': 'white'})
C_COMP = '#8C8C8C'

def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

def gvd(d, key='localization'):
    out = {}
    for k, r in d.items():
        try: out[k] = float(r.get(key, ''))
        except (TypeError, ValueError): pass
    return out

d1 = {r['ensg']: r for r in load('BMATID_varpart_invivo.csv')}
d3 = {r['ensg']: r for r in load('BMATID_varpart_invivo_comp.csv')}
M1, M3 = gvd(d1), gvd(d3)

fig = plt.figure(figsize=(5.9, 5.0))
# Print-size driven: the figure is placed at the 5.83 in text width, so the canvas
# must be close to that or every point size shrinks with it.  At 9.8 in the 33 gene
# labels of panel a printed at ~3.5 pt; at 6.4 in they print at ~5.5 pt.  The canvas
# is therefore narrow and tall rather than wide: panel a needs height (33 stacked
# labels), not width.
BOX = {
    'a': [0.076, 0.085, 0.600, 0.830],   # positive control needs the widest box
    'b': [0.735, 0.085, 0.255, 0.830],   # survival by tier
}

# ---------------- b  survival by tier (drawn first so the shared y-log axis is set) -------
ax = fig.add_axes(BOX['b'])
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
ax.set_title('b  Only 11% of the strongest\nsite genes survive', fontsize=8.6,
             loc='left', fontweight='bold')

# ---------------- a  positive control ------------------------------------------------
axc = fig.add_axes(BOX['a'])
# Panel a is a positive control over four marker groups.  The per-gene values
# are the site-attributable variance in vivo alone (M1) and after composition
# adjustment (M3), read from the merged table written by 50_fig1_panel_data.py (itself
# assembled from the three variance tables) -- none of it is transcribed into
# this source.  M1/M3 are drawn at 2 decimals, so the same rounding is applied.
GROUPS_A = {
 "Axis-patterning TF": dict(color="#C1272D", marker="o", genes=[
    "TBX5", "RSPO1", "DLX5", "BARX1", "TBX18", "PYGO1", "TBX3", "EN1",
    "SHOX2", "MEOX2", "HOXA6", "TBX15"]),
 "Haematopoietic": dict(color="#2E5F8A", marker="s", genes=[
    "LYZ", "PTPRC", "CD79A", "CD19", "MS4A1", "CD3E"]),
 "Bone / mineral.": dict(color="#1B7B4B", marker="^", genes=[
    "RUNX2", "DMP1", "MEPE", "SPP1", "SP7", "ALPL", "IBSP", "COL1A1"]),
 "Adipogenic core": dict(color="#E08214", marker="D", genes=[
    "FABP4", "CIDEC", "CFD", "PLIN4", "PLIN1", "ADIPOQ", "PPARG"]),
}
MER = os.path.join(T, 'BMATID_varpart_allgenes_merged.csv')
if not os.path.isfile(MER):
    raise SystemExit("missing %s -- run 50_fig1_panel_data.py first" % MER)
_val = {}
for _r in load('BMATID_varpart_allgenes_merged.csv'):
    if _r.get('symbol'):
        _val[_r['symbol']] = (round(float(_r['M1']), 2), round(float(_r['M3']), 2))
groups_c = {}
for _gn, _g in GROUPS_A.items():
    _genes = []
    for _sym in _g['genes']:
        if _sym not in _val:
            raise SystemExit("positive-control gene %s is not in %s"
                             % (_sym, os.path.basename(MER)))
        _genes.append((_sym, _val[_sym][0], _val[_sym][1]))
    groups_c[_gn] = dict(color=_g['color'], marker=_g['marker'], genes=_genes)
print('panel a: %d genes read from BMATID_varpart_allgenes_merged.csv'
      % sum(len(g['genes']) for g in groups_c.values()))
X0, X1 = 0.0, 1.0
axc.axvline(X0, color="#BBBBBB", lw=0.9); axc.axvline(X1, color="#BBBBBB", lw=0.9)
axc.axhspan(-0.014, 0.014, color="#F0F0F0", zorder=0)
span = 0.92
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
    lab_y = decollide(a_vals, 0.048)
    for (name, a, b), ly in zip(g['genes'], lab_y):
        axc.plot([X0 + o, X1 + o], [a, b], '-', color=g['color'], lw=1.0, alpha=0.75, zorder=3)
        axc.plot([X0 + o], [a], g['marker'], color=g['color'], ms=3.8, zorder=4)
        axc.plot([X1 + o], [b], g['marker'], color=g['color'], ms=3.8, zorder=4)
        if abs(ly - a) > 0.004:
            axc.plot([X0 + o - 0.055, X0 + o], [ly, a], '-', color=g['color'], lw=0.4,
                     alpha=0.45, zorder=2)
        # white halo: with four groups the labels of the right-hand group sit over
        # the data column of the group to its left, so the gene names need a
        # background to stay readable (pre-existing crowding, not a data change)
        axc.text(X0 + o - 0.075, ly, name, fontsize=5.7, color=g['color'], ha='right',
                 va='center', bbox=dict(facecolor='white', edgecolor='none', pad=0.12,
                                        alpha=0.88))
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
axc.set_xlim(-1.06, 2.10); axc.set_ylim(-0.06, 1.16)
axc.set_xticks([X0, X1])
axc.set_xticklabels(['site + donor\nonly', '+ composition\naxes'], fontsize=7.4)
axc.set_ylabel('Site-attributable variance (per gene)', fontsize=7.6)
axc.set_title('a  Positive control: known positional-identity genes\n'
              '     survive, composition markers do not', fontsize=8.6, loc='left', fontweight='bold')
axc.tick_params(axis='y', labelsize=7.2)

fig.savefig(os.path.join(FIG, 'BMATID_Fig4_positive_control.png'), bbox_inches='tight')
fig.savefig(os.path.join(FIG, 'BMATID_Fig4_positive_control.pdf'), bbox_inches='tight')
plt.close(fig)
print('wrote BMATID_Fig4_positive_control.png/.pdf')
print('panel a group medians:', {gn: (round(np.median([x[1] for x in g['genes']]), 3),
                                      round(np.median([x[2] for x in g['genes']]), 3))
                                 for gn, g in groups_c.items()})
print('panel b tier counts: in vivo alone', n_m1, '| after composition', n_m3)
