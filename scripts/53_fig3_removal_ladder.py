## Figure 3 - attributing the site signal to composition: a removal ladder
##             with count-matched random controls   [renumbered 2026-09-22: Fig.2 -> Fig.3]
## a: site-dependent gene count as a function of how many CD45-covarying genes are removed
## b: same-number random removal as specificity control
## Data: results/tables/BMATID_ladder_data.json, written by 19_contamination_filter_limma.R
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

import io, os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

T   = os.path.join(BASE, "results", "tables") + "/"
FIG = os.path.join(BASE, "results", "figures") + "/"
L   = os.path.join(BASE, "logs", "fig3_out.txt")
log = io.open(L, "w", encoding="utf-8", buffering=1)
def say(*a):
    s = " ".join(str(x) for x in a); log.write(s + "\n"); print(s)

# ---- data: the limma-voom ladder, written by 19_contamination_filter_limma.R -------
# model ~ donor + location; the composition fingerprint is co-variation with
# PTPRC/CD45 on the paired in vivo libraries; the equal-number random control is
# fixed by set.seed(42) in that script, so this figure is reproducible end to end.
DAT = T + "BMATID_ladder_data.json"
if not os.path.isfile(DAT):
    raise SystemExit("missing %s -- run 19_contamination_filter_limma.R first" % DAT)
D = json.load(io.open(DAT, encoding="utf-8"))
base_n = D["baseline_genes"]                 # genes entering the model (CPM >= 1 in >= 4)
base_k = D["baseline_site_dependent"]        # site-dependent genes, no removal
remove = D["removed_cd45"]                   # genes removed, one entry per |r(CD45)| cut
kept   = D["site_dependent_after"]           # site-dependent genes left after removal
rnd_x  = D["random_removed"]                 # equal-number random control: genes removed
rnd_y  = D["random_site_dependent"]          # equal-number random control: genes left
assert len(remove) == len(kept), "ladder arrays are not paired"
assert len(rnd_x) == len(rnd_y), "random-control arrays are not paired"

say("ladder (removed, site genes):", list(zip(remove, kept)))
say("random controls:", list(zip(rnd_x, rnd_y)))
say("relative to baseline: CD45", [round(k / base_k * 100, 1) for k in kept])
say("relative to baseline: random", [round(r / base_k * 100, 1) for r in rnd_y])

fig = plt.figure(figsize=(11.0, 4.3), dpi=300)
gs  = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.30,
                       left=0.085, right=0.975, top=0.86, bottom=0.15)

C_SITE = "#B2182B"      # red-ish (Chinese convention: increase = red)
C_RAND = "#4D4D4D"
C_BASE = "#2166AC"

# ---------- panel a ----------
ax = fig.add_subplot(gs[0, 0])
ax.axhline(base_k, ls=":", lw=1.3, color=C_BASE, zorder=1)
ax.text(base_n * 0.02 + 150, base_k + 90, "all genes: %d site-dependent" % base_k,
        color=C_BASE, fontsize=8.4, va="bottom")

ax.plot(remove, kept, "-o", color=C_SITE, lw=2.0, ms=7, mfc="white", mew=1.8,
        label="remove genes covarying with CD45 (composition fingerprint)", zorder=3)
for x, y in zip(remove, kept):
    ax.annotate("%d" % y, (x, y), textcoords="offset points", xytext=(0, -16),
                ha="center", fontsize=8.2, color=C_SITE, fontweight="bold")

ax.plot(rnd_x, rnd_y, "s--", color=C_RAND, lw=1.7, ms=7, mfc="white", mew=1.6,
        label="remove the same number of genes at random", zorder=3)
for x, y in zip(rnd_x, rnd_y):
    ax.annotate("%d" % y, (x, y), textcoords="offset points", xytext=(0, 9),
                ha="center", fontsize=8.2, color=C_RAND, fontweight="bold")

ax.set_xlabel("number of genes removed from the model", fontsize=9.5)
ax.set_ylabel("genes with a site effect\n(FDR < 0.05, limma-voom)", fontsize=9.5)
ax.set_title("a   Attribution of the in vivo site signal", fontsize=10.2, loc="left",
             fontweight="bold", pad=8)
ax.set_xlim(-350, 7600); ax.set_ylim(-140, 3200)
ax.tick_params(labelsize=8.6)
ax.spines[["top", "right"]].set_visible(False)
ax.legend(fontsize=7.9, frameon=False, loc="upper right")

# ---------- panel b ----------
ax2 = fig.add_subplot(gs[0, 1])
lbl = ["CD45\ncovarying", "random\n(equal n)", "CD45\ncovarying", "random\n(equal n)"]
vals = [kept[2] / base_k * 100, rnd_y[0] / base_k * 100,
        kept[4] / base_k * 100, rnd_y[1] / base_k * 100]
cols = [C_SITE, C_RAND, C_SITE, C_RAND]
xs = np.arange(4)
bars = ax2.bar(xs, vals, color=cols, width=0.62, alpha=0.9, zorder=2)
for b, v in zip(bars, vals):
    ax2.text(b.get_x() + b.get_width() / 2, v + 2.5, "%.0f%%" % v,
             ha="center", fontsize=9.0, fontweight="bold")
for i, (a, b) in enumerate([(0, 1), (2, 3)]):
    ax2.plot([xs[a], xs[b]], [max(vals[a], vals[b]) + 14] * 2, "-", color="#333333", lw=1.0)
    ax2.text((xs[a] + xs[b]) / 2, max(vals[a], vals[b]) + 16, "%.2f×" % (vals[b] / vals[a]),
             ha="center", fontsize=8.0)
ax2.set_xticks(xs); ax2.set_xticklabels(lbl, fontsize=8.0)
ax2.set_ylabel("site-dependent genes retained\n(% of the 2,712 baseline)", fontsize=9.5)
ax2.set_title("b   Same number of genes, different genes", fontsize=10.2, loc="left",
              fontweight="bold", pad=8)
ax2.set_ylim(0, 108)
ax2.tick_params(labelsize=8.6)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(0.5, -0.30, "removing 2,975 / 7,070 genes", transform=ax2.transAxes,
         ha="center", fontsize=7.6, color="#555555")

fig.savefig(FIG + "BMATID_Fig3_removal_ladder.png", dpi=300, facecolor="white")
fig.savefig(FIG + "BMATID_Fig3_removal_ladder.pdf", facecolor="white")
say("\nwrote BMATID_Fig3_removal_ladder.png/.pdf")

say("DONE")
log.close()
