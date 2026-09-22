## Figure 2 - which programmes survive adipogenic culture   [renumbered 2026-09-22: Fig.3 -> Fig.2]
## a: in vitro / in vivo effect-size retention by functional module
## b: in vivo vs in vitro median F, with the identity line
## Data: results/tables/BMATID_module_retention.csv, written by 40_functional_modules.py
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

import io
import csv as _csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = os.path.join(BASE, "results", "figures") + "/"
L   = os.path.join(BASE, "logs", "fig2_out.txt")
log = io.open(L, "w", encoding="utf-8", buffering=1)
def say(*a):
    s = " ".join(str(x) for x in a); log.write(s + "\n"); print(s)

# ---- data: results/tables/BMATID_module_retention.csv, written by
# 40_functional_modules.py. Nothing below is transcribed by hand: the per-module
# medians and retention are read from the same run the manuscript quotes.
T   = os.path.join(BASE, "results", "tables") + "/"
RET = T + "BMATID_module_retention.csv"
if not os.path.isfile(RET):
    raise SystemExit("missing %s -- run 40_functional_modules.py first" % RET)
_ret = {r["module"]: r for r in _csv.DictReader(io.open(RET, encoding="utf-8"))}
# display name per module id of 40_*; the plot order is the module definition order
ORDER = ["C_造血/浆细胞(污染指纹)", "A_骨/矿化", "F_MSC/干性",
         "B_造血龛", "D_成脂核心", "E_内皮"]
LABEL = {"C_造血/浆细胞(污染指纹)": "Haematopoietic /\nplasma cell",
         "A_骨/矿化":              "Bone /\nmineralisation",
         "F_MSC/干性":             "MSC /\nstemness",
         "B_造血龛":               "Haematopoietic\nniche",
         "D_成脂核心":             "Adipogenic\ncore",
         "E_内皮":                 "Endothelial"}
# module, n genes, median F in vivo, median F in vitro, retention %
mods = [("%s (%d)" % (LABEL[k], int(_ret[k]["n_genes"])),
         int(_ret[k]["n_genes"]),
         float(_ret[k]["median_F_in_vivo"]),
         float(_ret[k]["median_F_in_vitro"]),
         float(_ret[k]["retention_pct"])) for k in ORDER]
say("module retention table:")
for m in mods:
    say("  %-34s n=%2d  F_vivo=%6.2f  F_vitro=%5.2f  retention=%5.1f%%" %
        (m[0].replace("\n", " "), m[1], m[2], m[3], m[4]))

C_ADIPO = "#B2182B"
C_OTHER = "#4A6FA5"
C_TXT   = "#222222"

fig = plt.figure(figsize=(11.2, 4.6), dpi=300)
gs  = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0], wspace=0.34,
                       left=0.128, right=0.975, top=0.85, bottom=0.20)

# ---------- panel a: retention ----------
ax = fig.add_subplot(gs[0, 0])
order = sorted(mods, key=lambda m: m[4])
names = [m[0] for m in order]
vals  = [m[4] for m in order]
cols  = [C_ADIPO if "Adipo" in n else C_OTHER for n in names]
y = np.arange(len(names))
bars = ax.barh(y, vals, color=cols, height=0.62, alpha=0.92, zorder=2)
for b, v in zip(bars, vals):
    ax.text(v + 1.8, b.get_y() + b.get_height() / 2, "%.1f%%" % v,
            va="center", fontsize=8.8, fontweight="bold")
ax.axvline(100, ls=":", lw=1.1, color="#888888")
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8.1)
ax.set_xlim(0, 112)
ax.set_xlabel("in vitro / in vivo effect-size retention (%)", fontsize=9.5)
ax.set_title("a   What adipogenic culture preserves", fontsize=10.2, loc="left",
             fontweight="bold", pad=8)
ax.tick_params(axis="x", labelsize=8.6)
ax.spines[["top", "right"]].set_visible(False)
ax.text(0.985, 0.06, "identity line = no loss", transform=ax.transAxes,
        ha="right", fontsize=7.6, color="#777777")

# ---------- panel b: F vivo vs F vitro ----------
ax2 = fig.add_subplot(gs[0, 1])
lim = [0, 32]
ax2.plot(lim, lim, ls="--", lw=1.2, color="#999999", zorder=1)
ax2.fill_between(lim, lim, [lim[1]] * 2, color="#F2F2F2", zorder=0)
fx = [m[2] for m in mods]; fy = [m[3] for m in mods]
# per-module label placement: four modules sit within 11-13 on the x axis, and their
# default labels collided with each other and with the right edge of the panel
OFF = {
    "Adipogenic core":        (8, 2, "left"),
    "Endothelial":            (10, 5, "left"),
    "MSC / stemness":         (-6, -11, "right"),
    "Bone / mineralisation":  (7, 5, "left"),
    "Haematopoietic niche":   (-6, 8, "right"),
    "Haematopoietic / plasma cell": (6, -11, "left"),
}
for m in mods:
    isA = "Adipo" in m[0]
    ax2.scatter(m[2], m[3], s=76, color=C_ADIPO if isA else C_OTHER,
                edgecolor="white", lw=1.2, zorder=3)
    key = m[0].replace("\n", " ").split(" (")[0]
    dx, dy, ha = OFF.get(key, (6, -7, "left"))
    ax2.annotate(key, (m[2], m[3]),
                 textcoords="offset points", xytext=(dx, dy), ha=ha, fontsize=7.4,
                 color=C_TXT if isA else "#444444",
                 fontweight="bold" if isA else "normal")
ax2.set_xlim(*lim); ax2.set_ylim(-0.5, lim[1])
ax2.set_xlabel("median site-effect F, in vivo", fontsize=9.5)
ax2.set_ylabel("median site-effect F, in vitro", fontsize=9.5)
ax2.set_title("b   Site effect before and after culture", fontsize=10.2, loc="left",
              fontweight="bold", pad=8)
ax2.tick_params(labelsize=8.6)
ax2.spines[["top", "right"]].set_visible(False)
ax2.text(0.97, 0.93, "above the line =\npreserved in vitro", transform=ax2.transAxes,
         ha="right", va="top", fontsize=7.6, color="#666666")

fig.savefig(FIG + "BMATID_Fig2_module_retention.png", dpi=300, facecolor="white")
fig.savefig(FIG + "BMATID_Fig2_module_retention.pdf", facecolor="white")
say("\nwrote BMATID_Fig2_module_retention.png/.pdf")
say("DONE")
log.close()
