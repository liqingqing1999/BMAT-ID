# -*- coding: utf-8 -*-
"""Summarise the three variance-decomposition fits into one json.

Reads the per-gene site-attributable variance produced by
21_varpart_main_models.R and writes results/tables/BMATID_varpart_summary.json,
the machine-readable record of the headline quantities quoted in the paper
(mean / median site variance, fraction of genes collapsing to zero, survival of
the strongest site genes).

The headline zero fraction uses the manuscript's rule -- the site term is exactly
zero, not merely below a numerical tolerance. The 1e-9 variant is stored alongside
it because the Figure 1 legend shows the two tiers side by side.
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

T = os.path.join(BASE, "results", "tables")


def load(fn):
    with open(os.path.join(T, fn), newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def gv(rows, key="localization"):
    out = {}
    for r in rows:
        try:
            out[r["ensg"]] = float(r.get(key, ""))
        except (TypeError, ValueError):
            pass
    return out


M1 = gv(load("BMATID_varpart_invivo.csv"))        # in vivo, site only
M2 = gv(load("BMATID_varpart_invitro.csv"))       # in vitro, site only
M3 = gv(load("BMATID_varpart_invivo_comp.csv"))   # in vivo + composition axes

l1 = np.array(list(M1.values()))
l2 = np.array(list(M2.values()))
l3 = np.array(list(M3.values()))

hi1 = [g for g in M1 if M1[g] > .5]
surv = sum(1 for g in hi1 if g in M3 and M3[g] > .5)

out = dict(
    mean_pct_loc_vivo=float(100 * l1.mean()),
    mean_pct_loc_vitro=float(100 * l2.mean()),
    mean_pct_loc_comp=float(100 * l3.mean()),
    median_pct_loc_vivo=float(100 * np.median(l1)),
    median_pct_loc_comp=float(100 * np.median(l3)),
    frac_loc_zero_comp=float(np.mean(l3 == 0)),            # exact zero: the reported rule
    n_loc_zero_comp=int((l3 == 0).sum()),
    n_fitted_genes=int(len(l3)),                           # the fitted universe (20,041)
    frac_loc_le1e9_comp=float(np.mean(l3 <= 1e-9)),        # legend-only comparison
    n_loc_gt50_vivo=int((l1 > .5).sum()),
    n_loc_gt50_comp=int(surv),
    survival_of_top_site_genes=float(surv / len(hi1)),
    n_loc_gt80_vivo=int((l1 > .8).sum()),
    n_loc_gt80_comp=int(sum(1 for g in M1 if M1[g] > .8 and M3.get(g, 0) > .8)),
)
with open(os.path.join(T, "BMATID_varpart_summary.json"), "w", encoding="utf-8") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
