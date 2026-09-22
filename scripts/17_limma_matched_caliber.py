# -*- coding: utf-8 -*-
"""口径对齐：用原文阈值（FDR<0.05 & |log2FC|>2）复算 primary vs diff 的 DEG"""
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

import gzip, os
p = os.path.join(BASE, "data", "raw", "GSE291355_counts.tsv.gz")
LOG = os.path.join(BASE, "logs", "caliber_head.txt")
os.makedirs(os.path.dirname(LOG), exist_ok=True)
buf = []
def w(s=''):
    buf.append(str(s)); print(s)

with gzip.open(p, 'rt', encoding='utf-8', errors='ignore') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    w('ncol = %d' % len(hdr))
    w('header:')
    for i, h in enumerate(hdr):
        w('  [%d] %r' % (i, h))
    w('')
    w('first 3 data rows (first 6 cols):')
    for k in range(3):
        ln = fh.readline().rstrip('\n').split('\t')
        w('  ' + ' | '.join('%s' % x for x in ln[:6]))

with open(LOG, 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(buf))
print('\n[log]', LOG)
