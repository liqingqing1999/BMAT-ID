# -*- coding: utf-8 -*-
"""解析三轴标记基因 symbol -> Ensembl gene id（权威映射，避免手写错号）"""
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

import json, time, urllib.request, urllib.parse, ssl, os
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

AXES = {
 'Hemato': ['PTPRC','CD14','CD68','LYZ','CD3E','MS4A1','MZB1','JCHAIN','CSF1R','CD19','CD79A'],
 'Bone'  : ['BGLAP','SP7','ALPL','IBSP','SPP1','MEPE','DMP1','RUNX2','PHEX','COL1A1','TMEM119','ENPP1'],
 'Adipo' : ['ADIPOQ','PLIN1','PLIN4','FABP4','GPAM','CFD','LEP','LPL','PLIN2','CIDEC'],
}
syms = sorted({s for v in AXES.values() for s in v})
out = {}
for i in range(0, len(syms), 60):
    ch = syms[i:i+60]
    d = urllib.parse.urlencode({'q': ','.join(ch), 'scopes': 'symbol', 'fields': 'ensembl.gene,symbol',
                                'species': 'human'}).encode()
    for a in range(4):
        try:
            r = urllib.request.Request('https://mygene.info/v3/query', data=d,
                                       headers={'Content-Type': 'application/x-www-form-urlencoded'})
            res = json.loads(urllib.request.urlopen(r, context=ctx, timeout=90).read())
            for h in (res if isinstance(res, list) else res.get('hits', [])):
                q = h.get('query')
                if not q or h.get('notfound'):
                    continue
                e = h.get('ensembl')
                if isinstance(e, list): e = e[0]
                eid = e.get('gene') if isinstance(e, dict) else None
                if eid and q not in out:
                    out[q] = eid.split('.')[0]
            break
        except Exception as ex:
            time.sleep(2)

lines = ['axis\tsymbol\tensg']
for ax, ss in AXES.items():
    for s in ss:
        lines.append('%s\t%s\t%s' % (ax, s, out.get(s, 'NA')))
p = os.path.join(BASE, 'results/tables/BMATID_marker_ensg.tsv')
open(p, 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('saved', p)
for ax, ss in AXES.items():
    print(ax, '->', ', '.join('%s=%s' % (s, out.get(s, 'NA')) for s in ss))
