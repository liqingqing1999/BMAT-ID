# -*- coding: utf-8 -*-
"""Residual site genes -> symbols -> Enrichr enrichment
   (GO / KEGG / Reactome / Hallmark / WikiPathway)."""
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

import os, json, time, uuid, urllib.request, urllib.parse, ssl

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
log = open(os.path.join(BASE, 'logs/enrich_out.txt'), 'w', encoding='utf-8')
def w(s=''):
    print(s); log.write(str(s) + '\n'); log.flush()

LIBS = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023',
        'KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020', 'WikiPathway_2023_Human']

def read_tsv(path):
    rows = [l.rstrip('\n').split('\t') for l in open(path, encoding='utf-8') if l.strip()]
    hdr = rows[0]
    return [dict(zip(hdr, r)) for r in rows[1:]]

def mygene_symbols(ensgs):
    sym = {}
    ensgs = [e for e in ensgs if e]
    for i in range(0, len(ensgs), 400):
        chunk = ensgs[i:i+400]
        data = urllib.parse.urlencode({'q': ','.join(chunk), 'scopes': 'ensembl.gene',
                                       'fields': 'symbol', 'species': 'human'}).encode()
        for attempt in range(4):
            try:
                req = urllib.request.Request('https://mygene.info/v3/query', data=data,
                                             headers={'Content-Type': 'application/x-www-form-urlencoded'})
                res = json.loads(urllib.request.urlopen(req, context=ctx, timeout=90).read())
                for h in (res if isinstance(res, list) else res.get('hits', [])):
                    q = h.get('query')
                    if q and not h.get('notfound') and h.get('symbol'):
                        sym[q] = h['symbol']
                break
            except Exception as e:
                w('   mygene retry %d: %s' % (attempt + 1, e)); time.sleep(2)
    return sym

UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0'}

def multipart(fields):
    """Enrichr addList 只接受 multipart/form-data；urlencoded 一律 400。"""
    b = uuid.uuid4().hex
    body = b''
    for k, v in fields:
        body += ('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n' % (b, k, v)).encode('utf-8')
    body += ('--%s--\r\n' % b).encode()
    h = {'Content-Type': 'multipart/form-data; boundary=%s' % b}; h.update(UA)
    return body, h

def enrichr(genes, desc):
    payload, hdr = multipart([('list', '\n'.join(genes)), ('description', desc)])
    uid = None
    for attempt in range(4):
        try:
            r = urllib.request.urlopen(urllib.request.Request('https://maayanlab.cloud/Enrichr/addList',
                                                              data=payload, headers=hdr),
                                       context=ctx, timeout=90)
            uid = json.loads(r.read())['userListId']
            break
        except Exception as e:
            w('   addList retry %d: %s' % (attempt + 1, e)); time.sleep(3)
    if uid is None:
        return None
    out = {}
    for lib in LIBS:
        for attempt in range(3):
            try:
                u = ('https://maayanlab.cloud/Enrichr/enrich?userListId=%s&backgroundType=%s'
                     % (uid, urllib.parse.quote(lib)))
                res = json.loads(urllib.request.urlopen(u, context=ctx, timeout=120).read())
                out[lib] = res.get(lib, [])
                break
            except Exception as e:
                w('   %s retry %d: %s' % (lib, attempt + 1, e)); time.sleep(3)
        w('   %-34s terms=%d' % (lib, len(out.get(lib, []))))
    return out

def show(res, lib, n=15):
    t = res.get(lib) or []
    if not t:
        return
    w('')
    w('--- %s (top %d by adj.P) ---' % (lib, min(n, len(t))))
    t = sorted(t, key=lambda x: x[6])
    for row in t[:n]:
        w('  adjP=%.2e  n=%2d/%-4d  %s' % (row[6], len(row[5]), row[3], row[1]))
        if lib.startswith('GO') or 'KEGG' in lib:
            w('        genes: %s' % ', '.join(row[5][:14]))

def save_csv(res, tag):
    import csv
    p = os.path.join(BASE, 'results/tables/%s' % tag)
    with open(p, 'w', newline='', encoding='utf-8') as fh:
        cw = csv.writer(fh)
        cw.writerow(['library', 'term', 'overlap', 'pval', 'zscore', 'combined_score', 'adj_pval', 'genes'])
        for lib, rows in res.items():
            for row in sorted(rows, key=lambda x: x[6]):
                cw.writerow([lib, row[1], len(row[5]), row[2], row[3], row[4], row[6], ';'.join(row[5])])
    w('saved: %s' % p)

# ---------- 主集：494 ----------
rec = read_tsv(os.path.join(BASE, 'results/tables/BMATID_residual494_schemeA.tsv'))
w('residual schemeA genes: %d' % len(rec))
ens = [r['ensg'].split('.')[0] for r in rec]
sym = mygene_symbols(ens)
w('mygene mapped: %d / %d' % (len(sym), len(ens)))
genes = sorted({sym[e] for e in ens if e in sym})
w('unique symbols: %d' % len(genes))
w('symbols: %s' % ', '.join(genes[:60]) + (' ...' if len(genes) > 60 else ''))

json.dump(genes, open(os.path.join(BASE, 'results/tables/BMATID_residual494_symbols.json'), 'w'), indent=0)

w('')
w('===== Enrichr: 494 残差部位基因 =====')
resA = enrichr(genes, 'BMAT-ID: residual site genes after CD45 decontamination (n=494)')
if resA:
    for lib in LIBS:
        show(resA, lib, 15)
    save_csv(resA, 'BMATID_enrichment_residual494.csv')
    json.dump(resA, open(os.path.join(BASE, 'logs/enrich494_raw.json'), 'w'))

# ---------- 敏感性：1248 ----------
rec2 = read_tsv(os.path.join(BASE, 'results/tables/BMATID_residual1248_schemeB.tsv'))
ens2 = [r['ensg'].split('.')[0] for r in rec2]
sym2 = mygene_symbols([e for e in ens2 if e not in sym])
sym.update(sym2)
genes2 = sorted({sym[e] for e in ens2 if e in sym})
w('')
w('===== Enrichr: 1248 (敏感性口径, n=%d) =====' % len(genes2))
resB = enrichr(genes2, 'BMAT-ID: residual site genes, permissive filter (n=1248)')
if resB:
    for lib in ['GO_Biological_Process_2023', 'KEGG_2021_Human', 'MSigDB_Hallmark_2020']:
        show(resB, lib, 10)
    # NB: the 1248-gene enrichment table is written by 42_enrichr_sensitivity_genes.py;
    # this script only records its own run, so the two cannot overwrite each other.
    json.dump(resB, open(os.path.join(BASE, 'logs/enrich1248_raw.json'), 'w'))

# ---------- 对照：体内 2712 全集中最强的一组（污染未剔除）----------
log.close()
print('DONE')
