# -*- coding: utf-8 -*-
"""敏感性口径 1248 基因的 Enrichr 富集（只存 CSV，stdout 极简）"""
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

import os, json, time, uuid, csv, urllib.request, ssl
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
UA = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) research-tool/1.0'}
LIBS = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'KEGG_2021_Human',
        'Reactome_2022', 'MSigDB_Hallmark_2020', 'WikiPathway_2023_Human']

def multipart(fields):
    b = uuid.uuid4().hex
    body = b''
    for k, v in fields:
        body += ('--%s\r\nContent-Disposition: form-data; name="%s"\r\n\r\n%s\r\n' % (b, k, v)).encode('utf-8')
    body += ('--%s--\r\n' % b).encode()
    h = {'Content-Type': 'multipart/form-data; boundary=%s' % b}; h.update(UA)
    return body, h

rows = [l.rstrip('\n').split('\t') for l in open(os.path.join(BASE, 'results/tables/BMATID_residual1248_schemeB.tsv'), encoding='utf-8') if l.strip()]
hdr = rows[0]; ens = [r[hdr.index('ensg')].split('.')[0] for r in rows[1:]]

# symbol 映射（把 494 已映射的复用，其余现查）
sym = {}
p494 = os.path.join(BASE, 'results/tables/BMATID_residual494_symbols.json')
have = set(json.load(open(p494)))
def mygene(ids):
    out = {}
    ids = [i for i in ids if i]
    import urllib.parse
    for i in range(0, len(ids), 400):
        ch = ids[i:i+400]
        d = urllib.parse.urlencode({'q': ','.join(ch), 'scopes': 'ensembl.gene', 'fields': 'symbol',
                                    'species': 'human'}).encode()
        for a in range(4):
            try:
                r = urllib.request.Request('https://mygene.info/v3/query', data=d,
                                           headers={'Content-Type': 'application/x-www-form-urlencoded'})
                res = json.loads(urllib.request.urlopen(r, context=ctx, timeout=90).read())
                for h in (res if isinstance(res, list) else res.get('hits', [])):
                    if h.get('query') and not h.get('notfound') and h.get('symbol'):
                        out[h['query']] = h['symbol']
                break
            except Exception:
                time.sleep(2)
    return out
sym = mygene(ens)
genes = sorted({sym[e] for e in ens if e in sym})
print('mapped %d/%d -> %d symbols' % (len(sym), len(ens), len(genes)))

body, hd = multipart([('list', '\n'.join(genes)), ('description', 'BMAT-ID: residual site genes, permissive filter (n=1248)')])
uid = None
for a in range(4):
    try:
        uid = json.loads(urllib.request.urlopen(urllib.request.Request('https://maayanlab.cloud/Enrichr/addList', data=body, headers=hd), context=ctx, timeout=120).read())['userListId']
        break
    except Exception as e:
        print('retry', a, e); time.sleep(3)
assert uid
res = {}
for lib in LIBS:
    for a in range(3):
        try:
            import urllib.parse
            u = 'https://maayanlab.cloud/Enrichr/enrich?userListId=%s&backgroundType=%s' % (uid, urllib.parse.quote(lib))
            res[lib] = json.loads(urllib.request.urlopen(urllib.request.Request(u, headers=UA), context=ctx, timeout=180).read()).get(lib, [])
            break
        except Exception as e:
            time.sleep(3)
    print('  %-32s %d terms' % (lib, len(res.get(lib, []))))
p = os.path.join(BASE, 'results/tables/BMATID_enrichment_residual1248.csv')
with open(p, 'w', newline='', encoding='utf-8') as fh:
    cw = csv.writer(fh); cw.writerow(['library','term','overlap','pval','zscore','combined_score','adj_pval','genes'])
    for lib, rr in res.items():
        for r in sorted(rr, key=lambda x: x[6]):
            cw.writerow([lib, r[1], len(r[5]), r[2], r[3], r[4], r[6], ';'.join(r[5])])
print('saved', p)
