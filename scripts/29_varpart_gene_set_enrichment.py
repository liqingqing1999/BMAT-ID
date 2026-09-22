## 对 VP 的两个基因集做 Enrichr 富集（multipart 提交）
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

import os, csv, json, urllib.request, urllib.parse, ssl, time, random

T   = os.path.join(BASE, "results", "tables")
OUT = os.path.join(BASE, "logs", "varpart_enrich_out.txt")
con = open(OUT, 'w', encoding='utf-8')
def w(*a):
    s = ' '.join(str(x) for x in a); print(s); con.write(s + '\n'); con.flush()
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

LIBS = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023',
        'KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020', 'WikiPathway_2023_Human',
        'Jensen_TISSUES', 'GTEx_Tissue_Expression_Up']

def addlist(genes, desc):
    bnd = '----WebKitFormBoundary' + ''.join(random.choice('abcdef0123456789') for _ in range(16))
    body = []
    for name, val in [('list', '\n'.join(genes)), ('description', desc)]:
        body.append(('--' + bnd).encode())
        body.append(('Content-Disposition: form-data; name="%s"' % name).encode())
        body.append(b''); body.append(str(val).encode())
    body.append(('--' + bnd + '--').encode()); body.append(b'')
    data = b'\r\n'.join(body)
    req = urllib.request.Request('https://maayanlab.cloud/Enrichr/addList', data=data,
                                 headers={'Content-Type': 'multipart/form-data; boundary=' + bnd,
                                          'Content-Length': str(len(data))})
    return json.loads(urllib.request.urlopen(req, context=ctx, timeout=90).read())['userListId']

def run(genes, tag):
    genes = [g for g in genes if g]
    w('')
    w('=' * 96)
    w('  %s   n=%d' % (tag, len(genes)))
    w('=' * 96)
    if len(genes) < 15:
        w('  基因太少，跳过'); return
    try:
        uid = addlist(genes, tag)
    except Exception as e:
        w('  addList 失败:', e); return
    time.sleep(1.5)
    rows = []
    for lib in LIBS:
        try:
            u = 'https://maayanlab.cloud/Enrichr/enrich?userListId=%d&backgroundType=%s' % (uid, lib)
            j = json.loads(urllib.request.urlopen(u, context=ctx, timeout=90).read())
            for _, terms in j.items():
                for e in terms:
                    if e[6] < 0.05:
                        rows.append({'gene_set': lib, 'term': e[1], 'p': e[2], 'adj': e[6],
                                     'combined': e[4], 'n_genes': len(e[5]), 'genes': ';'.join(e[5][:10])})
        except Exception as e:
            w('  %s 失败: %s' % (lib, e))
    rows.sort(key=lambda r: r['adj'])
    outp = os.path.join(T, 'BMATID_varpart_enrichment_%s.csv' % tag.split()[0])
    with open(outp, 'w', newline='', encoding='utf-8') as f:
        wr = csv.DictWriter(f, fieldnames=['gene_set', 'term', 'p', 'adj', 'combined', 'n_genes', 'genes'])
        wr.writeheader()
        for r in rows:
            wr.writerow(r)
    w('  显著条目（adjP<0.05）：%d → %s' % (len(rows), os.path.basename(outp)))
    w('')
    w('  %-26s %-52s %10s %5s' % ('gene_set', 'term', 'adjP', 'n'))
    for r in rows[:22]:
        w('  %-26s %-52s %10.2e %5d' % (r['gene_set'][:26], r['term'][:52], r['adj'], r['n_genes']))
    return rows

# 基因集 1：233 个幸存者（M1>50% 且 M3>50%）
s = list(csv.DictReader(open(os.path.join(T, 'BMATID_varpart_survivors.csv'), encoding='utf-8')))
run([r['symbol'] for r in s], 'survivors233 组成校正后幸存强部位基因')

# 基因集 2：1591 个（M3 %loc>20%）
r2 = list(csv.DictReader(open(os.path.join(T, 'BMATID_varpart_residual.csv'), encoding='utf-8')))
run([r['symbol'] for r in r2], 'residual1591 组成校正后 %loc>20%')

con.close()
print('->', OUT)
