## 组成校正后幸存的强部位基因（M1>50% 且 M3>50%）是什么
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

import os, csv, json, urllib.request, urllib.parse, ssl, time

T = os.path.join(BASE, "results", "tables")
OUT = os.path.join(BASE, "logs", "varpart_survivors_out.txt")
con = open(OUT, 'w', encoding='utf-8')
def w(*a):
    s = ' '.join(str(x) for x in a); print(s); con.write(s + '\n'); con.flush()

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE
def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))
def gv(rows, key):
    o = {}
    for r in rows:
        try: o[r['ensg']] = float(r[key])
        except (TypeError, ValueError): pass
    return o

M1 = gv(load('BMATID_varpart_invivo.csv'), 'localization')
M2 = gv(load('BMATID_varpart_invitro.csv'), 'localization')
M3r = load('BMATID_varpart_invivo_comp.csv')
M3 = gv(M3r, 'localization')
H3 = gv(M3r, 'Hemato'); B3 = gv(M3r, 'Bone')

strong = [g for g in M1 if M1[g] > .5 and M3.get(g, 0) > .5]
w('=' * 92)
w('  组成校正后幸存的强部位基因：M1(部位)>50%% 且 M3(部位,加组成后)>50%%')
w('=' * 92)
w('  n = %d' % len(strong))
w('  其中体外部位方差(M2) > 20%% 的：%d (%.1f%%)' % (
    sum(1 for g in strong if M2.get(g, 0) > .2),
    100 * sum(1 for g in strong if M2.get(g, 0) > .2) / len(strong)))

def mygene(ens):
    out = {}
    for i in range(0, len(ens), 180):
        q = urllib.parse.urlencode({'q': ','.join(ens[i:i+180]), 'scopes': 'ensembl.gene',
                                    'fields': 'symbol,name,summary', 'species': 'human'}).encode()
        req = urllib.request.Request('https://mygene.info/v3/query', data=q,
                                     headers={'Content-Type': 'application/x-www-form-urlencoded'})
        try:
            j = json.loads(urllib.request.urlopen(req, context=ctx, timeout=60).read())
            for h in (j if isinstance(j, list) else j.get('hits', [])):
                if 'query' in h:
                    out[h['query'].split('.')[0]] = (h.get('symbol') or '', h.get('name') or '')
        except Exception as e:
            w('  mygene err:', e)
        time.sleep(0.3)
    return out

sym = mygene([g.split('.')[0] for g in strong])   # mygene 的 ensembl.gene 索引不含版本号
key = lambda g: g.split('.')[0]          # mygene 返回的 query 不带版本号
w('  symbol 映射成功 %d / %d' % (sum(1 for g in strong if sym.get(key(g), ('', ''))[0]), len(strong)))
w('')
rows = sorted(strong, key=lambda g: -(M1[g] - M3.get(g, 0)))
w('  %-11s %8s %8s %8s %8s %8s  %s' % ('symbol', '%M1', '%M3', '%M2', '%Hem', '%Bone', 'name'))
for g in rows:
    s, nm = sym.get(key(g), ('', ''))
    w('  %-11s %8.2f %8.2f %8.2f %8.2f %8.2f  %s' % (
        s or g[-10:], M1[g], M3.get(g, 0), M2.get(g, 0), H3.get(g, 0), B3.get(g, 0), nm[:40]))

w('')
w('  【家族计数】')
import collections
fam = collections.Counter()
for g in strong:
    s = sym.get(key(g), ('', ''))[0]
    if s.startswith('HOX'): fam['HOX cluster'] += 1
    elif s.startswith(('DLX', 'BARX', 'LHX', 'ALX', 'EMX', 'EN1', 'MSX', 'PITX', 'SIX', 'OTX')): fam['other homeobox'] += 1
    elif s.startswith('IGH') or s.startswith('IGK') or s.startswith('IGL') or s == 'JCHAIN': fam['Immunoglobulin'] += 1
    elif s.startswith('COL'): fam['Collagen'] += 1
    elif s.startswith('SLC'): fam['Solute carrier'] += 1
    elif s.startswith('ZNF') or s.startswith('ZBTB'): fam['Zinc finger'] += 1
    elif s.startswith('LINC') or s.endswith('-AS1'): fam['lncRNA'] += 1
    elif not s: fam['(no symbol)'] += 1
    else: fam['other'] += 1
for k, v in fam.most_common():
    w('    %-16s %d' % (k, v))

with open(os.path.join(T, 'BMATID_varpart_survivors.csv'), 'w', newline='', encoding='utf-8') as f:
    wr = csv.writer(f)
    wr.writerow(['ensg', 'symbol', 'name', 'pct_loc_M1', 'pct_loc_M3', 'pct_loc_M2', 'pct_Hemato', 'pct_Bone'])
    for g in rows:
        s, nm = sym.get(key(g), ('', ''))
        wr.writerow([g, s, nm, '%.4f' % M1[g], '%.4f' % M3.get(g, 0), '%.4f' % M2.get(g, 0),
                     '%.4f' % H3.get(g, 0), '%.4f' % B3.get(g, 0)])
w('')
w('  写出: BMATID_varpart_survivors.csv')
con.close()
print('->', OUT)
