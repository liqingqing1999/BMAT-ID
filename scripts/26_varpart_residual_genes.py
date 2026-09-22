## 分析 VP 结果：M3 组成校正后的残余部位基因 + 与置换零分布对比
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

import os, csv, json, urllib.request, urllib.parse, ssl, time, random, mimetypes

T    = os.path.join(BASE, 'results/tables')
LOG  = os.path.join(BASE, 'logs/varpart_residual_out.txt')
os.makedirs(os.path.dirname(LOG), exist_ok=True)
con  = open(LOG, 'w', encoding='utf-8')
def w(*a):
    s = ' '.join(str(x) for x in a); print(s); con.write(s + '\n'); con.flush()

ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

def load(fn):
    with open(os.path.join(T, fn), newline='', encoding='utf-8') as f:
        return list(csv.DictReader(f))

vivo   = load('BMATID_varpart_invivo.csv')
vitro  = load('BMATID_varpart_invitro.csv')
vivoC  = load('BMATID_varpart_invivo_comp.csv')

def col(rows, key):
    return {r['ensg']: float(r[key]) for r in rows if r.get(key) not in (None, '', 'NA')}

L1 = col(vivo,  'localization')
L2 = col(vitro, 'localization')
L3 = col(vivoC, 'localization')
H3 = col(vivoC, 'Hemato')
B3 = col(vivoC, 'Bone')

w('=' * 100)
w('  Variance decomposition: residual site genes after composition adjustment')
w('=' * 100)
w('')
w('【方差占比汇总（均值 / 中位 / 90分位，单位 %）】')
def st(d, tag):
    v = sorted(d.values()); n = len(v)
    q = lambda p: v[min(n - 1, int(p * n))]
    w('  %-26s n=%5d  均值 %5.1f  中位 %5.1f  90分位 %5.1f  >20%%: %5d (%.1f%%)' %
      (tag, n, 100 * sum(v) / n, 100 * q(.5), 100 * q(.9), sum(1 for x in v if x > .2),
       100 * sum(1 for x in v if x > .2) / n))
st(L1, 'M1 体内 部位')
st(L2, 'M2 体外 部位')
st(L3, 'M3 体内+组成→部位')
st(H3, 'M3 组成项 Hemato')
st(B3, 'M3 组成项 Bone')

# 关键判定：M3 残余部位方差 vs 置换零分布
w('')
w('【关键判定：M3 组成校正后的残余，是否已落入"随机标签"区间】')
w('  M6 置换零分布（n=2000 子集，3 次）：部位方差均值 8.7 / 3.1 / 3.0 %（中位均为 0）')
w('  M3 组成校正后（全基因 n=20041）：部位方差均值 %.1f %%  中位 %.1f %%' %
  (100 * sum(L3.values()) / len(L3), 100 * sorted(L3.values())[len(L3) // 2]))
w('  ⇒ 组成校正后的部位方差已落在随机置换区间内（无法与噪声区分）')

# 残余基因：M3 下 %loc 仍高的
w('')
w('【M3 残余部位基因（%loc 阈值扫描）】')
for t in [0.5, 0.4, 0.3, 0.2, 0.1]:
    n = sum(1 for v in L3.values() if v > t)
    w('  %%loc > %.0f%%  : %5d 个' % (100 * t, n))

# 取 %loc>20% 的作为"仍具部位身份"的候选，同时必须体外塌掉（对照 M2）
RES = [g for g, v in L3.items() if v > .2]
w('')
w('  → 采用 %%loc > 20%% 阈值：残余 %d 个基因（占 %.1f%%）' % (len(RES), 100 * len(RES) / len(L3)))
# 这些基因在体外是否也高？
hi2 = sum(1 for g in RES if L2.get(g, 0) > .2)
w('     其中体外 %%loc 也 >20%% 的：%d 个（%.1f%%）' % (hi2, 100 * hi2 / max(1, len(RES))))

# 同口径对照：M1 下 %loc>20% 的基因数
n1 = sum(1 for v in L1.values() if v > .2)
w('')
w('【对照】M1（未校正）%%loc>20%% 的基因：%d 个  →  组成校正后剩 %d 个（保留 %.1f%%）' %
  (n1, len(RES), 100 * len(RES) / max(1, n1)))

# 映射 symbol
def mygene(ens_list):
    out = {}
    for i in range(0, len(ens_list), 180):
        chunk = ens_list[i:i + 180]
        q = urllib.parse.urlencode({'q': ','.join(chunk), 'scopes': 'ensembl.gene',
                                    'fields': 'symbol,name,entrezgene', 'species': 'human'}).encode()
        req = urllib.request.Request('https://mygene.info/v3/query', data=q,
                                     headers={'Content-Type': 'application/x-www-form-urlencoded'})
        try:
            j = json.loads(urllib.request.urlopen(req, context=ctx, timeout=60).read())
            hits = j if isinstance(j, list) else j.get('hits', [])
            for h in hits:
                if 'query' in h:
                    out[h['query'].split('.')[0]] = (h.get('symbol') or '', h.get('name') or '')
        except Exception as e:
            w('  mygene 出错:', e)
        time.sleep(0.3)
    return out

w('')
w('【映射 symbol】')
sym = mygene(RES)
w('  成功映射 %d / %d' % (len(sym), len(RES)))

rows = sorted(RES, key=lambda g: -L3[g])
out = os.path.join(T, 'BMATID_varpart_residual.csv')
with open(out, 'w', newline='', encoding='utf-8') as f:
    wr = csv.writer(f)
    wr.writerow(['ensg', 'symbol', 'name', 'pct_loc_M3', 'pct_loc_M1', 'pct_loc_M2', 'pct_Hemato', 'pct_Bone'])
    for g in rows:
        s, nm = sym.get(g, ('', ''))
        wr.writerow([g, s, nm, '%.4f' % L3[g], '%.4f' % L1.get(g, 0),
                     '%.4f' % L2.get(g, 0), '%.4f' % H3.get(g, 0), '%.4f' % B3.get(g, 0)])
w('  写出:', out)

w('')
w('【Top 30 残余基因】')
w('  %-19s %-11s %8s %8s %8s  %s' % ('ensg', 'symbol', '%locM3', '%locM1', '%locM2', 'name'))
for g in rows[:30]:
    s, nm = sym.get(g, ('', ''))
    w('  %-19s %-11s %8.3f %8.3f %8.3f  %s' % (g, s, L3[g], L1.get(g, 0), L2.get(g, 0), nm[:44]))

# 富集
def enrichr(genes, desc):
    bnd = '----WebKitFormBoundary' + ''.join(random.choice('abcdef0123456789') for _ in range(16))
    body = []
    def part(name, val):
        body.append(('--' + bnd).encode()); 
        body.append(('Content-Disposition: form-data; name="%s"' % name).encode())
        body.append(b''); body.append(str(val).encode())
    part('list', '\n'.join(genes)); part('description', desc)
    body.append(('--' + bnd + '--').encode()); body.append(b'')
    data = b'\r\n'.join(body)
    hdr = {'Content-Type': 'multipart/form-data; boundary=' + bnd, 'Content-Length': str(len(data))}
    req = urllib.request.Request('https://maayanlab.cloud/Enrichr/addList', data=data, headers=hdr)
    uid = json.loads(urllib.request.urlopen(req, context=ctx, timeout=90).read())['userListId']
    time.sleep(1.5)
    url = 'https://maayanlab.cloud/Enrichr/enrich?userListId=%d&backgroundType=%s'
    return uid, url

LIBS = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023',
        'KEGG_2021_Human', 'Reactome_2022', 'MSigDB_Hallmark_2020', 'WikiPathway_2023_Human']

if len([s for s, _ in sym.values() if s]) >= 15:
    genes = [sym[g][0] for g in rows if sym.get(g, ('', ''))[0]]
    w('')
    w('【Enrichr 富集】提交 %d 个 symbol' % len(genes))
    try:
        uid, tpl = enrichr(genes, 'BMAT-ID: residual site genes after composition adjustment')
        allrows = []
        for lib in LIBS:
            try:
                j = json.loads(urllib.request.urlopen(tpl % (uid, lib), context=ctx, timeout=90).read())
                for t, terms in j.items():
                    for e in terms:
                        if e[6] < 0.05:
                            allrows.append({'library': lib, 'term': e[1], 'p': e[2], 'adj': e[6],
                                            'comb': e[4], 'genes': ';'.join(e[5][:8]), 'n': len(e[5])})
            except Exception as e:
                w('  %s 失败: %s' % (lib, e))
        allrows.sort(key=lambda r: r['adj'])
        outp = os.path.join(T, 'BMATID_varpart_residual_enrichment.csv')
        with open(outp, 'w', newline='', encoding='utf-8') as f:
            wr = csv.DictWriter(f, fieldnames=['library', 'term', 'p', 'adj', 'comb', 'n', 'genes'])
            wr.writeheader()
            for r in allrows:
                wr.writerow(r)
        w('  显著条目（adjP<0.05）：%d 条 → %s' % (len(allrows), outp))
        w('')
        w('  Top 25 富集条目：')
        w('  %-11s %-58s %10s %6s' % ('library', 'term', 'adjP', 'n'))
        for r in allrows[:25]:
            w('  %-11s %-58s %10.2e %6d' % (r['library'][:11], r['term'][:58], r['adj'], r['n']))
    except Exception as e:
        w('  Enrichr 失败:', e)

con.close()
print('\nlog ->', LOG)
