
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

import json, urllib.request, urllib.parse, ssl, sys
ctx = ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = ssl.CERT_NONE

# 来自 logs/deseq2_out.txt 第 19-34 行的体外残余显著基因（全基因 DESeq2 LRT，FDR<0.05）
rows = [
    ("ENSG00000123080",  0.63, 398.2, "9.32e-10"),
    ("ENSG00000077943", -2.10, 572.6, "2.04e-06"),
    ("ENSG00000183421", -2.29,  20.7, "7.68e-04"),
    ("ENSG00000182742",  0.28,  16.9, "1.00e-03"),
    ("ENSG00000162631",  2.09, 142.7, "3.28e-03"),
    ("ENSG00000038295", -0.12,  24.3, "3.28e-03"),
    ("ENSG00000082074",  3.14,  52.7, "3.28e-03"),
    ("ENSG00000151025",  2.67,  34.9, "3.28e-03"),
    ("ENSG00000185567",  2.54, 3694.2, "3.28e-03"),
    ("ENSG00000250802",  1.22, 181.6, "6.63e-03"),
    ("ENSG00000138678",  0.97,  16.4, "1.46e-02"),
    ("ENSG00000260027",  0.48,  19.3, "1.46e-02"),
    ("ENSG00000163545", -1.07, 111.6, "3.11e-02"),
    ("ENSG00000182732", -2.36,  66.7, "4.24e-02"),
    ("ENSG00000151790", -1.77,  95.1, "4.52e-02"),
    ("ENSG00000120093",  0.76,  46.6, "4.52e-02"),
]

out = []
out.append("=" * 96)
out.append("  体外残余部位基因（DESeq2 LRT FDR<0.05, n=16）symbol 映射  mygene.info")
out.append("=" * 96)
out.append("%-20s %-12s %-10s %9s %10s %-11s" % ("ensembl", "symbol", "type", "log2FC", "baseMean", "FDR"))
out.append("-" * 96)

info = {}
for e, lfc, bm, fdr in rows:
    url = "https://mygene.info/v3/query?" + urllib.parse.urlencode(
        {"q": e, "scopes": "ensembl.gene", "fields": "symbol,name,type_of_gene,genomic_pos,alias", "species": "human"})
    sym, name, typ, pos = "?", "", "", ""
    try:
        j = json.loads(urllib.request.urlopen(url, context=ctx, timeout=45).read())
        hits = j if isinstance(j, list) else j.get("hits", [])
        if hits:
            h = hits[0]
            sym = h.get("symbol", "?")
            name = h.get("name", "")
            typ = h.get("type_of_gene", "")
            g = h.get("genomic_pos")
            if isinstance(g, list): g = g[0]
            if isinstance(g, dict): pos = "%s:%s" % (g.get("chr"), g.get("start"))
    except Exception as ex:
        name = "ERR %s" % ex
    info[e] = dict(sym=sym, name=name, typ=typ, pos=pos, lfc=lfc, bm=bm, fdr=fdr)
    out.append("%-20s %-12s %-10s %+9.2f %10.1f %-11s" % (e, sym, typ, lfc, bm, fdr))

out.append("")
out.append("=" * 96)
out.append("  基因全名")
out.append("=" * 96)
for e, d in info.items():
    out.append("  %-12s %-20s %s" % (d["sym"], e, d["name"][:70]))
    if d["pos"]:
        out.append("  %-12s %-20s %s" % ("", "", d["pos"]))

out.append("")
out.append("=" * 96)
out.append("  按功能类别粗看")
out.append("=" * 96)
cat_order = [
    ("细胞周期/CDKN2A 家族", ["CDKN2A", "CDKN2B", "CDKN2A-AS1", "CDKN2B-AS1"]),
    ("HOX / 位置身份", ["HOXA", "HOXB", "HOXC", "HOXD"]),
    ("细胞黏附/ECM", ["ITGA", "ITGB", "COL", "LAMA", "LAMB", "FN1", "SPP1", "MMP", "ADAM"]),
    ("GPCR/信号", ["GPR", "ADGR", "RXFP", "PTH", "WNT", "FZD", "BMP"]),
    ("代谢/成脂", ["ADIPOQ", "PLIN", "FABP", "LPL", "PPARG", "CEBP", "LEP", "CIDEC", "GPAM"]),
    ("转录因子", ["SOX", "TWIST", "ZEB", "SNAI", "PRRX", "MEOX", "MSC"]),
]
used = set()
for label, keys in cat_order:
    hit = [(d["sym"], e, d["name"]) for e, d in info.items()
           if any(d["sym"].upper().startswith(k) for k in keys)]
    if hit:
        out.append("  [%s]" % label)
        for s, e, n in hit:
            out.append("     %-12s %-20s %s" % (s, e, n[:60]))
            used.add(e)
left = [(d["sym"], e, d["name"]) for e, d in info.items() if e not in used]
if left:
    out.append("  [未归类]")
    for s, e, n in left:
        out.append("     %-12s %-20s %s" % (s, e, n[:60]))

txt = "\n".join(out)
open(os.path.join(BASE, "logs", "deseq2_residual_symbols.txt"), "w", encoding="utf-8").write(txt)
print(txt)
