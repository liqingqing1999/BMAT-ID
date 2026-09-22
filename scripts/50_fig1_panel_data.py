## 1c prep: map all genes to symbols, merge the four variance-decomposition tables,
## and report what actually happens to three candidate gene groups.
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

import io, os, json, time
import numpy as np, pandas as pd

T = os.path.join(BASE, "results", "tables") + "/"
L = os.path.join(BASE, "logs", "figc_prep_out.txt")
os.makedirs(os.path.dirname(L), exist_ok=True)
log = io.open(L, "w", encoding="utf-8", buffering=1)
def say(*a):
    log.write(" ".join(str(x) for x in a) + "\n")

say("== load ==")
viv = pd.read_csv(T + "BMATID_varpart_invivo.csv")
vco = pd.read_csv(T + "BMATID_varpart_invivo_comp.csv")
vit = pd.read_csv(T + "BMATID_varpart_invitro.csv")
c45 = pd.read_csv(T + "BMATID_primary_cd45_annotated.csv")
sur = pd.read_csv(T + "BMATID_varpart_survivors.csv")
say("invivo", viv.shape, list(viv.columns))
say("invivo_comp", vco.shape, list(vco.columns))
say("invitro", vit.shape, list(vit.columns))
say("cd45", c45.shape, list(c45.columns))
say("survivors", sur.shape, list(sur.columns))

def base(s):
    return s.astype(str).str.split(".").str[0]

df = pd.DataFrame({"ensg": viv["ensg"].astype(str), "ensg_b": base(viv["ensg"])})
df["M1"] = viv["localization"].astype(float)
df = df.merge(pd.DataFrame({"ensg_b": base(vco["ensg"]), "M3": vco["localization"].astype(float),
                            "Hemato": vco["Hemato"].astype(float), "Bone": vco["Bone"].astype(float)}),
              on="ensg_b", how="left")
df = df.merge(pd.DataFrame({"ensg_b": base(vit["ensg"]), "M2": vit["localization"].astype(float)}),
              on="ensg_b", how="left")
df = df.merge(pd.DataFrame({"ensg_b": base(c45["ensg"]), "r_cd45": c45["r_cd45"].astype(float),
                            "F_prim": c45["F_prim"].astype(float), "FDR_prim": c45["FDR_prim"].astype(float)}),
              on="ensg_b", how="left")
say("merged", df.shape, "na:", df[["M1","M3","M2","r_cd45"]].isna().sum().to_dict())

## ---- symbol mapping (mygene, cached)
cache = T + "BMATID_allgenes_symbols.csv"
if os.path.exists(cache):
    sym = pd.read_csv(cache)
    say("symbol cache loaded", sym.shape)
else:
    import urllib.request, urllib.parse, ssl as _ssl
    ctx = _ssl.create_default_context(); ctx.check_hostname = False; ctx.verify_mode = _ssl.CERT_NONE
    ids = sorted(df["ensg_b"].unique().tolist())
    say("querying mygene for", len(ids), "ids")
    hits = {}
    for i in range(0, len(ids), 180):
        chunk = ids[i:i + 180]
        q = urllib.parse.urlencode({"q": ",".join(chunk), "scopes": "ensembl.gene",
                                    "fields": "symbol,name", "species": "human"}).encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=q,
                                     headers={"Content-Type": "application/x-www-form-urlencoded"})
        try:
            j = json.loads(urllib.request.urlopen(req, context=ctx, timeout=60).read())
            for h in (j if isinstance(j, list) else j.get("hits", [])):
                if "query" in h:
                    hits[h["query"].split(".")[0]] = (h.get("symbol"), h.get("name"))
        except Exception as e:
            say("  mygene err:", repr(e))
        time.sleep(0.25)
    rows = [(k, v[0], v[1]) for k, v in hits.items()]
    sym = pd.DataFrame(rows, columns=["ensg_b", "symbol", "name"])
    sym.to_csv(cache, index=False, encoding="utf-8")
    say("mapped:", sym["symbol"].notna().sum(), "/", len(sym))
df = df.merge(sym, on="ensg_b", how="left")
say("with symbol:", df["symbol"].notna().sum(), "/", len(df))
df.to_csv(T + "BMATID_varpart_allgenes_merged.csv", index=False, encoding="utf-8")
say("wrote BMATID_varpart_allgenes_merged.csv")

## ---- what happens to candidate groups
groups = {
 "axis_TF (BARX1,DLX5,TBX5,TBX15,TBX3,TBX18,MEOX2,SHOX2,HOXA6,HOXB3,HOXB4,HOXB7,EN1,PYGO1,RSPO1)":
   ["BARX1","DLX5","TBX5","TBX15","TBX3","TBX18","MEOX2","SHOX2","HOXA6","HOXB3","HOXB4","HOXB7","EN1","PYGO1","RSPO1"],
 "Ig/plasma (IGHG2,IGHM,IGHD,IGHA2,IGKC,IGKV1-5,IGLV3-25,IGLV2-14,IGLV4-69,MZB1,JCHAIN,SDC1,SLAMF7)":
   ["IGHG2","IGHM","IGHD","IGHA2","IGKC","IGKV1-5","IGLV3-25","IGLV2-14","IGLV4-69","MZB1","JCHAIN","SDC1","SLAMF7"],
 "adipo core (PPARG,CEBPA,PLIN1,PLIN4,FABP4,ADIPOQ,LPL,PLIN2,CIDEC,CFD)":
   ["PPARG","CEBPA","PLIN1","PLIN4","FABP4","ADIPOQ","LPL","PLIN2","CIDEC","CFD"],
 "bone/mineral (RUNX2,SP7,ALPL,IBSP,SPP1,COL1A1,BGLAP,DMP1,PHEX,MEPE,ENPP1)":
   ["RUNX2","SP7","ALPL","IBSP","SPP1","COL1A1","BGLAP","DMP1","PHEX","MEPE","ENPP1"],
 "haemato classic (PTPRC,CD14,CD68,LYZ,CD3E,MS4A1,CD19,CD79A,CSF1R)":
   ["PTPRC","CD14","CD68","LYZ","CD3E","MS4A1","CD19","CD79A","CSF1R"],
}
say("\n===== group behaviour: M1 (in vivo site) -> M3 (after composition) -> M2 (in vitro) =====")
for gname, genes in groups.items():
    sub = df[df["symbol"].isin(genes)]
    got = sorted(sub["symbol"].dropna().unique())
    miss = [g for g in genes if g not in got]
    say("\n-- ", gname)
    say("   found %d/%d  missing: %s" % (len(got), len(genes), ",".join(miss) or "-"))
    if len(sub):
        say("   median M1=%.3f  M3=%.3f  M2=%.3f  r_cd45=%.3f  |M3>0.5: %d/%d  |M3<0.1: %d/%d  |M2>0.2: %d/%d" % (
            sub["M1"].median(), sub["M3"].median(), sub["M2"].median(), sub["r_cd45"].median(),
            (sub["M3"] > 0.5).sum(), len(sub), (sub["M3"] < 0.1).sum(), len(sub),
            (sub["M2"] > 0.2).sum(), len(sub)))
        show = sub[["symbol","M1","M3","M2","Hemato","Bone","r_cd45"]].sort_values("M1", ascending=False)
        for _, r in show.iterrows():
            say("     %-12s M1=%.2f M3=%.2f M2=%.2f Hem=%.2f Bone=%.2f rCD45=%+.2f" % (
                r["symbol"], r["M1"], r["M3"], r["M2"], r["Hemato"], r["Bone"], r["r_cd45"]))

## ---- global relation: r_cd45 -> survival ratio
d2 = df.dropna(subset=["M1","M3","r_cd45"])
d2 = d2[d2["M1"] > 0.5]
d2["ratio"] = d2["M3"] / d2["M1"]
say("\n===== among genes with M1>0.5 (n=%d), M3/M1 ratio by r_cd45 decile =====" % len(d2))
d2["bin"] = pd.cut(d2["r_cd45"], bins=[-1,-0.7,-0.5,-0.3,0,0.3,0.5,0.7,1])
say(d2.groupby("bin", observed=True)["ratio"].agg(["count","median"]).to_string())
say("\nmedian ratio for |r|>0.7: %.3f ; for |r|<0.3: %.3f" % (
    d2[d2["r_cd45"].abs() > 0.7]["ratio"].median(), d2[d2["r_cd45"].abs() < 0.3]["ratio"].median()))

say("\nDONE")
log.close()
print("ok")
