"""POST HOC (written after L-C4 was seen): calibrate L-C4 against random gene sets.

Can only limit the C4 claim. 1,000 random 100-gene sets drawn from genes measured in all three
untouched cohorts, excluding HALLMARK_HYPOXIA genes, FLT1 and HIF2dep genes, and matched to the
HIF2dep mean-expression decile distribution in GSE35574. Each set is scored, the same partial Spearman
(| hyp_all) within PE + FGR-only is computed in each cohort, and the results are pooled by the same
Fisher-z rule. The output is HIF2dep's percentile among the random sets.
Also reports overlap of HIF2dep with the F1 module and with Ackerman VT Blue / VT Black.

Usage: python l4_calibration_posthoc.py DIR
"""

import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
src = (HERE / "l3_confirm.py").read_text()
ns = {"__file__": str(HERE / "l3_confirm.py"), "__name__": "l3"}
exec(compile(src[: src.index("def scores(")], "l3_prefix", "exec"), ns)
DIR = Path(sys.argv[1])
read_series, read_gpl, gene_matrix = ns["read_series"], ns["read_gpl"], ns["gene_matrix"]
HYP, HIF2 = ns["HYP"], ns["HIF2"]

cohort_z, groups = {}, {}
m, e, _ = read_series(DIR / "GSE35574_series_matrix.txt.gz")
base = m["!Sample_title"].str.replace(r"\.rep$", "", regex=True)
g = gene_matrix(e, read_gpl(DIR / "GPL6102.txt", "Symbol")).T.groupby(base.values).mean().T
mm = m.assign(base=base.values).drop_duplicates("base").set_index("base").loc[g.columns]
cohort_z["GSE35574"] = g
groups["GSE35574"] = mm["classification"].map({"PE": "PE", "IUGR": "FGR", "CONTROL": "control"})
m, e, _ = read_series(DIR / "GSE24129_series_matrix.txt.gz")
cohort_z["GSE24129"] = gene_matrix(e, read_gpl(DIR / "GPL6244.txt", "gene_assignment"))
groups["GSE24129"] = m["disease state"].map(lambda v: "PE" if v.startswith("Pre-eclampsia") else ("FGR" if v.startswith("Fetal growth") else "control"))
m, e, _ = read_series(DIR / "GSE147776_series_matrix.txt.gz")
cohort_z["GSE147776"] = gene_matrix(e, read_gpl(DIR / "GPL20844.txt", "GENE_SYMBOL"))
groups["GSE147776"] = m["disease state"].map({"Preeclampsia with severe features": "PE",
                                              "Intrauterine growth restriction associated with severe preeclampsia": "PE",
                                              "Intrauterine growth restriction": "FGR", "Control group": "control"})

Z, mean_expr = {}, {}
for c, genes in cohort_z.items():
    mean_expr[c] = genes.mean(axis=1)
    Z[c] = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T
common = set.intersection(*[set(z.columns) for z in Z.values()])
pool = sorted(common - HYP - HIF2 - {"FLT1"})
hif2_common = sorted(HIF2 & common)
deciles = pd.qcut(mean_expr["GSE35574"].reindex(sorted(common)), 10, labels=False)
need = deciles.reindex(hif2_common).value_counts()
pool_by_decile = {d: np.array([g for g in pool if deciles.get(g) == d]) for d in range(10)}


def pooled_partial(gene_list):
    zs, ws = [], []
    for c in Z:
        sel = groups[c].isin(["PE", "FGR"]).values
        z = Z[c][sel]
        score = z[gene_list].mean(axis=1).values
        hyp = z[sorted(g for g in HYP if g in z.columns)].mean(axis=1).values
        X = np.column_stack([np.ones(len(hyp)), hyp])
        resid = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        rho = stats.spearmanr(resid(score), resid(z["FLT1"].values)).statistic
        zs.append(np.arctanh(rho)); ws.append(sel.sum() - 3)
    zs, ws = np.array(zs), np.array(ws)
    return float(np.tanh((ws * zs).sum() / ws.sum()))


observed = pooled_partial(hif2_common)
rng = np.random.default_rng(0)
null = []
for _ in range(1000):
    pick = np.concatenate([rng.choice(pool_by_decile[d], n, replace=False) for d, n in need.items()])
    null.append(pooled_partial(list(pick)))
null = np.array(null)

f1 = set(pd.read_csv(ROOT / "notebooks/results/F1_module_genes.csv")["gene"])
overlap = {"F1_module": sorted(HIF2 & f1)}
s4 = pd.read_excel(ROOT / "data_external/NIHMS1883393-supplement-Table_S4.xlsx", header=None)
hdr = s4.index[s4.apply(lambda r: r.astype(str).str.upper().eq("MODULE").any() and r.astype(str).str.upper().eq("SYMBOL").any(), axis=1)][0]
s4.columns = s4.loc[hdr].astype(str).str.upper()
s4 = s4.loc[hdr + 1:]
for module in ["BLUE", "BLACK"]:
    genes_in = set(s4.loc[s4["MODULE"].astype(str).str.upper().str.contains(module), "SYMBOL"].astype(str))
    overlap[f"Ackerman_VT_{module.lower()}_of_{len(genes_in)}"] = sorted(HIF2 & genes_in)
result = {"posthoc": True, "hif2dep_genes_common": len(hif2_common), "random_pool": len(pool),
          "observed_pooled_partial_rho": observed, "null_median": float(np.median(null)),
          "null_95th": float(np.percentile(null, 95)), "percentile": float((null < observed).mean() * 100),
          "empirical_p": float(((null >= observed).sum() + 1) / 1001),
          "overlap_counts": {k: (len(v) if isinstance(v, list) else v) for k, v in overlap.items()},
          "overlap_genes": overlap}
(HERE / "L4_calibration_posthoc.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps({k: v for k, v in result.items() if k != "overlap_genes"}, indent=2))
