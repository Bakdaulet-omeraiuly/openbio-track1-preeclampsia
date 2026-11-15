"""Addendum K, step 4: one-time Family A + B1 tests in GSE75010 (rules in docs/analysis-plan.md, Addendum K).

Usage: python k_confirm_gse75010.py SERIES_MATRIX.txt GPL6244.txt
Inputs:
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE75nnn/GSE75010/matrix/GSE75010_series_matrix.txt.gz (gunzip)
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6244&targ=self&form=text&view=data
"""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
series_path, gpl_path = Path(sys.argv[1]), Path(sys.argv[2])

# ---- metadata and expression table from the series matrix
meta_rows, table, in_table = {}, [], False
for line in series_path.open():
    if line.startswith("!series_matrix_table_begin"):
        in_table = True
        continue
    if line.startswith("!series_matrix_table_end"):
        break
    if in_table:
        table.append(line.rstrip("\n").split("\t"))
    elif line.startswith("!Sample_geo_accession"):
        meta_rows["gsm"] = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
    elif line.startswith("!Sample_characteristics_ch1"):
        values = [v.strip('"') for v in line.rstrip("\n").split("\t")[1:]]
        meta_rows[values[0].split(":")[0]] = [v.split(": ", 1)[-1] for v in values]
meta = pd.DataFrame(meta_rows).set_index("gsm")
header = [h.strip('"') for h in table[0]]
expr = pd.DataFrame([r[1:] for r in table[1:]], index=[r[0].strip('"') for r in table[1:]],
                    columns=header[1:]).apply(pd.to_numeric, errors="coerce")
assert list(expr.columns) == list(meta.index) and expr.shape[1] == 157
log_applied = bool(expr.max().max() > 100)
if log_applied:
    expr = np.log2(expr + 1)

# ---- probeset -> first gene symbol (GPL6244 gene_assignment), mean per gene
def read_gpl(path):
    lines, inside = [], False
    for line in path.open():
        if line.startswith("!platform_table_begin"):
            inside = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if inside:
            lines.append(line.rstrip("\n").split("\t"))
    return pd.DataFrame(lines[1:], columns=lines[0])


gpl = read_gpl(gpl_path)
gpl = gpl[gpl["ID"].str.match(r"^\d+$", na=False)]
symbol = gpl.set_index("ID")["gene_assignment"].map(
    lambda s: s.split(" // ")[1].strip() if isinstance(s, str) and " // " in s else None)
symbol = symbol.reindex(expr.index)
genes = expr[symbol.notna()].groupby(symbol[symbol.notna()].values).mean()
z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T  # samples x genes

gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
winner = json.loads((HERE / "K_screen_winner.json").read_text())["winner"]
SETS = {"hyp_only": HYP - GLY, "gly_only": GLY - HYP, "hyp_all": HYP,
        "B_winner": set(gmt[winner["program"]]) - HYP}
d = pd.DataFrame({k: z[sorted(g for g in v if g in z.columns)].mean(axis=1) for k, v in SETS.items()})
n_genes = {k: int(sum(g in z.columns for g in v)) for k, v in SETS.items()}
d["FLT1"] = z["FLT1"]
d["PE"] = (meta["diagnosis"] == "PE").astype(int)
d["IUGR"] = (meta["iugr diagnosis"] == "Yes").astype(int)
d["ga"] = meta["ga (week)"].astype(float) + meta["ga (day)"].astype(float) / 7
d["sex"] = meta["infant gender"]
d["labour"] = (meta["attempted vaginal delivery"] == "Yes").astype(int)
d["chorio"] = (meta["chorioamnionitis diagnosis"] == "Yes").astype(int)


def one_sided(model, term, sign=1):
    t = float(model.tvalues[term])
    p2 = float(model.pvalues[term])
    return p2 / 2 if np.sign(t) == sign else 1 - p2 / 2


def contrast(frame, y, cov="ga + C(sex) + labour", sign=1):
    model = smf.ols(f"{y} ~ PE + {cov}", data=frame).fit()
    low, high = model.conf_int().loc["PE"]
    return {"n": int(model.nobs), "n_PE": int(frame.PE.sum()), "estimate": float(model.params["PE"]),
            "CI_low": float(low), "CI_high": float(high), "p_one_sided": one_sided(model, "PE", sign)}


early = d[(d.ga < 34) & ((d.PE == 1) | (d.IUGR == 0))]
early = early[~((early.PE == 0) & (early.IUGR == 1))]
pe_vs_fgr = d[(d.PE == 1) | ((d.PE == 0) & (d.IUGR == 1))]

rows = []
rows.append({"test": "A1 hyp_only: PE<34wk vs non-PE non-IUGR <34wk", **contrast(early, "hyp_only")})
rows.append({"test": "A2 gly_only: same contrast (descriptive; prediction = not raised)", **contrast(early, "gly_only")})
pe = d[d.PE == 1]
rho, p = stats.spearmanr(pe.hyp_all, pe.FLT1)
rows.append({"test": "A3 Spearman(hyp_all, FLT1) within PE", "n": len(pe), "estimate": float(rho),
             "p_one_sided": p / 2 if rho > 0 else 1 - p / 2})
rows.append({"test": "A4 hyp_only: all PE vs non-PE IUGR", **contrast(pe_vs_fgr, "hyp_only")})
rows.append({"test": f"B1 {winner['program']} minus hypoxia genes: all PE vs non-PE IUGR ({winner['direction']})",
             **contrast(pe_vs_fgr, "B_winner", sign=winner["sign"])})
primary = [0, 2, 3, 4]
holm = multipletests([rows[i]["p_one_sided"] for i in primary], method="holm")[1]
for i, q in zip(primary, holm):
    rows[i]["holm_p"] = float(q)
    rows[i]["supported_after_holm"] = bool(q < 0.05 and (i != 2 or rows[i]["estimate"] > 0))
# sensitivity (descriptive)
no_chorio = lambda f: f[f.chorio == 0]
rows.append({"test": "S A1 excluding chorioamnionitis", **contrast(no_chorio(early), "hyp_only")})
rows.append({"test": "S A1 without labour covariate", **contrast(early, "hyp_only", cov="ga + C(sex)")})
rows.append({"test": "S A4 excluding chorioamnionitis", **contrast(no_chorio(pe_vs_fgr), "hyp_only")})
rows.append({"test": "S A4 without labour covariate", **contrast(pe_vs_fgr, "hyp_only", cov="ga + C(sex)")})
rows.append({"test": "S B1 excluding chorioamnionitis", **contrast(no_chorio(pe_vs_fgr), "B_winner", sign=winner["sign"])})
rows.append({"test": "S B1 without labour covariate", **contrast(pe_vs_fgr, "B_winner", cov="ga + C(sex)", sign=winner["sign"])})
out = pd.DataFrame(rows)
out.to_csv(HERE / "K_gse75010_results.csv", index=False)
info = {"log2_applied": log_applied, "genes_after_collapse": int(genes.shape[0]), "set_genes_measured": n_genes,
        "FLT1_present": "FLT1" in z.columns, "n_early_PE": int(early.PE.sum()), "n_early_control": int((early.PE == 0).sum()),
        "n_PE": int(d.PE.sum()), "n_nonPE_IUGR": int(((d.PE == 0) & (d.IUGR == 1)).sum())}
(HERE / "K_gse75010_info.json").write_text(json.dumps(info, indent=2) + "\n")
pd.set_option("display.width", 250)
print(json.dumps(info, indent=2))
print(out.round(4).to_string(index=False))
