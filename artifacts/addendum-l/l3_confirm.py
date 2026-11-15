"""Addendum L step 3: one-time confirmation in untouched cohorts GSE35574, GSE24129, GSE147776.

Inputs (series matrices and GPL annotation tables from GEO):
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE35nnn/GSE35574/matrix/GSE35574_series_matrix.txt.gz
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24nnn/GSE24129/matrix/GSE24129_series_matrix.txt.gz
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147776/matrix/GSE147776_series_matrix.txt.gz
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL6102&targ=self&form=text&view=data  (also GPL6244, GPL20844)
Usage: python l3_confirm.py DIR   (DIR holds the three .txt.gz matrices and GPL6102.txt, GPL6244.txt, GPL20844.txt)
"""

import gzip
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DIR = Path(sys.argv[1])
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
HIF2 = set((HERE / "L1_HIF2dep_genes.txt").read_text().split())
SETS = {"hyp_only": HYP - GLY, "hyp_all": HYP, "hif2dep": HIF2}


def read_series(path):
    meta, table, inside = {}, [], False
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                inside = True
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            values = [v.strip('"') for v in line.rstrip("\n").split("\t")]
            if inside:
                table.append(values)
            elif values[0] == "!Sample_characteristics_ch1":
                meta[values[1].split(":")[0]] = [v.split(": ", 1)[-1] for v in values[1:]]
            elif values[0] in ("!Sample_geo_accession", "!Sample_description", "!Sample_title"):
                meta.setdefault(values[0], values[1:])
    m = pd.DataFrame(meta).set_index("!Sample_geo_accession")
    expr = pd.DataFrame([r[1:] for r in table[1:]], index=[r[0] for r in table[1:]], columns=table[0][1:])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    assert list(expr.columns) == list(m.index)
    log_applied = bool(expr.max().max() > 100)
    if log_applied:
        expr = np.log2(expr.clip(lower=0) + 1)
    return m, expr, log_applied


def read_gpl(path, column):
    rows, inside = [], False
    for line in path.open():
        if line.startswith("!platform_table_begin"):
            inside = True
            continue
        if line.startswith("!platform_table_end"):
            break
        if inside:
            rows.append(line.rstrip("\n").split("\t"))
    gpl = pd.DataFrame(rows[1:], columns=rows[0][: len(rows[1])] if len(rows[0]) > len(rows[1]) else rows[0])
    sym = gpl.set_index("ID")[column]
    if column == "gene_assignment":
        sym = sym.map(lambda s: s.split(" // ")[1].strip() if isinstance(s, str) and " // " in s else None)
    sym = sym.where(sym.notna() & (sym.astype(str).str.strip() != "") & (sym.astype(str) != "---"))
    return sym


def gene_matrix(expr, sym):
    sym = sym.reindex(expr.index)
    genes = expr[sym.notna()].groupby(sym[sym.notna()].values).mean()
    return genes


def scores(genes):
    z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T
    s = pd.DataFrame({k: z[sorted(g for g in v if g in z.columns)].mean(axis=1) for k, v in SETS.items()})
    s = (s - s.mean()) / s.std()
    s["FLT1"] = z["FLT1"]
    return s, {k: int(sum(g in z.columns for g in v)) for k, v in SETS.items()}


cohorts = {}
# GSE35574
m, e, lg = read_series(DIR / "GSE35574_series_matrix.txt.gz")
# Technical repeats are marked by a ".rep" title suffix (8 arrays); see the Addendum L amendment.
base = m["!Sample_title"].str.replace(r"\.rep$", "", regex=True)
assert base.nunique() == 86 and (m["!Sample_title"].str.endswith(".rep")).sum() == 8
genes = gene_matrix(e, read_gpl(DIR / "GPL6102.txt", "Symbol"))
genes = genes.T.groupby(base.values).mean().T
mm = m.assign(base=base.values).drop_duplicates("base").set_index("base").loc[genes.columns]
s, ng = scores(genes)
d = s.join(pd.DataFrame({"group": mm["classification"].map({"PE": "PE", "IUGR": "FGR", "CONTROL": "control"}),
                         "GA": mm["gestational age"].astype(float), "sex": mm["gender"],
                         "labour": mm["induction of labor"], "batch": mm["batch"]}))
cohorts["GSE35574"] = {"data": d, "cov": " + GA + C(sex) + C(labour) + C(batch)", "log2": lg, "genes": ng,
                       "arrays": int(e.shape[1]), "placentas": int(len(d))}
# GSE24129
m, e, lg = read_series(DIR / "GSE24129_series_matrix.txt.gz")
s, ng = scores(gene_matrix(e, read_gpl(DIR / "GPL6244.txt", "gene_assignment")))
grp = m["disease state"].map(lambda v: "PE" if v.startswith("Pre-eclampsia") else ("FGR" if v.startswith("Fetal growth") else "control"))
cohorts["GSE24129"] = {"data": s.join(grp.rename("group")), "cov": "", "log2": lg, "genes": ng,
                       "arrays": int(e.shape[1]), "placentas": int(e.shape[1])}
# GSE147776
m, e, lg = read_series(DIR / "GSE147776_series_matrix.txt.gz")
s, ng = scores(gene_matrix(e, read_gpl(DIR / "GPL20844.txt", "GENE_SYMBOL")))
grp = m["disease state"].map({"Preeclampsia with severe features": "PE",
                              "Intrauterine growth restriction associated with severe preeclampsia": "PE",
                              "Intrauterine growth restriction": "FGR", "Control group": "control"})
cohorts["GSE147776"] = {"data": s.join(grp.rename("group")), "cov": "", "log2": lg, "genes": ng,
                        "arrays": int(e.shape[1]), "placentas": int(e.shape[1])}


def coef(frame, y, other, cov):
    f = frame[frame.group.isin(["PE", other])].copy()
    f["PE"] = (f.group == "PE").astype(int)
    model = smf.ols(f"{y} ~ PE{cov}", data=f).fit()
    return float(model.params["PE"]), float(model.bse["PE"]), int(f.PE.sum()), int((f.PE == 0).sum())


def partial(frame):
    X = np.column_stack([np.ones(len(frame)), frame["hyp_all"].values])
    resid = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return float(stats.spearmanr(resid(frame["hif2dep"].values), resid(frame["FLT1"].values)).statistic)


def pool_coef(items):
    est = np.array([i[0] for i in items]); se = np.array([i[1] for i in items]); w = 1 / se**2
    fe = float((w * est).sum() / w.sum()); fe_se = float(1 / np.sqrt(w.sum()))
    q = float((w * (est - fe) ** 2).sum()); k = len(items)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum()))
    wr = 1 / (se**2 + tau2); re_ = float((wr * est).sum() / wr.sum()); re_se = float(1 / np.sqrt(wr.sum()))
    return {"fixed": fe, "fixed_se": fe_se, "CI_low": fe - 1.96 * fe_se, "CI_high": fe + 1.96 * fe_se,
            "z": fe / fe_se, "random": re_, "random_se": re_se, "tau2": tau2}


def pool_rho(items):
    rho = np.array([i[0] for i in items]); n = np.array([i[1] for i in items]); w = n - 3
    zz = float((w * np.arctanh(rho)).sum() / w.sum()); se = float(1 / np.sqrt(w.sum()))
    return {"pooled_rho": float(np.tanh(zz)), "CI_low": float(np.tanh(zz - 1.96 * se)),
            "CI_high": float(np.tanh(zz + 1.96 * se)), "z": zz / se}


per, results = [], {}
names = list(cohorts)


def run(selection, cov_override=None):
    out = {}
    c1 = [coef(cohorts[c]["data"], "hyp_only", "FGR", cov_override if (cov_override is not None and c == "GSE35574") else cohorts[c]["cov"]) for c in selection]
    c2 = [coef(cohorts[c]["data"], "hyp_only", "control", cov_override if (cov_override is not None and c == "GSE35574") else cohorts[c]["cov"]) for c in selection]
    c3 = []
    c4 = []
    for c in selection:
        d = cohorts[c]["data"]
        pe = d[d.group == "PE"]
        c3.append((float(stats.spearmanr(pe.hyp_all, pe.FLT1).statistic), len(pe)))
        dis = d[d.group.isin(["PE", "FGR"])]
        c4.append((partial(dis), len(dis)))
    out["C1"] = pool_coef(c1); out["C1"]["p_two_sided"] = float(2 * (1 - stats.norm.cdf(abs(out["C1"]["z"]))))
    out["C2"] = pool_coef(c2); out["C2"]["p_one_sided"] = float(1 - stats.norm.cdf(out["C2"]["z"]))
    out["C3"] = pool_rho(c3); out["C3"]["p_one_sided"] = float(1 - stats.norm.cdf(out["C3"]["z"]))
    out["C4"] = pool_rho(c4); out["C4"]["p_one_sided"] = float(1 - stats.norm.cdf(out["C4"]["z"]))
    out["per_cohort"] = {c: {"C1_PE_vs_FGR": c1[i], "C2_PE_vs_control": c2[i], "C3_rho_n": c3[i], "C4_partial_rho_n": c4[i]}
                         for i, c in enumerate(selection)}
    return out


main = run(names)
holm = multipletests([main["C1"]["p_two_sided"], main["C2"]["p_one_sided"], main["C3"]["p_one_sided"],
                      main["C4"]["p_one_sided"]], method="holm")[1]
for key, q in zip(["C1", "C2", "C3", "C4"], holm):
    main[key]["holm_p"] = float(q)
result = {
    "cohorts": {c: {k: v for k, v in cohorts[c].items() if k != "data"} | {
        "group_counts": cohorts[c]["data"].group.value_counts().to_dict()} for c in names},
    "main": main,
    "sensitivity_GSE35574_without_labour_batch": run(names, cov_override=" + GA + C(sex)"),
    "leave_one_out": {f"without_{c}": {k: v for k, v in run([x for x in names if x != c]).items() if k != "per_cohort"} for c in names},
}
(HERE / "L3_confirmation.json").write_text(json.dumps(result, indent=2, default=float) + "\n")
print(json.dumps({"cohorts": result["cohorts"], "main": main}, indent=2, default=float))
