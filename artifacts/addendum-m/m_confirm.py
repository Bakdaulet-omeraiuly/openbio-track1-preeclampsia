"""Addendum M: one-time confirmation in GSE44711, GSE74341, GSE66273 (rules in plan.md).

Inputs (GEO series matrices and GPL annotation tables):
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE44nnn/GSE44711/matrix/GSE44711_series_matrix.txt.gz
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74341/matrix/GSE74341_series_matrix.txt.gz
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE66nnn/GSE66273/matrix/GSE66273_series_matrix.txt.gz
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558&targ=self&form=text&view=data  (also GPL16699, GPL4133)
Usage: python m_confirm.py DIR
"""

import gzip
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
DIR = Path(sys.argv[1])
NPERM, SEED = 10_000, 42
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
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
            elif values[0] in ("!Sample_characteristics_ch1", "!Sample_characteristics_ch2"):
                ch = values[0][-3:]
                meta[f"{ch}:{values[1].split(':')[0]}"] = [v.split(": ", 1)[-1] for v in values[1:]]
            elif values[0] in ("!Sample_geo_accession", "!Sample_title", "!Sample_label_ch1", "!Sample_label_ch2"):
                meta.setdefault(values[0], values[1:])
    m = pd.DataFrame(meta).set_index("!Sample_geo_accession")
    expr = pd.DataFrame([r[1:] for r in table[1:]], index=[r[0] for r in table[1:]], columns=table[0][1:])
    expr = expr.apply(pd.to_numeric, errors="coerce")
    assert list(expr.columns) == list(m.index)
    return m, expr


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
    gpl = pd.DataFrame(rows[1:], columns=rows[0][: len(rows[1])])
    sym = gpl.set_index("ID")[column].astype(str).str.strip()
    return sym.where((sym != "") & (sym != "nan"))


def gene_matrix(expr, sym):
    sym = sym.reindex(expr.index)
    keep = sym.notna()
    return expr[keep].groupby(sym[keep].values).mean().dropna(how="all")


def scores(genes):
    genes = genes[genes.notna().all(axis=1)]
    z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T
    s = pd.DataFrame({k: z[sorted(g for g in v if g in z.columns)].mean(axis=1) for k, v in SETS.items()})
    s = (s - s.mean()) / s.std()
    s["FLT1"] = z["FLT1"]
    return s, {k: int(sum(g in z.columns for g in v)) for k, v in SETS.items()}


def maybe_log(expr):
    applied = bool(expr.max().max() > 100)
    return (np.log2(expr.clip(lower=0) + 1) if applied else expr), applied


cohorts = {}
# GSE44711: 8 EOPET vs 8 GA-matched controls, chorionic villi, Illumina HT-12 v4
m, e = read_series(DIR / "GSE44711_series_matrix.txt.gz")
e, lg = maybe_log(e)
s, ng = scores(gene_matrix(e, read_gpl(DIR / "GPL10558.txt", "Symbol")))
grp = m["ch1:condition"].map({"EOPET": "PE", "Control": "control"})
assert grp.value_counts().to_dict() == {"PE": 8, "control": 8}, grp.value_counts()
d = s.join(pd.DataFrame({"group": grp, "GA": m["ch1:gestational age (weeks)"].astype(float)}))
cohorts["GSE44711"] = {"data": d, "cov": " + GA", "log2_applied": lg, "genes": ng, "tissue": "chorionic villi"}

# GSE74341: two-channel, VALUE = log2(Cy3/Cy5) = sample (ch2) / pooled reference (ch1)
m, e = read_series(DIR / "GSE74341_series_matrix.txt.gz")
assert (m["!Sample_label_ch2"] == "Cy3").all() and (m["!Sample_label_ch1"] == "Cy5").all()
s, ng = scores(gene_matrix(e, read_gpl(DIR / "GPL16699.txt", "GENE_SYMBOL")))
status = m["ch2:subject status"]
grp = status.map(lambda v: "PE" if v.startswith("patients with early-onset") else
                 ("control" if v.startswith("controls who delivered preterm") else
                  ("LOPE" if v.startswith("patients with late-onset") else "term_control")))
assert grp.value_counts().to_dict() == {"LOPE": 8, "PE": 7, "control": 5, "term_control": 5}, grp.value_counts()
cohorts["GSE74341"] = {"data": s.join(grp.rename("group")), "cov": "", "log2_applied": False, "genes": ng,
                       "tissue": "placental tissue"}

# GSE66273: 12 preterm PE (6 with HELLP) vs 5 preterm controls, chorionic plate
m, e = read_series(DIR / "GSE66273_series_matrix.txt.gz")
e, lg = maybe_log(e)
s, ng = scores(gene_matrix(e, read_gpl(DIR / "GPL4133.txt", "GENE_SYMBOL")))
cond = m["ch1:condition"]
grp = cond.map(lambda v: "control" if v.startswith("Preterm control") else "PE")
assert grp.value_counts().to_dict() == {"PE": 12, "control": 5}, grp.value_counts()
d = s.join(pd.DataFrame({"group": grp, "GA": m["ch1:gestational weeks"].astype(float),
                         "sex": m["ch1:baby gender"], "HELLP": cond.str.contains("HELLP")}))
cohorts["GSE66273"] = {"data": d, "cov": " + GA + C(sex)", "log2_applied": lg, "genes": ng, "tissue": "chorionic plate"}


def coef(f, cov):
    f = f[f.group.isin(["PE", "control"])].copy()
    f["PE"] = (f.group == "PE").astype(int)
    model = smf.ols(f"hyp_only ~ PE{cov}", data=f).fit()
    return float(model.params["PE"]), float(model.bse["PE"]), int(f.PE.sum()), int((f.PE == 0).sum())


def pool_coef(items):
    est = np.array([i[0] for i in items]); se = np.array([i[1] for i in items]); w = 1 / se**2
    fe = float((w * est).sum() / w.sum()); fe_se = float(1 / np.sqrt(w.sum()))
    q = float((w * (est - fe) ** 2).sum()); k = len(items)
    tau2 = max(0.0, (q - (k - 1)) / (w.sum() - (w**2).sum() / w.sum())) if k > 1 else 0.0
    return {"fixed": fe, "CI_low": fe - 1.96 * fe_se, "CI_high": fe + 1.96 * fe_se, "z": fe / fe_se,
            "p_normal_one_sided": float(1 - stats.norm.cdf(fe / fe_se)), "tau2": tau2}


def partial(frame):
    X = np.column_stack([np.ones(len(frame)), frame["hyp_all"].values])
    resid = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return float(stats.spearmanr(resid(frame["hif2dep"].values), resid(frame["FLT1"].values)).statistic)


def pool_rho(items):
    rho = np.array([i[0] for i in items]); n = np.array([i[1] for i in items]); w = n - 3
    zz = float((w * np.arctanh(rho)).sum() / w.sum()); se = float(1 / np.sqrt(w.sum()))
    return {"pooled_rho": float(np.tanh(zz)), "CI_low": float(np.tanh(zz - 1.96 * se)),
            "CI_high": float(np.tanh(zz + 1.96 * se)), "z": zz / se,
            "p_normal_one_sided": float(1 - stats.norm.cdf(zz / se))}


def m1(frames, covs):
    return pool_coef([coef(frames[c], covs[c]) for c in frames])


def m2_items(frames):
    out = []
    for f in frames.values():
        pe = f[f.group == "PE"]
        out.append((float(stats.spearmanr(pe.hyp_all, pe.FLT1).statistic), len(pe)))
    return out


def m1_perm(frames, covs, rng):
    obs = m1(frames, covs)["z"]
    null = np.empty(NPERM)
    for i in range(NPERM):
        perm = {}
        for c, f in frames.items():
            g = f.copy()
            mask = g.group.isin(["PE", "control"])
            g.loc[mask, "group"] = rng.permutation(g.loc[mask, "group"].values)
            perm[c] = g
        null[i] = m1(perm, covs)["z"]
    return float((1 + (null >= obs).sum()) / (NPERM + 1))


def m2_perm(frames, rng):
    obs = pool_rho(m2_items(frames))["z"]
    null = np.empty(NPERM)
    for i in range(NPERM):
        perm = {}
        for c, f in frames.items():
            pe = f[f.group == "PE"].copy()
            pe["FLT1"] = rng.permutation(pe["FLT1"].values)
            perm[c] = pe
        null[i] = pool_rho(m2_items(perm))["z"]
    return float((1 + (null >= obs).sum()) / (NPERM + 1))


frames = {c: v["data"] for c, v in cohorts.items()}
covs = {c: v["cov"] for c, v in cohorts.items()}
rng = np.random.default_rng(SEED)

main = {"M1": m1(frames, covs), "M2": pool_rho(m2_items(frames))}
main["M1"]["p_perm_one_sided"] = m1_perm(frames, covs, rng)
main["M2"]["p_perm_one_sided"] = m2_perm(frames, rng)
holm = multipletests([main["M1"]["p_perm_one_sided"], main["M2"]["p_perm_one_sided"]], method="holm")[1]
for key, q in zip(["M1", "M2"], holm):
    main[key]["holm_p"] = float(q)
    est = main[key]["fixed"] if key == "M1" else main[key]["pooled_rho"]
    main[key]["supported"] = bool(q < 0.05 and est > 0)
m3_items = [(partial(f[f.group == "PE"]), int((f.group == "PE").sum())) for f in frames.values()]
main["M3_descriptive"] = pool_rho(m3_items)
per_cohort = {c: {"M1_PE_vs_control (est, se, nPE, nCtl)": coef(frames[c], covs[c]),
                  "M2_rho_n": m2_items({c: frames[c]})[0], "M3_partial_rho_n": m3_items[i]}
              for i, c in enumerate(frames)}

# pre-specified sensitivity analyses (normal-theory p; they do not change the verdict)
sens = {}
for c in frames:
    rest = {k: v for k, v in frames.items() if k != c}
    sens[f"without_{c}"] = {"M1": m1(rest, covs), "M2": pool_rho(m2_items(rest))}
f = frames["GSE66273"]
noh = dict(frames, GSE66273=f[~f.HELLP.astype(bool)])
sens["GSE66273_without_HELLP"] = {"M1": m1(noh, covs), "M2": pool_rho(m2_items(noh))}
f = frames["GSE44711"]
early = dict(frames, GSE44711=f[(f.group == "PE") | (f.GA < 34)])
sens["GSE44711_controls_lt34wk"] = {"M1": m1(early, covs),
                                    "n_controls_GSE44711": int(((f.group == "control") & (f.GA < 34)).sum())}
sens["M1_no_covariates"] = m1(frames, {c: "" for c in frames})

result = {
    "cohorts": {c: {k: v for k, v in cohorts[c].items() if k != "data"} | {
        "group_counts": cohorts[c]["data"].group.value_counts().to_dict()} for c in cohorts},
    "main": main, "per_cohort": per_cohort, "sensitivity": sens,
    "settings": {"n_perm": NPERM, "seed": SEED},
}
(HERE / "M_confirmation.json").write_text(json.dumps(result, indent=2, default=float) + "\n")
print(json.dumps({"main": main, "per_cohort": per_cohort}, indent=2, default=float))
