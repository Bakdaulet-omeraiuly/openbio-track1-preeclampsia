"""Addendum L, POST HOC (written 2026-11-15 after all L results were seen; it can only limit the claim).

Reviewer concern: L-C4 pools PE and FGR-only placentas and removes only hyp_all, so part of the HIF2dep-FLT1 link
could come from the difference between diagnosis groups. Sensitivities (same scores and pooling as L3):
  (a) residualise on [1, hyp_all, PE indicator]   (diagnosis-adjusted)
  (b) PE placentas only, residualise on [1, hyp_all]
  (c) FGR-only placentas only
  (d) GSE35574 only: + GA, sex, labour, batch (the only cohort with these covariates), PE + FGR-only
Pre-stated reading: if (a) and (b) stay clearly positive, the link is not explained by diagnosis-group differences.
Usage: python l8_diagnosis_sensitivity_posthoc.py DIR   (same DIR as l3_confirm.py)
"""
import json
from pathlib import Path
import sys

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
src = (HERE / "l3_confirm.py").read_text().split("\ndef coef(")[0]  # loader part only: builds `cohorts`
ns = {"__file__": str(HERE / "l3_confirm.py")}
exec(compile(src, "l3_confirm_loader", "exec"), ns)
cohorts = ns["cohorts"]


def partial(frame, extra=()):
    cols = [np.ones(len(frame)), frame["hyp_all"].values] + [np.asarray(c, float) for c in extra]
    X = np.column_stack(cols)
    r = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    return float(stats.spearmanr(r(frame["hif2dep"].values), r(frame["FLT1"].values)).statistic)


def pool(items):
    rho = np.array([i[0] for i in items]); n = np.array([i[1] for i in items]); w = n - 3 - np.array([i[2] for i in items])
    z = float((w * np.arctanh(rho)).sum() / w.sum()); se = float(1 / np.sqrt(w.sum()))
    return {"pooled_rho": float(np.tanh(z)), "CI": [float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))],
            "p_one_sided": float(1 - stats.norm.cdf(z / se)), "per_cohort": {c: round(i[0], 3) for c, i in zip(cohorts, items)},
            "n": {c: int(i[1]) for c, i in zip(cohorts, items)}}


out = {}
items = []
for c, v in cohorts.items():
    d = v["data"]; dis = d[d.group.isin(["PE", "FGR"])]
    items.append((partial(dis), len(dis), 0))
out["baseline_L_C4"] = pool(items)
items = []
for c, v in cohorts.items():
    d = v["data"]; dis = d[d.group.isin(["PE", "FGR"])]
    items.append((partial(dis, [(dis.group == "PE").astype(float)]), len(dis), 1))
out["a_diagnosis_adjusted"] = pool(items)
for grp, key in (("PE", "b_PE_only"), ("FGR", "c_FGR_only")):
    items = []
    for c, v in cohorts.items():
        d = v["data"]; sub = d[d.group == grp]
        items.append((partial(sub), len(sub), 0))
    out[key] = pool(items)
d = cohorts["GSE35574"]["data"]; dis = d[d.group.isin(["PE", "FGR"])].dropna(subset=["GA"])
extra = [(dis.group == "PE").astype(float), dis.GA, (dis.sex == dis.sex.iloc[0]).astype(float),
         (dis.labour == dis.labour.iloc[0]).astype(float), (dis.batch == dis.batch.iloc[0]).astype(float)]
r = partial(dis, extra)
out["d_GSE35574_full_covariates"] = {"rho": r, "n": int(len(dis)), "unadjusted_same_set": partial(dis)}
(HERE / "L8_diagnosis_sensitivity_posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=1))
