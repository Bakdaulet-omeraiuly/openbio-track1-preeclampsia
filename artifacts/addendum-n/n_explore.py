"""Addendum N (exploratory; plan in plan.md): clinical anchors of the hypoxia score in GSE75010.

Usage: python n_explore.py GSE75010_series_matrix.txt GPL6244.txt
Reuses the Addendum K loading/scoring code unchanged.
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
KS = ROOT / "artifacts/addendum-k/k_confirm_gse75010.py"
src = KS.read_text()
cut = src.index('d["sex"] = meta["infant gender"]') + len('d["sex"] = meta["infant gender"]')
ns = {"__file__": str(KS), "__name__": "k"}
exec(compile(src[:cut], "k_prefix", "exec"), ns)
d, meta = ns["d"].copy(), ns["meta"]


def num(col):
    return pd.to_numeric(meta[col].replace("NA", np.nan), errors="coerce")


d["labour"] = meta["attempted vaginal delivery"]
d["uPI"], d["umbPI"] = num("mean uterine pi"), num("mean umbilical pi")
d["sbp"], d["bwz"], d["pwz"] = num("maximum systolic bp"), num("newborn weight z-score"), num("placental weight z-score")
d["hellp"] = (meta["hellp diagnosis"] == "Yes").astype(float)
d["early"] = (d.ga < 34).astype(int)
ctrl = (d.PE == 0) & (d.IUGR == 0)


def ols_pe(frame):
    f = smf.ols("hyp_only ~ PE + ga + C(sex) + C(labour)", data=frame).fit()
    ci = f.conf_int().loc["PE"]
    return {"estimate": float(f.params["PE"]), "CI_low": float(ci[0]), "CI_high": float(ci[1]),
            "p_two_sided": float(f.pvalues["PE"]), "n_PE": int(frame.PE.sum()), "n_control": int((frame.PE == 0).sum())}


def partial(frame, x, y, sign):
    f = frame[[x, y, "ga"]].dropna()
    X = np.column_stack([np.ones(len(f)), f.ga.values])
    r = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    rho, p2 = stats.spearmanr(r(f[x].values), r(f[y].values))
    p1 = p2 / 2 if np.sign(rho) == sign else 1 - p2 / 2
    return {"rho": float(rho), "n": int(len(f)), "p_one_sided": float(p1), "predicted_sign": sign}


res = {}
sub = d[(d.PE == 1) | ctrl]
res["N1_early"] = ols_pe(sub[sub.early == 1])
res["N1_late"] = ols_pe(sub[sub.early == 0])
fi = smf.ols("hyp_only ~ PE * early + ga + C(sex) + C(labour)", data=sub).fit()
ci = fi.conf_int().loc["PE:early"]
res["N1_interaction_PE_x_early"] = {"estimate": float(fi.params["PE:early"]), "CI_low": float(ci[0]),
                                    "CI_high": float(ci[1]), "p_two_sided": float(fi.pvalues["PE:early"])}

res["N2_hyp_only_vs_uterinePI"] = partial(d, "hyp_only", "uPI", +1)
res["N2s_FLT1_vs_uterinePI"] = partial(d, "FLT1", "uPI", +1)
res["N2s_hyp_only_vs_umbilicalPI"] = partial(d, "hyp_only", "umbPI", +1)
res["N2s_uPI_by_group_n"] = {"PE": int(d.loc[d.PE == 1, "uPI"].notna().sum()),
                             "non-PE": int(d.loc[d.PE == 0, "uPI"].notna().sum())}
pe = d[d.PE == 1]
for key, col, sign in [("N3_sbp", "sbp", +1), ("N3_newborn_weight_z", "bwz", -1),
                       ("N3_placental_weight_z", "pwz", -1), ("N3_hellp", "hellp", +1)]:
    res[key] = partial(pe, "hyp_only", col, sign)
fam = ["N2_hyp_only_vs_uterinePI", "N3_sbp", "N3_newborn_weight_z", "N3_placental_weight_z", "N3_hellp"]
for k, q in zip(fam, multipletests([res[k]["p_one_sided"] for k in fam], method="holm")[1]):
    res[k]["holm_p"] = float(q)
res["exploratory"] = True
(HERE / "N_results.json").write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps(res, indent=2))
