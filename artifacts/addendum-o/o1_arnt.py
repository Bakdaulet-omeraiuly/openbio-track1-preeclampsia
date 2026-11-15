"""Addendum O1 (plan in plan.md): ARNT knockdown vs PT2385-defined HIF2dep program, GSE80995.

Usage: python o1_arnt.py GSE80995_series_matrix.txt.gz GPL10558.txt
Inputs: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE80nnn/GSE80995/matrix/GSE80995_series_matrix.txt.gz
        https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10558&targ=self&form=text&view=data
"""

import itertools
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
src = (ROOT / "artifacts/addendum-l/l3_confirm.py").read_text()
ns = {"__file__": str(ROOT / "artifacts/addendum-l/l3_confirm.py"), "__name__": "l3"}
sys_argv = sys.argv[:]
sys.argv = [sys.argv[0], "."]
exec(compile(src[: src.index("def scores(")], "l3_prefix", "exec"), ns)
sys.argv = sys_argv
read_series, read_gpl, gene_matrix, HYP, HIF2 = ns["read_series"], ns["read_gpl"], ns["gene_matrix"], ns["HYP"], ns["HIF2"]

m, e, log_applied = read_series(Path(sys.argv[1]))
genes = gene_matrix(e, read_gpl(Path(sys.argv[2]), "Symbol"))
title = m["!Sample_title"]
info = pd.DataFrame({"sh": m["treatment"], "o2": m["oxygen"].str.startswith("Hypoxia").map({True: "hyp", False: "norm"}),
                     "donor": title.str.extract(r"_(\d+P)_")[0]}, index=m.index)
cult = info[info.sh.isin(["shScramble", "shARNT"])]
assert len(cult) == 20 and cult.donor.nunique() == 5
g = genes[cult.index]
g = g[g.notna().all(axis=1)]
z = g.sub(g.mean(axis=1), axis=0).div(g.std(axis=1), axis=0).T


def score(gene_list):
    cols = sorted(set(gene_list) & set(z.columns))
    return z[cols].mean(axis=1), len(cols)


def delta(s, o2):
    d = {}
    for dn in sorted(cult.donor.unique()):
        sel = cult[(cult.donor == dn) & (cult.o2 == o2)]
        d[dn] = float(s[sel.index[sel.sh == "shARNT"]].iloc[0] - s[sel.index[sel.sh == "shScramble"]].iloc[0])
    return pd.Series(d)


def signflip_p_less(x):
    x = np.asarray(x); obs = x.mean()
    null = [np.mean(x * np.array(sg)) for sg in itertools.product([1, -1], repeat=len(x))]
    return float(np.mean(np.array(null) <= obs + 1e-12))


res = {"log2_applied_by_reader": log_applied, "n_genes_after_complete_filter": int(z.shape[1])}
single = {}
for gname in ["ARNT", "EPAS1", "HIF1A", "FLT1"]:
    if gname in z.columns:
        single[gname] = {"delta_hyp": delta(z[gname], "hyp").to_dict(), "delta_norm": delta(z[gname], "norm").to_dict()}
res["single_genes"] = single
hyp_s, n_hyp = score(HYP)
hif_s, n_hif = score(HIF2)
res["n_genes"] = {"hyp_all": n_hyp, "hif2dep": n_hif}
dA = delta(z["ARNT"], "hyp") if "ARNT" in z.columns else None
res["gate_ARNT_delta_hyp_mean"] = None if dA is None else float(dA.mean())
res["gate_pass"] = bool(dA is not None and dA.mean() < 0)
dH = delta(hyp_s, "hyp")
res["positive_control_hyp_all"] = {"delta_mean": float(dH.mean()), "per_donor": dH.to_dict(), "p_signflip": signflip_p_less(dH)}
res["positive_control_pass"] = bool(dH.mean() < 0)
d1 = delta(hif_s, "hyp")
d1n = delta(hif_s, "norm")
res["O1_primary"] = {"delta_mean": float(d1.mean()), "per_donor": d1.to_dict(), "p_signflip_one_sided": signflip_p_less(d1),
                     "n_donors_negative": int((d1 < 0).sum())}
res["O1_specificity_hyp_minus_norm"] = float((d1 - d1n).mean())
res["O1_delta_norm_mean"] = float(d1n.mean())
res["O1_beyond_hyp_all"] = float((d1 - dH).mean())
if "FLT1" in z.columns:
    dF = delta(z["FLT1"], "hyp"); res["FLT1_delta_hyp"] = {"mean": float(dF.mean()), "per_donor": dF.to_dict()}

mean_expr = g.mean(axis=1)
pool = sorted(set(z.columns) - HYP - HIF2 - {"FLT1"})
hif_in = sorted(HIF2 & set(z.columns))
dec = pd.qcut(mean_expr.reindex(sorted(z.columns)), 10, labels=False)
need = dec.reindex(hif_in).value_counts()
by = {d: np.array([x for x in pool if dec.get(x) == d]) for d in range(10)}
rng = np.random.default_rng(0)
null = np.array([delta(z[list(np.concatenate([rng.choice(by[d], n, replace=False) for d, n in need.items()]))].mean(axis=1), "hyp").mean()
                 for _ in range(1000)])
res["O1_calibration"] = {"null_median": float(np.median(null)), "null_5th": float(np.percentile(null, 5)),
                         "percentile_of_observed": float((null < d1.mean()).mean() * 100)}
ok = res["gate_pass"] and res["positive_control_pass"]
res["reading"] = ("uninterpretable (gate or positive control failed)" if not ok else
                  ("PT2385-defined program also depends on the HIF complex genetically"
                   if (d1.mean() < 0 and res["O1_primary"]["p_signflip_one_sided"] <= 0.05 and d1.mean() < res["O1_calibration"]["null_5th"])
                   else "not shown"))
(HERE / "O1_results.json").write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps({k: v for k, v in res.items() if k != "single_genes"}, indent=2))
print("FLT1/ARNT/EPAS1 hyp deltas:", {k: round(float(np.mean(list(v["delta_hyp"].values()))), 3) for k, v in single.items()})
