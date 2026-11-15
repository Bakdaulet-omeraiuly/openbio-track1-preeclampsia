"""Addendum V, POST HOC (written after V1-V4 were seen; it can only limit the claim).

Question: does EPAS1 loss block syncytial differentiation, so that FLT1 falls because the cells are less
STB-like, rather than through HIF-2 regulation?

At 3% O2, null vs wild type, per study: STB/syncytial markers (CGA, CGB3, CYP19A1, ERVW-1, ERVFRD-1, GCM1, PSG4,
KISS1, LGALS13, HSD3B1, SDC1) and an STB marker score (mean z). Then FLT1 adjusted for the STB score (OLS) in
EPAS1-null vs wild type.
Pre-stated reading: if the STB score falls strongly with EPAS1 loss and FLT1's drop vanishes after adjustment,
the FLT1 effect may be a differentiation effect. If FLT1 still falls after adjustment, it is not explained by this
marker score. With n = 12 the adjustment is weak; this is descriptive.
Usage: python v2_differentiation_posthoc.py DIR
"""

import json
from pathlib import Path
import sys

import numpy as np
import statsmodels.formula.api as smf

sys.path.insert(0, str(Path(__file__).resolve().parent))
import v_isoform as V  # noqa: E402

MARKERS = ["CGA", "CGB3", "CYP19A1", "ERVW-1", "ERVFRD-1", "GCM1", "PSG4", "KISS1", "LGALS13", "HSD3B1", "SDC1"]
d = Path(sys.argv[1])
out = {}
for study in V.URL:
    c, meta = V.load(d, study)
    cpm = c / c.sum() * 1e6
    lg = np.log2(cpm[(cpm >= 1).sum(axis=1) >= 3] + 1)
    meta = meta.set_index("col")
    hyp = meta.index[meta.condition == "hypoxia"]
    sub = lg[hyp]
    null = (meta.loc[hyp, "gene"] != "Not applicable").values
    present = [m for m in MARKERS if m in sub.index]
    z = sub.loc[present].sub(sub.loc[present].mean(1), axis=0).div(sub.loc[present].std(1), axis=0)
    stb = z.mean()
    res = {"markers_null_minus_wt_log2": {m: float(sub.loc[m, null].mean() - sub.loc[m, ~null].mean()) for m in present},
           "stb_score_null_minus_wt": float(stb[null].mean() - stb[~null].mean())}
    df = np.array([sub.loc["FLT1"].values, null.astype(float), stb.values]).T
    import pandas as pd
    frame = pd.DataFrame(df, columns=["FLT1", "null", "stb"])
    fit = smf.ols("FLT1 ~ null + stb", data=frame).fit()
    res["FLT1_null_effect_adjusted_for_stb"] = {"coef": float(fit.params["null"]), "se": float(fit.bse["null"]),
                                               "p_two_sided": float(fit.pvalues["null"])}
    res["FLT1_null_effect_unadjusted"] = float(frame.FLT1[frame.null == 1].mean() - frame.FLT1[frame.null == 0].mean())
    out[study] = res
(Path(__file__).resolve().parent / "V_differentiation_posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
for s, r in out.items():
    print("==", s, V.TARGET[s], "STB score null-wt %.2f" % r["stb_score_null_minus_wt"])
    print("   markers", {k: round(v, 2) for k, v in r["markers_null_minus_wt_log2"].items()})
    print("   FLT1 unadjusted %.2f | adjusted for STB score %.2f (se %.2f, p %.3g)" % (
        r["FLT1_null_effect_unadjusted"], r["FLT1_null_effect_adjusted_for_stb"]["coef"],
        r["FLT1_null_effect_adjusted_for_stb"]["se"], r["FLT1_null_effect_adjusted_for_stb"]["p_two_sided"]))
