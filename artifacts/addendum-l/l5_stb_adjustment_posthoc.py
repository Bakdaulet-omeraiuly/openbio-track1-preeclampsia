"""POST HOC (written after L-C4 and L4 were seen; committed before this script was run).

Question: is the L-C4 link (HIF2dep score vs FLT1 | hyp_all, PE + FGR-only) explained by
syncytiotrophoblast (STB) content? This can only LIMIT the C4 claim.

STB covariate: the official datathon marker matrix (track01_marker_matrix.csv), class
"syncytiotrophoblast" = CGA, CGB, ERVW-1, PSG1. None overlap HIF2dep; FLT1 is not among them.
CGB is a multi-copy family whose symbol differs by platform, so a fixed alias rule is used:
the CGB value is the mean z of whichever of CGB, CGB3, CGB5, CGB8 the platform measures.
Measured per platform (symbol check, done before writing this script):
  GPL6102 (GSE35574): CGA, CGB5, CGB8, PSG1 (no ERVW-1) -> 3 markers
  GPL6244 (GSE24129): CGA, CGB, CGB5, ERVW-1, PSG1     -> 4 markers
  GPL20844 (GSE147776): CGA, CGB, ERVW-1, PSG1         -> 4 markers
With 3-4 genes this is a noisy proxy for STB content.

Primary (one): partial Spearman(HIF2dep, FLT1 | hyp_all + STB) in PE + FGR-only, per cohort,
pooled by Fisher z with weights n - 3 - k (k = number of covariates).
Secondary (one): | hyp_all + STB + EVT (EVT markers CEACAM1, HLA-G, ITGA1, MMP2, as measured).
Baseline must reproduce L3 (0.34 / 0.44 / 0.49, pooled 0.39); asserted below.
Calibration: the L4 random-set null rerun under the same primary adjustment (same seed, same
decile matching, same 1,000 draws).

Reading rule (fixed before the run):
  * adjusted pooled 95% CI includes 0 -> C4 cannot be separated from STB marker content;
  * adjusted pooled rho clearly > 0 and above the adjusted null 95th percentile -> "not explained
    by this STB marker score" (NOT "free of composition": 3-4 markers, and STB content may be a
    mediator, since HIF-2 is reported to impair syncytialisation).

Usage: python l5_stb_adjustment_posthoc.py DIR
"""

import json
from pathlib import Path
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

mk = pd.read_csv(ROOT / "openbio-datathon-2026/track-01/data/interpretation/track01_marker_matrix.csv")
STB = sorted(mk.loc[mk.cell_class == "syncytiotrophoblast", "gene_symbol"].astype(str))
EVT = sorted(mk.loc[mk.cell_class == "extravillous_trophoblast", "gene_symbol"].astype(str))
assert STB == ["CGA", "CGB", "ERVW-1", "PSG1"] and not (set(STB) | set(EVT)) & (HIF2 | {"FLT1"})
CGB_ALIASES = ["CGB", "CGB3", "CGB5", "CGB8"]

# Same loading as l3_confirm.py / l4_calibration_posthoc.py
genes_by, groups = {}, {}
m, e, _ = read_series(DIR / "GSE35574_series_matrix.txt.gz")
base = m["!Sample_title"].str.replace(r"\.rep$", "", regex=True)
g = gene_matrix(e, read_gpl(DIR / "GPL6102.txt", "Symbol")).T.groupby(base.values).mean().T
mm = m.assign(base=base.values).drop_duplicates("base").set_index("base").loc[g.columns]
genes_by["GSE35574"] = g
groups["GSE35574"] = mm["classification"].map({"PE": "PE", "IUGR": "FGR", "CONTROL": "control"})
m, e, _ = read_series(DIR / "GSE24129_series_matrix.txt.gz")
genes_by["GSE24129"] = gene_matrix(e, read_gpl(DIR / "GPL6244.txt", "gene_assignment"))
groups["GSE24129"] = m["disease state"].map(lambda v: "PE" if v.startswith("Pre-eclampsia") else ("FGR" if v.startswith("Fetal growth") else "control"))
m, e, _ = read_series(DIR / "GSE147776_series_matrix.txt.gz")
genes_by["GSE147776"] = gene_matrix(e, read_gpl(DIR / "GPL20844.txt", "GENE_SYMBOL"))
groups["GSE147776"] = m["disease state"].map({"Preeclampsia with severe features": "PE",
                                              "Intrauterine growth restriction associated with severe preeclampsia": "PE",
                                              "Intrauterine growth restriction": "FGR", "Control group": "control"})

Z, mean_expr, used = {}, {}, {}
for c, genes in genes_by.items():
    mean_expr[c] = genes.mean(axis=1)
    Z[c] = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T


def marker_score(z, markers):
    cols, parts = [], []
    for gname in markers:
        if gname == "CGB":
            al = [a for a in CGB_ALIASES if a in z.columns]
            if al:
                parts.append(z[al].mean(axis=1)); cols.append("CGB(" + "+".join(al) + ")")
        elif gname in z.columns:
            parts.append(z[gname]); cols.append(gname)
    return pd.concat(parts, axis=1).mean(axis=1), cols


def partial(x, y, covs):
    X = np.column_stack([np.ones(len(x))] + covs)
    resid = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    return float(stats.spearmanr(resid(x), resid(y)).statistic)


def pool(items, k):
    rho = np.array([i[0] for i in items]); n = np.array([i[1] for i in items]); w = n - 3 - k
    zz = float((w * np.arctanh(rho)).sum() / w.sum()); se = float(1 / np.sqrt(w.sum()))
    return {"pooled_rho": float(np.tanh(zz)), "CI_low": float(np.tanh(zz - 1.96 * se)),
            "CI_high": float(np.tanh(zz + 1.96 * se)), "p_one_sided": float(1 - stats.norm.cdf(zz / se)),
            "per_cohort": {c: {"rho": r, "n": int(nn)} for c, (r, nn) in zip(Z, items)}}


def run(gene_list, which):
    """which: 'base' | 'stb' | 'stb_evt'. Gene z is over all placentas, as in L3 (set-score rescaling is rank-neutral)."""
    out = []
    for c in Z:
        sel = groups[c].isin(["PE", "FGR"]).values
        z = Z[c][sel]
        gl = sorted(g for g in gene_list if g in z.columns)
        score = z[gl].mean(axis=1)
        hyp = z[sorted(g for g in HYP if g in z.columns)].mean(axis=1)
        covs = [hyp.values]
        if which in ("stb", "stb_evt"):
            s, cols = marker_score(z, STB); covs.append(s.values); used.setdefault(c, {})["STB"] = cols
        if which == "stb_evt":
            ev, cols = marker_score(z, EVT); covs.append(ev.values); used.setdefault(c, {})["EVT"] = cols
        out.append((partial(score.values, z["FLT1"].values, covs), int(sel.sum())))
    return out


K = {"base": 1, "stb": 2, "stb_evt": 3}
# Baseline reproduction of L3 C4 (per-cohort measured HIF2dep genes, | hyp_all)
baseline = pool(run(HIF2, "base"), 0)  # L3 used n - 3 weights
got = [round(baseline["per_cohort"][c]["rho"], 2) for c in Z]
assert got == [0.34, 0.44, 0.49], got
assert round(baseline["pooled_rho"], 2) == 0.39, baseline["pooled_rho"]

primary = pool(run(HIF2, "stb"), K["stb"] - 1)
secondary = pool(run(HIF2, "stb_evt"), K["stb_evt"] - 1)

# STB score vs HIF2dep / FLT1 (descriptive), PE + FGR-only
desc = {}
for c in Z:
    sel = groups[c].isin(["PE", "FGR"]).values
    z = Z[c][sel]
    s, _ = marker_score(z, STB)
    h = z[sorted(g for g in HIF2 if g in z.columns)].mean(axis=1)
    desc[c] = {"rho_STB_HIF2dep": float(stats.spearmanr(s, h).statistic),
               "rho_STB_FLT1": float(stats.spearmanr(s, z["FLT1"]).statistic)}

# Null under the same primary adjustment (same draws as L4)
common = set.intersection(*[set(z.columns) for z in Z.values()])
pool_genes = sorted(common - HYP - HIF2 - {"FLT1"} - set(STB) - set(EVT) - set(CGB_ALIASES))
hif2_common = sorted(HIF2 & common)
deciles = pd.qcut(mean_expr["GSE35574"].reindex(sorted(common)), 10, labels=False)
need = deciles.reindex(hif2_common).value_counts()
by_dec = {d: np.array([g for g in pool_genes if deciles.get(g) == d]) for d in range(10)}
obs_common = pool(run(hif2_common, "stb"), K["stb"] - 1)["pooled_rho"]
rng = np.random.default_rng(0)
null = np.array([pool(run(list(np.concatenate([rng.choice(by_dec[d], n, replace=False) for d, n in need.items()])), "stb"),
                      K["stb"] - 1)["pooled_rho"] for _ in range(1000)])

primary_ci_includes_0 = primary["CI_low"] <= 0
reading = ("C4 cannot be separated from STB marker content" if primary_ci_includes_0 else
           ("C4 not explained by this STB marker score" if obs_common > np.percentile(null, 95)
            else "adjusted rho positive but not above the adjusted random-set null"))
result = {"posthoc": True, "markers_used": used, "baseline_L3_C4": baseline,
          "primary_hyp_plus_STB": primary, "secondary_hyp_plus_STB_plus_EVT": secondary,
          "descriptive_STB_correlations": desc,
          "calibration_adjusted": {"observed_common_genes": obs_common, "n_common_hif2dep": len(hif2_common),
                                   "null_median": float(np.median(null)), "null_95th": float(np.percentile(null, 95)),
                                   "empirical_p": float(((null >= obs_common).sum() + 1) / 1001)},
          "reading": reading}
(HERE / "L5_stb_adjustment_posthoc.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
