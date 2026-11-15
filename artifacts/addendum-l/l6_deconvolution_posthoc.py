"""POST HOC (written after L-C4, L4 and L5 were seen; committed before this script was run).

Question: does L-C4 (HIF2dep score vs FLT1 | hyp_all, PE + FGR-only) survive adjustment for cell
proportions estimated by deconvolution against a placental single-cell reference? This can only
LIMIT the C4 claim. It is NOT the pre-registered new-cohort test named in results.md: the cohorts
are the same three already used for C4.

Reference: the datathon single-cell pseudobulk (track01_pseudobulk_counts.csv, 14 cell classes,
26 donors; villous + decidual compartments as supplied). Method = the existing Addendum B5 rule,
unchanged except where the arrays force a change:
  * class profile = mean CPM over all pseudobulks of the class (as B5);
  * signature = per class, the top 50 genes by log2(CPM+1) specificity (own class minus max of the
    others) with own log2(CPM+1) > 3 (as B5);
  * CHANGES: candidate genes are restricted to genes measured on all three array platforms, and
    HIF2dep genes, FLT1 and HALLMARK_HYPOXIA genes are removed before selection, so the
    signature cannot absorb the score, the outcome or the adjuster;
  * bulk input = 2**(log2 array value) per gene, NNLS per placenta, proportions renormalised to 1
    (as B5). Arrays are not linear in abundance, so these are rough relative estimates.

Validation (reported, no threshold): Spearman(estimated STB proportion, L5 STB marker score).

Primary (one): partial Spearman(HIF2dep, FLT1 | hyp_all + STB + villous CTB proportions), the B5
covariate pair, pooled by Fisher z with weights n - 3 - k (k = covariates beyond hyp_all = 2).
Secondary (one): | hyp_all + STB + CTB + EVT proportions (k = 3).
Calibration: random 100-gene sets as in L5 (decile-matched, seed 0, 1,000 draws), same primary
adjustment.
Reading rule (fixed before the run): pooled 95% CI includes 0 -> C4 cannot be separated from
estimated composition; clearly > 0 and above the adjusted null 95th percentile -> "not explained by
deconvolution-estimated STB/CTB proportions" (still not "free of composition": array NNLS is rough,
and composition may mediate).

Usage: python l6_deconvolution_posthoc.py DIR
"""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import nnls

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
DATA = ROOT / "openbio-datathon-2026/track-01/data"
src = (HERE / "l5_stb_adjustment_posthoc.py").read_text()
ns = {"__file__": str(HERE / "l5_stb_adjustment_posthoc.py"), "__name__": "l5"}
sys.argv = [sys.argv[0], sys.argv[1]]
exec(compile(src[: src.index("K = {")], "l5_prefix", "exec"), ns)  # loading, marker_score, partial, pool
genes_by, groups, Z, mean_expr = ns["genes_by"], ns["groups"], ns["Z"], ns["mean_expr"]
HYP, HIF2, STB = ns["HYP"], ns["HIF2"], ns["STB"]
marker_score, partial, pool = ns["marker_score"], ns["partial"], ns["pool"]

pb = pd.read_csv(DATA / "extension/pseudobulk/track01_pseudobulk_counts.csv", index_col=0)
pbm = pd.read_csv(DATA / "extension/pseudobulk/track01_pseudobulk_metadata.csv", index_col=0).loc[pb.columns]
classes = sorted(pbm.cell_class.unique())
cpm = pb / pb.sum() * 1e6
prof = pd.DataFrame({k: cpm[pbm.index[pbm.cell_class == k]].mean(axis=1) for k in classes})
common = set.intersection(*[set(g.index) for g in genes_by.values()])
cand = sorted((set(prof.index) & common) - HIF2 - HYP - {"FLT1"})
lprof = np.log2(prof.loc[cand] + 1)
sig = set()
for k in classes:
    spec = lprof[k] - lprof.drop(columns=k).max(axis=1)
    sig |= set(spec[lprof[k] > 3].sort_values(ascending=False).index[:50])
sig = sorted(sig)
P = prof.loc[sig].values

props = {}
for c, g in genes_by.items():
    bulk = 2 ** g.loc[sig]
    est = pd.DataFrame({s: nnls(P, bulk[s].values)[0] for s in bulk.columns}, index=classes).T
    props[c] = est.div(est.sum(axis=1), axis=0)

COVS = {"primary": ["syncytiotrophoblast", "villous_cytotrophoblast"],
        "secondary": ["syncytiotrophoblast", "villous_cytotrophoblast", "extravillous_trophoblast"]}


def run(gene_list, which):
    out = []
    for c in Z:
        sel = groups[c].isin(["PE", "FGR"]).values
        z = Z[c][sel]
        score = z[sorted(g for g in gene_list if g in z.columns)].mean(axis=1)
        hyp = z[sorted(g for g in HYP if g in z.columns)].mean(axis=1)
        covs = [hyp.values] + ([props[c].loc[z.index, k].values for k in COVS[which]] if which != "base" else [])
        out.append((partial(score.values, z["FLT1"].values, covs), int(sel.sum())))
    return out


baseline = pool(run(HIF2, "base"), 0)
assert [round(baseline["per_cohort"][c]["rho"], 2) for c in Z] == [0.34, 0.44, 0.49]
primary = pool(run(HIF2, "primary"), len(COVS["primary"]))
secondary = pool(run(HIF2, "secondary"), len(COVS["secondary"]))

summary = {}
for c in Z:
    sel = groups[c].isin(["PE", "FGR"]).values
    z = Z[c][sel]
    pr = props[c].loc[z.index]
    stb_mk, _ = marker_score(z, STB)
    h = z[sorted(g for g in HIF2 if g in z.columns)].mean(axis=1)
    summary[c] = {"median_props_PE_FGR": pr.median().round(3).to_dict(),
                  "rho_STBprop_STBmarker": float(stats.spearmanr(pr["syncytiotrophoblast"], stb_mk).statistic),
                  "rho_STBprop_HIF2dep": float(stats.spearmanr(pr["syncytiotrophoblast"], h).statistic),
                  "rho_STBprop_FLT1": float(stats.spearmanr(pr["syncytiotrophoblast"], z["FLT1"]).statistic),
                  "rho_CTBprop_HIF2dep": float(stats.spearmanr(pr["villous_cytotrophoblast"], h).statistic)}

pool_genes = sorted(common - HYP - HIF2 - {"FLT1"} - set(sig))
hif2_common = sorted(HIF2 & common)
deciles = pd.qcut(mean_expr["GSE35574"].reindex(sorted(common)), 10, labels=False)
need = deciles.reindex(hif2_common).value_counts()
by_dec = {d: np.array([g for g in pool_genes if deciles.get(g) == d]) for d in range(10)}
obs_common = pool(run(hif2_common, "primary"), 2)["pooled_rho"]
rng = np.random.default_rng(0)
null = np.array([pool(run(list(np.concatenate([rng.choice(by_dec[d], n, replace=False) for d, n in need.items()])), "primary"),
                      2)["pooled_rho"] for _ in range(1000)])

reading = ("C4 cannot be separated from estimated composition" if primary["CI_low"] <= 0 else
           ("C4 not explained by deconvolution-estimated STB/CTB proportions" if obs_common > np.percentile(null, 95)
            else "adjusted rho positive but not above the adjusted random-set null"))
result = {"posthoc": True, "n_signature_genes": len(sig), "signature_per_class_cap": 50,
          "baseline_L3_C4": baseline, "primary_hyp_STB_CTB": primary, "secondary_hyp_STB_CTB_EVT": secondary,
          "composition_summary": summary,
          "calibration_adjusted": {"observed_common_genes": obs_common, "null_median": float(np.median(null)),
                                   "null_95th": float(np.percentile(null, 95)),
                                   "empirical_p": float(((null >= obs_common).sum() + 1) / 1001)},
          "reading": reading}
(HERE / "L6_deconvolution_posthoc.json").write_text(json.dumps(result, indent=2) + "\n")
for c, p in props.items():
    p.to_csv(HERE / f"L6_estimated_proportions_{c}.csv")
print(json.dumps(result, indent=2))
