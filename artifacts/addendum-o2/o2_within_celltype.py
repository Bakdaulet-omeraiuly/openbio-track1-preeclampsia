"""Addendum O2: HIF2dep–FLT1 partial correlation (given hyp_all) within single trophoblast cell types.

Rules in plan.md. Data: datathon pseudobulk (openbio-datathon-2026/track-01/data/extension/pseudobulk/).
Usage: python o2_within_celltype.py
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
PB = ROOT / "openbio-datathon-2026/track-01/data/extension/pseudobulk"
NPERM, SEED, NNULL, NULL_SEED = 10_000, 42, 1000, 0

gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
assert "FLT1" not in HYP and "FLT1" not in HIF2

counts = pd.read_csv(PB / "track01_pseudobulk_counts.csv", index_col=0)
meta = pd.read_csv(PB / "track01_pseudobulk_metadata.csv", index_col=0)


def prep(cell_class, donors=None):
    m = meta[(meta.cell_class == cell_class) & (meta.n_cells >= 20)]
    assert m.donor_id.is_unique
    if donors is not None:
        m = m[m.donor_id.isin(donors)]
    c = counts[m.index]
    c = c.loc[(c >= 10).sum(axis=1) >= 3]
    lc = np.log2(c / c.sum() * 1e6 + 1)
    z = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1).replace(0, np.nan), axis=0).dropna()
    return m, lc, z


def resid(y, X):
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def partial(a, b, X):
    return float(stats.spearmanr(resid(a, X), resid(b, X)).statistic)


def design(m, hyp, clinical=False):
    cols = [np.ones(len(m)), hyp]
    if clinical:
        cols += [m.ga_delivery_weeks.astype(float).values, (m.pe_status == "PE").astype(float).values,
                 (m.onset_group == "early").astype(float).values]
    return np.column_stack(cols)


def run(cell_class, clinical=False, pe_only=False, perm=False, null=False):
    m, lc, z = prep(cell_class)
    if pe_only:
        m = m[m.pe_status == "PE"]
        m, lc, z = prep(cell_class, donors=set(m.donor_id))
    hif = z.loc[z.index.intersection(sorted(HIF2))].mean().values
    hyp = z.loc[z.index.intersection(sorted(HYP))].mean().values
    flt = lc.loc["FLT1"].values
    X = design(m, hyp, clinical)
    rho = partial(hif, flt, X)
    out = {"cell_class": cell_class, "n": int(len(m)), "groups": m.clinical_group.value_counts().to_dict(),
           "n_hif2dep_genes": int(len(z.index.intersection(sorted(HIF2)))),
           "n_hyp_genes": int(len(z.index.intersection(sorted(HYP)))), "partial_rho": rho}
    n = len(m); zf = np.arctanh(rho); se = 1 / np.sqrt(n - 3 - (X.shape[1] - 2))
    out["approx_CI"] = [float(np.tanh(zf - 1.96 * se)), float(np.tanh(zf + 1.96 * se))]
    if perm:
        rng = np.random.default_rng(SEED)
        ra, rb = resid(hif, X), resid(flt, X)
        null_r = np.array([stats.spearmanr(ra, rng.permutation(rb)).statistic for _ in range(NPERM)])
        out["p_perm_one_sided"] = float((1 + (null_r >= rho).sum()) / (NPERM + 1))
    if null:
        rng = np.random.default_rng(NULL_SEED)
        mean_expr = lc.loc[z.index].mean(axis=1)
        dec = pd.qcut(mean_expr.rank(method="first"), 10, labels=False)
        target = dec[z.index.intersection(sorted(HIF2))].value_counts()
        pool = dec[~dec.index.isin(HYP | HIF2 | {"FLT1"})]
        rs = []
        for _ in range(NNULL):
            genes = np.concatenate([rng.choice(pool.index[pool == d], size=k, replace=False) for d, k in target.items()])
            rs.append(partial(z.loc[genes].mean().values, flt, X))
        rs = np.array(rs)
        out["random_null"] = {"median": float(np.median(rs)), "p95": float(np.percentile(rs, 95)),
                              "percentile_of_observed": float((rs < rho).mean() * 100),
                              "empirical_p": float((1 + (rs >= rho).sum()) / (NNULL + 1))}
    return out


result = {
    "primary_STB": run("syncytiotrophoblast", perm=True),
    "S1_STB_clinical_covariates": run("syncytiotrophoblast", clinical=True),
    "S2_VCT": run("villous_cytotrophoblast"),
    "S2_EVT": run("extravillous_trophoblast"),
    "S3_STB_random_null": run("syncytiotrophoblast", null=True),
    "S4_STB_PE_only": run("syncytiotrophoblast", pe_only=True),
}
p = result["primary_STB"]
p["supported"] = bool(p["p_perm_one_sided"] < 0.05 and p["partial_rho"] > 0)
(HERE / "O2_results.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
