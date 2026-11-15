"""Addendum P (plan in plan.md): GSE186257, HIF2dep-FLT1 replication with cell-state adjustment.

Usage: python p_confirm.py GSE186257_series_matrix.txt.gz GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz
Inputs: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186257/matrix/GSE186257_series_matrix.txt.gz
        https://ftp.ncbi.nlm.nih.gov/geo/series/GSE186nnn/GSE186257/suppl/GSE186257_Plac_samples_DESEq2_norm_filtered.txt.gz
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
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
STB, EVT = ["CGA", "CGB", "ERVW-1", "PSG1"], ["CEACAM1", "HLA-G", "ITGA1", "MMP2"]
CGB_ALIASES = ["CGB", "CGB3", "CGB5", "CGB8"]

meta = {}
with gzip.open(sys.argv[1], "rt") as h:
    for line in h:
        v = [x.strip('"') for x in line.rstrip("\n").split("\t")]
        if v[0] == "!Sample_title":
            meta["title"] = v[1:]
        elif v[0] == "!Sample_characteristics_ch1":
            meta[v[1].split(":")[0]] = [x.split(": ", 1)[-1] for x in v[1:]]
meta = pd.DataFrame(meta).set_index("title")
tab = pd.read_csv(sys.argv[2], sep="\t", index_col=0)
assert set(tab.columns) == set(meta.index) and tab.shape[1] == 44
tab = tab[meta.index]
sym = tab.index.str.split("$").str[1]
genes = np.log2(tab.clip(lower=0) + 1).groupby(sym.values).mean()
z = genes.sub(genes.mean(axis=1), axis=0).div(genes.std(axis=1), axis=0).T
z = z.loc[:, z.notna().all()]


def mscore(gl):
    cols = []
    for g in gl:
        if g == "CGB":
            al = [a for a in CGB_ALIASES if a in z.columns]
            if al:
                cols.append(z[al].mean(axis=1))
        elif g in z.columns:
            cols.append(z[g])
    return pd.concat(cols, axis=1).mean(axis=1)


def sscore(gs):
    return z[sorted(set(gs) & set(z.columns))].mean(axis=1)


d = pd.DataFrame({"hif2dep": sscore(HIF2), "hyp_all": sscore(HYP), "hyp_only": sscore(HYP - GLY),
                  "FLT1": z["FLT1"], "stb": mscore(STB), "evt": mscore(EVT)})
d["group"] = np.where(meta["disease"] == "Severe preeclampsia", "PE", np.where(meta["sga"] == "yes", "SGA", "control"))
d["sex"] = meta["gender"].values
n_genes = {"hif2dep": len(HIF2 & set(z.columns)), "hyp_all": len(HYP & set(z.columns)),
           "stb_markers": [g for g in STB + CGB_ALIASES if g in z.columns], "evt_markers": [g for g in EVT if g in z.columns]}
rng = np.random.default_rng(0)


def resid(frame, y, covs):
    X = np.column_stack([np.ones(len(frame))] + [frame[c].values for c in covs])
    return frame[y].values - X @ np.linalg.lstsq(X, frame[y].values, rcond=None)[0]


def partial(frame, x, y, covs, nperm=10000):
    rx, ry = resid(frame, x, covs), resid(frame, y, covs)
    rho = stats.spearmanr(rx, ry).statistic
    null = np.array([stats.spearmanr(rx, rng.permutation(ry)).statistic for _ in range(nperm)])
    n, k = len(frame), len(covs)
    se = 1 / np.sqrt(max(n - 3 - k, 1)); zz = np.arctanh(rho)
    return {"rho": float(rho), "n": int(n), "CI_low": float(np.tanh(zz - 1.96 * se)), "CI_high": float(np.tanh(zz + 1.96 * se)),
            "p_perm_one_sided": float(((null >= rho).sum() + 1) / (nperm + 1))}


dis = d[d.group.isin(["PE", "SGA"])]
pe = d[d.group == "PE"]
res = {"n_genes": n_genes, "group_counts": d.group.value_counts().to_dict()}
res["P1"] = partial(dis, "hif2dep", "FLT1", ["hyp_all"])
res["P2"] = partial(pe, "hyp_all", "FLT1", [])
for k, q in zip(["P1", "P2"], multipletests([res["P1"]["p_perm_one_sided"], res["P2"]["p_perm_one_sided"]], method="holm")[1]):
    res[k]["holm_p"] = float(q)
res["P1_cs_STB_EVT"] = partial(dis, "hif2dep", "FLT1", ["hyp_all", "stb", "evt"])
res["P1_PE_only"] = partial(pe, "hif2dep", "FLT1", ["hyp_all"])
res["desc_rho_hif2dep_stb_in_disease"] = float(stats.spearmanr(dis.hif2dep, dis.stb).statistic)
res["desc_rho_hif2dep_evt_in_disease"] = float(stats.spearmanr(dis.hif2dep, dis.evt).statistic)


def contrast(other):
    f = d[d.group.isin(["PE", other])].copy(); f["PE"] = (f.group == "PE").astype(int)
    m = smf.ols("hyp_only ~ PE + C(sex)", data=f).fit(); ci = m.conf_int().loc["PE"]
    return {"estimate": float(m.params["PE"]), "CI_low": float(ci[0]), "CI_high": float(ci[1]), "p_two_sided": float(m.pvalues["PE"])}


res["desc_hyp_only_PE_vs_control"] = contrast("control")
res["desc_hyp_only_PE_vs_SGA"] = contrast("SGA")

mean_expr = genes.reindex(z.columns).mean(axis=1)
pool = sorted(set(z.columns) - HYP - HIF2 - {"FLT1"})
hin = sorted(HIF2 & set(z.columns))
dec = pd.qcut(mean_expr.reindex(sorted(z.columns)), 10, labels=False)
need = dec.reindex(hin).value_counts()
by = {k: np.array([g for g in pool if dec.get(g) == k]) for k in range(10)}
rx_cov = ["hyp_all"]
null = []
for _ in range(1000):
    pick = np.concatenate([rng.choice(by[k], n, replace=False) for k, n in need.items()])
    tmp = dis.assign(rs=z.loc[dis.index, list(pick)].mean(axis=1))
    null.append(stats.spearmanr(resid(tmp, "rs", rx_cov), resid(tmp, "FLT1", rx_cov)).statistic)
null = np.array(null)
res["desc_P1_random_set_null"] = {"median": float(np.median(null)), "p95": float(np.percentile(null, 95)),
                                  "percentile_of_observed": float((null < res["P1"]["rho"]).mean() * 100)}
res["reading"] = {"P1": "replicates" if res["P1"]["holm_p"] <= 0.05 else "not replicated",
                  "P1_cs": "not explained by marker scores" if res["P1_cs_STB_EVT"]["CI_low"] > 0 else "cell state could explain it here"}
(HERE / "P_results.json").write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps(res, indent=2))
