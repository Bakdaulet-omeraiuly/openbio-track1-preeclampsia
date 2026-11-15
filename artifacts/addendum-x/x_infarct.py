"""Addendum X: placental infarction vs HIF2dep program / FLT1 / generic hypoxia in ELGAN (GSE154829). Rules in plan.md.

  python x_infarct.py fetch DIR
  python x_infarct.py run DIR
"""

import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE154nnn/GSE154829/suppl/GSE154829_mRNA_countdata.csv.gz"
NNULL, SEED = 1000, 0
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
META = pd.read_csv(HERE / "samples.tsv", sep="\t", dtype=str)
FORMULA = "y ~ infarct + C(sex) + C(chorio) + C(fvni) + C(smoke) + mage + bmi"


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    out = d / URL.rsplit("/", 1)[1]
    if not out.exists():
        urllib.request.urlretrieve(URL, out)
    print("ok", out.name)


def load(d):
    c = pd.read_csv(d / URL.rsplit("/", 1)[1], index_col=0)
    c = c.select_dtypes(include=[np.number])
    col = {}
    for _, r in META.iterrows():
        hits = [x for x in c.columns if x in (r.title, r.gsm, r.title.replace("Subject_", ""))]
        assert len(hits) == 1, (r.title, hits[:3])
        col[r.gsm] = hits[0]
    c = c[[col[g] for g in META.gsm]]
    c.columns = META.gsm.tolist()
    c.index = [str(i).split("$", 1)[1] if "$" in str(i) else str(i) for i in c.index]
    if "FLT1" not in c.index:
        mp = pd.read_csv(ROOT / "data_external/ensembl_gene_id_to_name_GRCh37.87.tsv", sep="\t", header=None, names=["id", "sym"])
        mp = mp.drop_duplicates("id").set_index("id")["sym"]
        sym = pd.Series(pd.Index(c.index).astype(str).str.split(".").str[0], index=c.index).map(mp)
        c = c[sym.notna().values]; c.index = sym.dropna().values
    c = c.groupby(c.index).sum()
    assert "FLT1" in c.index
    return c


def prep(c, sel):
    sub = c[META.gsm[sel]]
    cpm = sub / sub.sum() * 1e6
    lg = np.log2(cpm[(cpm >= 1).sum(axis=1) >= 0.1 * cpm.shape[1]] + 1)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
    return lg, z


def score(z, genes):
    s = z.loc[z.index.intersection(sorted(genes))].mean()
    return (s - s.mean()) / s.std()


def model(y, m, formula=FORMULA, cluster=None):
    df = m.assign(y=y.values)
    for k in ("infarct", "mage", "bmi"):
        df[k] = pd.to_numeric(df[k], errors="coerce")
    df = df.dropna(subset=["y", "infarct", "mage", "bmi", "sex", "chorio", "fvni", "smoke"])
    df = df[~df[["sex", "chorio", "fvni", "smoke"]].isin(["NA", "2"]).any(axis=1)]
    if cluster is None:
        r = smf.ols(formula, data=df).fit(cov_type="HC3")
    else:
        r = smf.ols(formula, data=df).fit(cov_type="cluster", cov_kwds={"groups": pd.factorize(df[cluster])[0]})
    return float(r.params["infarct"]), float(r.bse["infarct"]), float(1 - stats.norm.cdf(r.params["infarct"] / r.bse["infarct"])), int(df.infarct.sum()), int((1 - df.infarct).sum())


def run(c, sel, formula=FORMULA, cluster=None, null=True):
    m = META[sel].set_index("gsm").copy()
    m["mother"] = m.title.str[:-1]
    lg, z = prep(c, sel)
    sc = {"hif2dep": score(z, HIF2), "hyp_only": score(z, HYP - GLY), "FLT1": (z.loc["FLT1"] - z.loc["FLT1"].mean()) / z.loc["FLT1"].std()}
    out = {k: dict(zip(("coef", "se", "p_one_sided", "n_infarct", "n_none"), model(v, m, formula, cluster))) for k, v in sc.items()}
    if not null:
        return out
    dec = pd.qcut(lg.mean(axis=1).rank(method="first"), 10, labels=False)
    pool = dec[~dec.index.isin(HYP | HIF2 | {"FLT1"})]
    for k, genes in (("hif2dep", HIF2), ("hyp_only", HYP - GLY)):
        need = dec[z.index.intersection(sorted(genes))].value_counts()
        rng = np.random.default_rng(SEED)
        nl = np.array([model(score(z, set(np.concatenate([rng.choice(pool.index[pool == q], size=n, replace=False)
                                                            for q, n in need.items()]))), m, formula, cluster)[0]
                       for _ in range(NNULL)])
        out[k].update({"random_median": float(np.median(nl)), "random_p95": float(np.percentile(nl, 95)),
                       "random_percentile": float((nl < out[k]["coef"]).mean() * 100)})
    same = pool.index[pool == dec["FLT1"]]
    rng = np.random.default_rng(SEED)
    zf = lambda g: (z.loc[g] - z.loc[g].mean()) / z.loc[g].std()
    nl = np.array([model(zf(g), m, formula, cluster)[0] for g in rng.choice(same, size=min(NNULL, len(same)), replace=False)])
    out["FLT1"].update({"random_median": float(np.median(nl)), "random_p95": float(np.percentile(nl, 95)),
                        "random_percentile": float((nl < out["FLT1"]["coef"]).mean() * 100)})
    return out


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    c = load(d)
    valid = META.infarct.isin(["0", "1"])
    single = valid & (META.fused2 == "0")
    res = {"primary_singletons": run(c, single.values)}
    p = res["primary_singletons"]
    holm = multipletests([p["hif2dep"]["p_one_sided"], p["FLT1"]["p_one_sided"]], method="holm")[1]
    for k, q in zip(("hif2dep", "FLT1"), holm):
        p[k]["holm_p"] = float(q)
        p[k]["supported"] = bool(q < 0.05 and p[k]["coef"] > p[k]["random_p95"])
    res["sens_a_all_births_cluster_by_mother"] = run(c, valid.values, cluster="mother", null=False)
    res["sens_b_no_inflammation_covariates"] = run(c, single.values, formula="y ~ infarct + C(sex) + C(smoke) + mage + bmi", null=False)
    (HERE / "X_results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
