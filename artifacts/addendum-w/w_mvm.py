"""Addendum W: MVM vs hypoxia program in GSE303463 (rules in plan.md).

  python w_mvm.py fetch DIR
  python w_mvm.py run DIR
"""

import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE303nnn/GSE303463/suppl/GSE303463_raw_counts.csv.gz"
NNULL, SEED = 1000, 0
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP, GLY = set(gmt["HALLMARK_HYPOXIA"]), set(gmt["HALLMARK_GLYCOLYSIS"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
META = pd.read_csv(HERE / "samples.tsv", sep="\t")
assert len(META) == 168 and (META.pathology == "MVM").sum() == 49


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    out = d / URL.rsplit("/", 1)[1]
    if not out.exists():
        urllib.request.urlretrieve(URL, out)
    print("ok", out.name)


def load(d):
    c = pd.read_csv(d / URL.rsplit("/", 1)[1], index_col=0)
    num = c.select_dtypes(include=[np.number])
    # sample columns must match GSM titles or GSM ids (asserted)
    m = META.copy()
    col = {}
    for _, r in m.iterrows():
        hits = [x for x in num.columns if x == r.title or x == r.gsm]
        assert len(hits) == 1, (r.title, r.gsm, hits[:3])
        col[r.gsm] = hits[0]
    num = num[[col[g] for g in m.gsm]]
    num.columns = m.gsm.tolist()
    if not set(["FLT1"]) & set(num.index):  # Ensembl IDs: map with the repo table
        mp = pd.read_csv(ROOT / "data_external/ensembl_gene_id_to_name_GRCh37.87.tsv", sep="\t", header=None, names=["id", "sym"])
        mp = mp.drop_duplicates("id").set_index("id")["sym"]
        sym = pd.Series(num.index.astype(str).str.split(".").str[0], index=num.index).map(mp)
        num = num[sym.notna().values].groupby(sym.dropna().values).sum()
    num = num.groupby(num.index).sum()
    assert "FLT1" in num.index
    return num, m.set_index("gsm")


def score(z, genes):
    s = z.loc[z.index.intersection(sorted(genes))].mean()
    return (s - s.mean()) / s.std()


def fit(y, m, disease=True, subset=None):
    df = m.assign(y=y.values, MVM=(m.pathology == "MVM").astype(int), GA=m.GA.astype(float))
    if subset is not None:
        df = df[subset]
    f = "y ~ MVM + GA + C(sex) + C(batch)" + (" + C(disease)" if disease and df.disease.nunique() > 1 else "")
    r = smf.ols(f, data=df).fit(cov_type="HC3")
    return float(r.params["MVM"]), float(r.bse["MVM"]), float(1 - stats.norm.cdf(r.params["MVM"] / r.bse["MVM"])), int(df.MVM.sum()), int((1 - df.MVM).sum())


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    counts, m = load(d)
    cpm = counts / counts.sum() * 1e6
    keep = (cpm >= 1).sum(axis=1) >= 0.1 * cpm.shape[1]
    lg = np.log2(cpm[keep] + 1)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
    sc = {"hyp_only": score(z, HYP - GLY), "hyp_all": score(z, HYP), "hif2dep": score(z, HIF2), "FLT1": z.loc["FLT1"]}
    res = {"n_genes": {k: int(len(z.index.intersection(sorted(v)))) for k, v in
                       {"hyp_only": HYP - GLY, "hyp_all": HYP, "hif2dep": HIF2}.items()}}
    est, se, p, n1, n0 = fit(sc["hyp_only"], m)
    # random-set calibration
    mean_expr = lg.mean(axis=1)
    dec = pd.qcut(mean_expr.rank(method="first"), 10, labels=False)
    target = dec[z.index.intersection(sorted(HYP - GLY))].value_counts()
    pool = dec[~dec.index.isin(HYP | HIF2 | {"FLT1"})]
    rng = np.random.default_rng(SEED); null = []
    for _ in range(NNULL):
        genes = np.concatenate([rng.choice(pool.index[pool == k], size=n, replace=False) for k, n in target.items()])
        null.append(fit(score(z, set(genes)), m)[0])
    null = np.array(null)
    res["W1"] = {"MVM_coef_SD": est, "se_HC3": se, "p_one_sided": p, "n_MVM": n1, "n_nonMVM": n0,
                 "random_median": float(np.median(null)), "random_p95": float(np.percentile(null, 95)),
                 "random_percentile": float((null < est).mean() * 100)}
    res["W1"]["supported"] = bool(p < 0.05 and est > res["W1"]["random_p95"])
    hdp = (m.disease != "Normal").values
    res["sens_a_HDP_only"] = dict(zip(("coef", "se", "p_one_sided", "n_MVM", "n_nonMVM"), fit(sc["hyp_only"], m, subset=hdp)))
    res["sens_b_no_disease_covariate"] = dict(zip(("coef", "se", "p_one_sided", "n_MVM", "n_nonMVM"), fit(sc["hyp_only"], m, disease=False)))
    res["sens_c_hyp_all"] = dict(zip(("coef", "se", "p_one_sided", "n_MVM", "n_nonMVM"), fit(sc["hyp_all"], m)))
    for k in ("hif2dep", "FLT1"):
        res[f"W2_{k}"] = dict(zip(("coef", "se", "p_one_sided", "n_MVM", "n_nonMVM"), fit(sc[k], m)))
    mvm = (m.pathology == "MVM").values
    res["W3_coupling"] = {g: float(stats.spearmanr(sc["hyp_all"][sel], sc["FLT1"][sel]).statistic)
                          for g, sel in (("MVM", mvm), ("nonMVM", ~mvm))}
    (HERE / "W_results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
