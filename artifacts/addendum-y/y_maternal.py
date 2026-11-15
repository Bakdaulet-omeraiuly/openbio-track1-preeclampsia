"""Addendum Y: placental FLT1 vs maternal disease (rules in plan.md).

Usage: python y_maternal.py EXT_DIR W_DIR
  EXT_DIR: GSE75010_series_matrix.txt.gz + GPL6244.txt (fetch_external_cohorts.sh)
  W_DIR:   GSE303463_raw_counts.csv.gz (Addendum W)
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
sys.path.insert(0, str(ROOT / "artifacts/addendum-w"))
import w_mvm as W  # noqa: E402


def read_gse75010(ext):
    meta, table, inside = {}, [], False
    with gzip.open(ext / "GSE75010_series_matrix.txt.gz", "rt") as h:
        for line in h:
            if line.startswith("!series_matrix_table_begin"):
                inside = True; continue
            if line.startswith("!series_matrix_table_end"):
                break
            v = [x.strip('"') for x in line.rstrip("\n").split("\t")]
            if inside:
                table.append(v)
            elif v[0] == "!Sample_characteristics_ch1":
                meta[v[1].split(":")[0]] = [x.split(": ", 1)[-1] for x in v[1:]]
            elif v[0] == "!Sample_geo_accession":
                meta["gsm"] = v[1:]
    m = pd.DataFrame(meta).set_index("gsm")
    e = pd.DataFrame([r[1:] for r in table[1:]], index=[r[0] for r in table[1:]], columns=table[0][1:]).apply(pd.to_numeric, errors="coerce")
    if e.max().max() > 100:
        e = np.log2(e.clip(lower=0) + 1)
    rows, ins = [], False
    for line in (ext / "GPL6244.txt").open():
        if line.startswith("!platform_table_begin"):
            ins = True; continue
        if line.startswith("!platform_table_end"):
            break
        if ins:
            rows.append(line.rstrip("\n").split("\t"))
    gpl = pd.DataFrame(rows[1:], columns=rows[0][: len(rows[1])]).set_index("ID")["gene_assignment"]
    sym = gpl.map(lambda s: s.split(" // ")[1].strip() if isinstance(s, str) and " // " in s else None).reindex(e.index)
    g = e[sym.notna()].groupby(sym[sym.notna()].values).mean()
    z = g.sub(g.mean(1), axis=0).div(g.std(1), axis=0)
    return m, z.loc["FLT1"]


def partial_spearman(x, y, cov):
    X = np.column_stack([np.ones(len(x)), cov])
    r = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]
    return stats.spearmanr(r(x), r(y))


def main():
    ext, wdir = Path(sys.argv[1]), Path(sys.argv[2])
    res = {}
    m, flt = read_gse75010(ext)
    m = m.assign(FLT1=flt.reindex(m.index).values, GA=m["ga (week)"].astype(float) + m["ga (day)"].astype(float) / 7)
    pe = m[m.diagnosis == "PE"]
    out, ps = {}, []
    for name, col, f in (("max_systolic_BP", "maximum systolic bp", float), ("max_diastolic_BP", "maximum diastolic bp", float),
                         ("proteinuria_mode", "mode proteinuria", float), ("HELLP", "hellp diagnosis", lambda v: float(v == "Yes"))):
        v = pd.to_numeric(pe[col], errors="coerce") if f is float else pe[col].map(f)
        ok = v.notna() & pe.FLT1.notna()
        r = partial_spearman(pe.FLT1[ok].values, v[ok].values.astype(float), pe.GA[ok].values)
        p1 = float(r.pvalue / 2 if r.statistic > 0 else 1 - r.pvalue / 2)
        out[name] = {"partial_rho": float(r.statistic), "n": int(ok.sum()), "p_one_sided": p1}; ps.append(p1)
    for (k, q) in zip(out, multipletests(ps, method="holm")[1]):
        out[k]["holm_p"] = float(q); out[k]["supported"] = bool(q < 0.05 and out[k]["partial_rho"] > 0)
    res["Y1_GSE75010_within_PE"] = out

    counts, wm = W.load(wdir)
    cpm = counts / counts.sum() * 1e6
    lg = np.log2(cpm[(cpm >= 1).sum(axis=1) >= 0.1 * cpm.shape[1]] + 1)
    z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
    dec = pd.qcut(lg.mean(axis=1).rank(method="first"), 10, labels=False)

    def fit(gene, sel, case):
        df = wm[sel].assign(y=z.loc[gene][sel].values, case=case[sel].astype(int), GA=wm.GA[sel].astype(float))
        r = smf.ols("y ~ case + GA + C(sex) + C(batch)", data=df).fit(cov_type="HC3")
        return float(r.params["case"]), float(r.bse["case"]), float(1 - stats.norm.cdf(r.params["case"] / r.bse["case"]))

    pe_mask = wm.disease.isin(["PE", "siPE"]).values
    htn_mask = wm.disease.isin(["gHTN", "cHTN"]).values
    norm_mask = (wm.disease == "Normal").values
    sel2 = pe_mask | htn_mask
    est, se, p = fit("FLT1", sel2, pd.Series(pe_mask, index=wm.index))
    same = dec.index[(dec == dec["FLT1"]) & (dec.index != "FLT1")]
    rng = np.random.default_rng(0)
    null = np.array([fit(g, sel2, pd.Series(pe_mask, index=wm.index))[0] for g in rng.choice(same, size=min(1000, len(same)), replace=False)])
    res["Y2_GSE303463_PE_vs_HTN"] = {"coef_SD": est, "se": se, "p_one_sided": p, "n_PE": int(pe_mask.sum()), "n_HTN": int(htn_mask.sum()),
                                     "random_median": float(np.median(null)), "random_p95": float(np.percentile(null, 95)),
                                     "random_percentile": float((null < est).mean() * 100)}
    res["Y2_GSE303463_PE_vs_HTN"]["supported"] = bool(p < 0.05 and est > np.percentile(null, 95))
    sel3 = norm_mask | htn_mask
    e3, s3, p3 = fit("FLT1", sel3, pd.Series(htn_mask, index=wm.index))
    res["Y3_GSE303463_HTN_vs_Normal_descriptive"] = {"coef_SD": e3, "se": s3, "p_one_sided": p3, "n_HTN": int(htn_mask.sum()), "n_Normal": int(norm_mask.sum())}
    (HERE / "Y_results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
