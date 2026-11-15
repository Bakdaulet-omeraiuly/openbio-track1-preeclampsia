"""Addendum R: within-donor, within-STB hypoxia–FLT1 link in GSE319236 single-cell data (rules in plan.md).

Download (per-library CellRanger filtered matrices, URLs in libraries.tsv):
  python r_within_donor.py fetch DIR
Run once:
  python r_within_donor.py run DIR
"""

import gzip
import itertools
import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd
from scipy import io, sparse, stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
LIBS = pd.read_csv(HERE / "libraries.tsv", sep="\t")
PE = {"FCM", "FAM", "FGS", "FEP", "FVS", "FFM", "FRH", "FSH", "FND", "FMM"}  # FND taken as "FD" (plan.md)
CONTROL = {"FVQ", "FLR", "FJJ", "FVB", "FLJ", "FJD", "FVT", "FYI", "FRK", "FNS"}
STB_MARKERS = ["CGA", "CYP19A1", "PSG1", "PSG4", "KISS1", "LGALS13"]
MIN_GENES, MAX_MITO, MIN_STB_CELLS, MIN_DONORS = 500, 0.20, 30, 8
NPERM, SEED, NNULL, NULL_SEED = 10_000, 42, 1000, 0

gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
assert "FLT1" not in HYP and "FLT1" not in HIF2
assert set(LIBS.donor) == PE | CONTROL


def fetch(d: Path):
    d.mkdir(parents=True, exist_ok=True)
    for _, r in LIBS.iterrows():
        for col in ("barcodes", "features", "matrix"):
            url = r[col].replace("ftp://", "https://")
            out = d / url.rsplit("/", 1)[1]
            if not out.exists() or out.stat().st_size == 0:
                urllib.request.urlretrieve(url, out)
            print("ok", out.name)


def load_library(d: Path, r):
    name = lambda col: d / r[col].rsplit("/", 1)[1]
    m = io.mmread(name("matrix")).tocsc()  # genes x cells
    feats = pd.read_csv(name("features"), sep="\t", header=None)
    genes = feats[1].astype(str).values
    if feats.shape[1] > 2:  # keep gene expression rows only
        keep = (feats[2] == "Gene Expression").values
        m, genes = m[keep], genes[keep]
    return m, genes


def stb_cells(m, genes, min_markers=3):
    """QC, then STB = at least min_markers of the 6 markers detected. Uses marker presence only."""
    total = np.asarray(m.sum(axis=0)).ravel()
    ngenes = np.asarray((m > 0).sum(axis=0)).ravel()
    mito = np.isin(np.arange(len(genes)), np.where(pd.Series(genes).str.startswith("MT-"))[0])
    mito_frac = np.asarray(m[mito].sum(axis=0)).ravel() / np.maximum(total, 1)
    qc = (ngenes >= MIN_GENES) & (mito_frac < MAX_MITO)
    idx = [np.where(genes == g)[0] for g in STB_MARKERS]
    present = sum((np.asarray(m[i[0]].todense()).ravel() > 0).astype(int) if len(i) else 0 for i in idx)
    sel = qc & (present >= min_markers)
    return sel, total, ngenes, mito_frac


def donor_tables(d: Path, min_markers=3):
    out = {}
    for donor, rows in LIBS.groupby("donor"):
        mats, tot, ng, mf, gene_ref = [], [], [], [], None
        for _, r in rows.iterrows():
            m, genes = load_library(d, r)
            if gene_ref is None:
                gene_ref = genes
            assert (genes == gene_ref).all()
            sel, total, ngenes, mito_frac = stb_cells(m, genes, min_markers)
            mats.append(m[:, np.where(sel)[0]]); tot.append(total[sel]); ng.append(ngenes[sel]); mf.append(mito_frac[sel])
        out[donor] = {"m": sparse.hstack(mats).tocsr(), "genes": gene_ref, "total": np.concatenate(tot),
                      "ngenes": np.concatenate(ng), "mito": np.concatenate(mf)}
    return out


def resid(y, X):
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def sign_flip_p(rhos, rng):
    rhos = np.asarray(rhos); obs = rhos.mean(); k = len(rhos)
    if k <= 16:
        null = np.array([(rhos * np.array(s)).mean() for s in itertools.product([1, -1], repeat=k)])
        return float((null >= obs - 1e-12).mean())
    null = np.array([(rhos * rng.choice([1, -1], size=k)).mean() for _ in range(NPERM)])
    return float((1 + (null >= obs - 1e-12).sum()) / (NPERM + 1))


def analyse(tables, donors, label):
    # gene universe: detected in >= 5% of pooled STB cells of the eligible donors
    genes = tables[donors[0]]["genes"]
    pooled_det = sum(np.asarray((tables[x]["m"] > 0).sum(axis=1)).ravel() for x in donors)
    n_cells = sum(tables[x]["m"].shape[1] for x in donors)
    uni = pooled_det / n_cells >= 0.05
    flt1 = int(np.where(genes == "FLT1")[0][0])
    uni[flt1] = False
    uni_idx = np.where(uni)[0]
    uni_genes = genes[uni_idx]
    hyp_cols = np.where(np.isin(uni_genes, sorted(HYP)))[0]
    hif_cols = np.where(np.isin(uni_genes, sorted(HIF2)))[0]
    per = {}
    for x in donors:
        t = tables[x]
        m = t["m"]
        norm = sparse.diags(1e4 / np.maximum(t["total"], 1)) @ m.T  # cells x genes
        lg = norm.tocsc()
        flt = np.log1p(np.asarray(lg[:, flt1].todense()).ravel())
        Z = np.log1p(lg[:, uni_idx].toarray())
        sd = Z.std(axis=0)
        Z = (Z - Z.mean(axis=0)) / np.where(sd > 0, sd, np.nan)
        cov = np.column_stack([np.ones(len(flt)), np.log10(t["total"]), np.log10(t["ngenes"]), t["mito"]])
        per[x] = {"Z": Z, "flt_r": resid(flt, cov), "cov": cov, "n": int(len(flt))}

    def mean_rho(cols, extra=None):
        rs = []
        for x in donors:
            p = per[x]
            s = np.nanmean(p["Z"][:, cols], axis=1)
            cov = p["cov"] if extra is None else np.column_stack([p["cov"], np.nanmean(p["Z"][:, extra], axis=1)])
            fr = p["flt_r"] if extra is None else resid(p["flt_r"] + 0, cov)
            rs.append(float(stats.spearmanr(resid(s, cov), fr).statistic))
        return rs

    rng = np.random.default_rng(SEED)
    r1 = mean_rho(hyp_cols)
    r2 = mean_rho(hif_cols, extra=hyp_cols)
    # random-set null, decile-matched on pooled mean log-normalised expression
    pooled_mean = np.zeros(len(uni_idx))
    for x in donors:
        t = tables[x]
        pooled_mean += np.asarray((sparse.diags(1e4 / np.maximum(t["total"], 1)) @ t["m"].T).tocsc()[:, uni_idx].log1p().sum(axis=0)).ravel()
    dec = pd.qcut(pd.Series(pooled_mean).rank(method="first"), 10, labels=False).values
    pool_ok = ~np.isin(uni_genes, sorted(HYP | HIF2))

    def null_draws(target_cols, extra=None):
        nrng = np.random.default_rng(NULL_SEED)
        need = pd.Series(dec[target_cols]).value_counts()
        out = []
        for _ in range(NNULL):
            cols = np.concatenate([nrng.choice(np.where((dec == d_) & pool_ok)[0], size=k, replace=False)
                                   for d_, k in need.items()])
            out.append(float(np.mean(mean_rho(cols, extra))))
        return np.array(out)

    null1 = null_draws(hyp_cols)
    null2 = null_draws(hif_cols, extra=hyp_cols)
    res = {
        "label": label, "donors": donors, "n_stb_cells": {x: per[x]["n"] for x in donors},
        "n_universe_genes": int(len(uni_idx)), "n_hyp_genes": int(len(hyp_cols)), "n_hif2dep_genes": int(len(hif_cols)),
        "R1": {"per_donor_rho": dict(zip(donors, r1)), "mean_rho": float(np.mean(r1)),
               "n_positive": int(sum(r > 0 for r in r1)), "sign_flip_p": sign_flip_p(r1, rng),
               "random_null_median": float(np.median(null1)), "random_null_p95": float(np.percentile(null1, 95)),
               "random_percentile": float((null1 < np.mean(r1)).mean() * 100)},
        "R2_descriptive": {"per_donor_rho": dict(zip(donors, r2)), "mean_rho": float(np.mean(r2)),
                           "n_positive": int(sum(r > 0 for r in r2)), "sign_flip_p": sign_flip_p(r2, rng),
                           "random_null_median": float(np.median(null2)), "random_null_p95": float(np.percentile(null2, 95)),
                           "random_percentile": float((null2 < np.mean(r2)).mean() * 100)},
        "donor_means": {x: {"group": "PE" if x in PE else "control",
                            "hyp_all_mean_lognorm": float(np.nanmean(np.log1p((sparse.diags(1e4 / np.maximum(tables[x]["total"], 1)) @ tables[x]["m"].T).tocsc()[:, uni_idx[hyp_cols]].toarray()))),
                            "FLT1_mean_lognorm": float(np.log1p((sparse.diags(1e4 / np.maximum(tables[x]["total"], 1)) @ tables[x]["m"].T).tocsc()[:, flt1].toarray()).mean())}
                        for x in donors},
    }
    res["R1"]["supported"] = bool(res["R1"]["sign_flip_p"] < 0.05 and np.mean(r1) > res["R1"]["random_null_p95"])
    return res


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    tables = donor_tables(d)
    counts = {x: int(t["m"].shape[1]) for x, t in tables.items()}
    eligible = sorted(x for x, n in counts.items() if n >= MIN_STB_CELLS)
    gate = {"stb_cells_per_donor": counts, "eligible_donors": eligible, "passed": len(eligible) >= MIN_DONORS}
    result = {"gate": gate}
    if gate["passed"]:
        main_res = analyse(tables, eligible, "primary")
        dm = main_res["donor_means"]
        for key in ("hyp_all_mean_lognorm", "FLT1_mean_lognorm"):
            a = [v[key] for v in dm.values() if v["group"] == "PE"]; b = [v[key] for v in dm.values() if v["group"] == "control"]
            main_res[f"R3_{key}_PE_vs_control"] = {"median_PE": float(np.median(a)), "median_control": float(np.median(b)),
                                                   "MWU_p_two_sided": float(stats.mannwhitneyu(a, b).pvalue) if a and b else None}
        result["primary"] = main_res
        sens = {}
        if "FND" in eligible:
            sens["a_without_FND"] = analyse(tables, [x for x in eligible if x != "FND"], "without FND")["R1"]
        for g, members in (("PE", PE), ("control", CONTROL)):
            sub = [x for x in eligible if x in members]
            if len(sub) >= 3:
                sens[f"b_{g}_only"] = analyse(tables, sub, f"{g} only")["R1"]
        t4 = donor_tables(d, min_markers=4)
        el4 = sorted(x for x, t in t4.items() if t["m"].shape[1] >= MIN_STB_CELLS)
        if len(el4) >= 3:
            sens["c_markers_ge4"] = analyse(t4, el4, ">=4 markers")["R1"] | {"n_donors": len(el4)}
        result["sensitivity"] = sens
    (HERE / "R_results.json").write_text(json.dumps(result, indent=2, default=float) + "\n")
    print(json.dumps({"gate": gate, "R1": result.get("primary", {}).get("R1"),
                      "R2": result.get("primary", {}).get("R2_descriptive")}, indent=2, default=float))


if __name__ == "__main__":
    main()
