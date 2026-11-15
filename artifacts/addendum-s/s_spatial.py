"""Addendum S: within-donor hypoxia–FLT1 link across Visium spots of intact placenta (GSE319236). Rules in plan.md.

  python s_spatial.py fetch DIR
  python s_spatial.py run DIR
"""

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
SEC = pd.read_csv(HERE / "sections.tsv", sep="\t")
SEC["donor"] = SEC["donor"].replace({"ND": "FND"})
PE = {"FCM", "FAM", "FGS", "FEP", "FVS", "FFM", "FRH", "FSH", "FND", "FMM"}
MIN_UMI, MIN_SPOTS, MIN_DONORS = 500, 200, 8
NPERM, SEED, NNULL, NULL_SEED = 10_000, 42, 1000, 0
MARKERS = {"STB": ["CGA", "CYP19A1", "PSG1", "PSG4", "KISS1", "LGALS13"], "endo": ["PECAM1", "VWF", "CDH5"],
           "stroma": ["COL1A1", "COL3A1", "LUM"], "immune": ["PTPRC", "CD14", "CD68"], "rbc": ["HBB", "HBA1"]}
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
ALL_MARKERS = {g for v in MARKERS.values() for g in v}
assert not (ALL_MARKERS & (HYP | HIF2 | {"FLT1"})), ALL_MARKERS & (HYP | HIF2 | {"FLT1"})
assert "FLT1" not in HYP | HIF2


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    for _, r in SEC.iterrows():
        for c in ("barcodes", "features", "matrix", "positions"):
            url = r[c].replace("ftp://", "https://"); out = d / url.rsplit("/", 1)[1]
            if not out.exists() or out.stat().st_size == 0:
                urllib.request.urlretrieve(url, out)
            print("ok", out.name)


def load(d, r):
    f = lambda c: d / r[c].rsplit("/", 1)[1]
    m = io.mmread(f("matrix")).tocsc()
    feats = pd.read_csv(f("features"), sep="\t", header=None)
    genes = feats[1].astype(str).values
    if feats.shape[1] > 2:
        keep = (feats[2] == "Gene Expression").values; m, genes = m[keep], genes[keep]
    bc = pd.read_csv(f("barcodes"), header=None)[0].astype(str).values
    pos = pd.read_csv(f("positions"), header=None)
    if not str(pos.iloc[0, 1]).lstrip("-").isdigit():  # header row present
        pos = pd.read_csv(f("positions"))
    pos = pos.set_axis(range(pos.shape[1]), axis=1)
    in_tissue = set(pos.loc[pos[1].astype(int) == 1, 0].astype(str))
    total = np.asarray(m.sum(axis=0)).ravel()
    keep = np.array([b in in_tissue for b in bc]) & (total >= MIN_UMI)
    return m[:, np.where(keep)[0]], genes


def donor_tables(d):
    out = {}
    for donor, rows in SEC.groupby("donor"):
        mats, secs, gref = [], [], None
        for i, (_, r) in enumerate(rows.iterrows()):
            m, genes = load(d, r)
            if gref is None:
                gref = genes
            assert (genes == gref).all()
            mats.append(m); secs.append(np.full(m.shape[1], float(i)))
        m = sparse.hstack(mats).tocsr()
        total = np.asarray(m.sum(axis=0)).ravel(); ng = np.asarray((m > 0).sum(axis=0)).ravel()
        mito = np.asarray(m[np.where(pd.Series(gref).str.startswith("MT-"))[0]].sum(axis=0)).ravel() / total
        out[donor] = {"m": m, "genes": gref, "total": total, "ngenes": ng, "mito": mito, "section": np.concatenate(secs)}
    return out


def resid(y, X):
    return y - X @ np.linalg.lstsq(X, y, rcond=None)[0]


def sign_flip_p(r, rng):
    r = np.asarray(r); obs = r.mean(); k = len(r)
    if k <= 16:
        null = np.array([(r * np.array(s)).mean() for s in itertools.product([1, -1], repeat=k)])
        return float((null >= obs - 1e-12).mean())
    null = np.array([(r * rng.choice([1, -1], size=k)).mean() for _ in range(NPERM)])
    return float((1 + (null >= obs - 1e-12).sum()) / (NPERM + 1))


def analyse(tables, donors, composition=True, stb_rich=False):
    genes = tables[donors[0]]["genes"]
    det = sum(np.asarray((tables[x]["m"] > 0).sum(axis=1)).ravel() for x in donors)
    n = sum(tables[x]["m"].shape[1] for x in donors)
    uni = det / n >= 0.05
    gi = {g: i for i, g in enumerate(genes)}
    flt1 = gi["FLT1"]; uni[flt1] = False
    for g in ALL_MARKERS:
        if g in gi:
            uni[gi[g]] = False
    uni_idx = np.where(uni)[0]; ug = genes[uni_idx]
    hyp_cols = np.where(np.isin(ug, sorted(HYP)))[0]; hif_cols = np.where(np.isin(ug, sorted(HIF2)))[0]
    per, pooled_mean = {}, np.zeros(len(uni_idx))
    for x in donors:
        t = tables[x]
        lg = (sparse.diags(1e4 / t["total"]) @ t["m"].T).tocsc().log1p()
        Z = lg[:, uni_idx].toarray(); pooled_mean += Z.sum(axis=0)
        sd = Z.std(0); Z = (Z - Z.mean(0)) / np.where(sd > 0, sd, np.nan)
        flt = np.asarray(lg[:, flt1].todense()).ravel()
        comp = []
        for name, ms in MARKERS.items():
            idx = [gi[g] for g in ms if g in gi]
            M = lg[:, idx].toarray(); s_ = M.std(0)
            comp.append(np.nanmean((M - M.mean(0)) / np.where(s_ > 0, s_, np.nan), axis=1))
        comp = np.nan_to_num(np.column_stack(comp))
        cols = [np.ones(len(flt)), np.log10(t["total"]), np.log10(t["ngenes"]), t["mito"]]
        if np.ptp(t["section"]) > 0:
            cols.append(t["section"])
        if composition:
            cols += list(comp.T)
        keep = np.ones(len(flt), bool)
        if stb_rich:
            keep = comp[:, 0] >= np.median(comp[:, 0])
        cov = np.column_stack(cols)[keep]
        per[x] = {"Z": Z[keep], "flt_r": resid(flt[keep], cov), "cov": cov, "n": int(keep.sum())}

    def rhos(cols, extra=None):
        out = []
        for x in donors:
            p = per[x]; s = np.nanmean(p["Z"][:, cols], axis=1)
            cov = p["cov"] if extra is None else np.column_stack([p["cov"], np.nanmean(p["Z"][:, extra], axis=1)])
            fr = p["flt_r"] if extra is None else resid(p["flt_r"], cov)
            out.append(float(stats.spearmanr(resid(s, cov), fr).statistic))
        return out

    dec = pd.qcut(pd.Series(pooled_mean).rank(method="first"), 10, labels=False).values
    ok = ~np.isin(ug, sorted(HYP | HIF2))

    def null(target, extra=None):
        rng = np.random.default_rng(NULL_SEED); need = pd.Series(dec[target]).value_counts(); vals = []
        for _ in range(NNULL):
            cols = np.concatenate([rng.choice(np.where((dec == d_) & ok)[0], size=k, replace=False) for d_, k in need.items()])
            vals.append(float(np.mean(rhos(cols, extra))))
        return np.array(vals)

    rng = np.random.default_rng(SEED)
    r1, r2 = rhos(hyp_cols), rhos(hif_cols, hyp_cols)
    n1, n2 = null(hyp_cols), null(hif_cols, hyp_cols)
    summ = lambda r, nl: {"per_donor_rho": dict(zip(donors, r)), "mean_rho": float(np.mean(r)),
                          "n_positive": int(sum(v > 0 for v in r)), "sign_flip_p": sign_flip_p(r, rng),
                          "random_null_median": float(np.median(nl)), "random_null_p95": float(np.percentile(nl, 95)),
                          "random_percentile": float((nl < np.mean(r)).mean() * 100)}
    s1, s2 = summ(r1, n1), summ(r2, n2)
    s1["supported"] = bool(s1["sign_flip_p"] < 0.05 and s1["mean_rho"] > s1["random_null_p95"])
    return {"donors": donors, "n_spots": {x: per[x]["n"] for x in donors}, "n_universe": int(len(uni_idx)),
            "n_hyp": int(len(hyp_cols)), "n_hif2dep": int(len(hif_cols)), "S1": s1, "S2_descriptive": s2}


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    tables = donor_tables(d)
    spots = {x: int(t["m"].shape[1]) for x, t in tables.items()}
    eligible = sorted(x for x, v in spots.items() if v >= MIN_SPOTS)
    res = {"gate": {"spots_per_donor": spots, "eligible": eligible, "passed": len(eligible) >= MIN_DONORS}}
    if res["gate"]["passed"]:
        res["primary"] = analyse(tables, eligible)
        sens = {}
        if "FND" in eligible:
            sens["a_without_ND"] = analyse(tables, [x for x in eligible if x != "FND"])["S1"]
        sens["b_no_composition"] = analyse(tables, eligible, composition=False)["S1"]
        sens["c_STB_rich_half"] = analyse(tables, eligible, stb_rich=True)["S1"]
        res["sensitivity"] = sens
    (HERE / "S_results.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
    p = res.get("primary", {})
    print(json.dumps({"gate": res["gate"], "S1": p.get("S1"), "S2": p.get("S2_descriptive")}, indent=1, default=float))
    for k, v in res.get("sensitivity", {}).items():
        print(k, "mean %.3f pos %d/%d p %.4f p95 %.3f supported %s" % (v["mean_rho"], v["n_positive"], len(v["per_donor_rho"]),
                                                                       v["sign_flip_p"], v["random_null_p95"], v["supported"]))


if __name__ == "__main__":
    main()
