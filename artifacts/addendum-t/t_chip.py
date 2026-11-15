"""Addendum T: ARNT ChIP-seq binding at the HIF2dep program, Hallmark hypoxia and FLT1 (rules in plan.md).

  python t_chip.py fetch DIR   # GSE209936 optimal peaks + UCSC hg38 ncbiRefSeq
  python t_chip.py run DIR
"""

import gzip
import json
from pathlib import Path
import sys
import urllib.request

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
URLS = {"peaks.bed.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE209nnn/GSE209936/suppl/GSE209936_output_optimal.bed.gz",
        "ncbiRefSeq.txt.gz": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeq.txt.gz"}
WINDOW, SENS = 10_000, (5_000, 50_000)
NNULL, SEED = 1000, 0


def fetch(d):
    d.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        if not (d / name).exists():
            urllib.request.urlretrieve(url, d / name)
        print("ok", name)


def load(d):
    peaks = pd.read_csv(d / "peaks.bed.gz", sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    ref = pd.read_csv(d / "ncbiRefSeq.txt.gz", sep="\t", header=None, usecols=[1, 2, 3, 4, 5, 12],
                      names=["tx", "chr", "strand", "txStart", "txEnd", "gene"])
    ref = ref[ref.tx.str.startswith("NM_")]
    ref["tss"] = np.where(ref.strand == "+", ref.txStart, ref.txEnd)
    return peaks, ref


def bound_genes(peaks, ref, window):
    out = set()
    by_chr = {c: g.sort_values("start") for c, g in peaks.groupby("chr")}
    for c, g in ref.groupby("chr"):
        if c not in by_chr:
            continue
        p = by_chr[c]; starts, ends = p.start.values, p.end.values
        order = np.argsort(ends); ends_sorted = ends[order]
        for gene, tss in zip(g.gene.values, g.tss.values):
            lo, hi = tss - window, tss + window
            # any peak with start <= hi and end >= lo
            i = np.searchsorted(starts, hi, side="right")
            if i and (ends[:i] >= lo).any():
                out.add(gene)
    return out


def enrichment(target, universe_expr, bound, exclude, rng):
    target = [g for g in target if g in universe_expr.index]
    dec = pd.qcut(universe_expr.rank(method="first"), 10, labels=False)
    obs = float(np.mean([g in bound for g in target]))
    pool = dec[~dec.index.isin(exclude)]
    need = dec[target].value_counts()
    null = []
    for _ in range(NNULL):
        draw = np.concatenate([rng.choice(pool.index[pool == k], size=n, replace=False) for k, n in need.items()])
        null.append(np.mean([g in bound for g in draw]))
    null = np.array(null)
    return {"n_genes": len(target), "fraction_bound": obs, "null_median": float(np.median(null)),
            "null_p95": float(np.percentile(null, 95)), "empirical_p": float((1 + (null >= obs).sum()) / (NNULL + 1)),
            "fold_vs_null_median": float(obs / np.median(null)) if np.median(null) > 0 else None}


def main():
    cmd, d = sys.argv[1], Path(sys.argv[2])
    if cmd == "fetch":
        return fetch(d)
    peaks, ref = load(d)
    expr = pd.read_csv(ROOT / "artifacts/addendum-l/L1_GSE216799_PT2385_vs_DMSO.csv", index_col=0)["baseMean"]
    expr = expr[(expr >= 10) & expr.index.isin(set(ref.gene))]
    gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
    hyp = set(gmt["HALLMARK_HYPOXIA"])
    hif2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
    res = {"n_peaks": int(len(peaks)), "n_universe": int(len(expr))}
    for w in (WINDOW,) + SENS:
        b = bound_genes(peaks, ref, w)
        bu = b & set(expr.index)
        rng = np.random.default_rng(SEED)
        t0 = enrichment(sorted(hyp), expr, bu, hyp | hif2 | {"FLT1"}, rng)
        rng = np.random.default_rng(SEED)
        t1 = enrichment(sorted(hif2), expr, bu, hyp | hif2 | {"FLT1"}, rng)
        res[f"window_{w}"] = {"base_rate_universe": float(len(bu) / len(expr)), "T0_hypoxia": t0, "T1_hif2dep": t1,
                              "FLT1_bound": "FLT1" in b}
    f = ref[ref.gene == "FLT1"]
    pc = peaks[peaks.chr == f.chr.iloc[0]]
    dist = min(int(np.min(np.maximum(0, np.maximum(pc.start.values - t, t - pc.end.values)))) for t in f.tss.unique())
    res["T2_FLT1"] = {"nearest_peak_bp_from_any_TSS": dist, "in_universe": bool("FLT1" in expr.index)}
    m = res[f"window_{WINDOW}"]
    m["T0_passed"] = bool(m["T0_hypoxia"]["empirical_p"] < 0.05)
    m["T1_supported"] = bool(m["T0_passed"] and m["T1_hif2dep"]["empirical_p"] < 0.05)
    (HERE / "T_results.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
