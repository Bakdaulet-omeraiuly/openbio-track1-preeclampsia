"""Addendum Q (plan in plan.md): HIF1A-KO vs EPAS1-KO in hypoxic human T cells (GSE324038), PT2385-defined program.

Usage: python q_hif1_vs_hif2.py GSE324038_Project_Knm32.rsem_genes_expected_count.txt.gz
Input: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE324nnn/GSE324038/suppl/GSE324038_Project_Knm32.rsem_genes_expected_count.txt.gz
"""

import itertools
import json
from pathlib import Path
import re
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((ROOT / "artifacts/addendum-l/L1_HIF2dep_genes.txt").read_text().split())
emap = pd.read_csv(ROOT / "data_external/ensembl_gene_id_to_name_GRCh37.87.tsv", sep="\t", header=None, index_col=0)[1]

raw = pd.read_csv(sys.argv[1], sep="\t", index_col=0).drop(columns=["transcript_id(s)"])
raw.index = raw.index.str.replace(r"\.\d+$", "", regex=True)
sym = emap.reindex(raw.index)
counts = raw[sym.notna().values].groupby(sym[sym.notna()].values).sum()


def parse(col):
    m = re.match(r"Sample_LRS[_-](\d+)_(HIF1a_KO|HIF2a_KO|NC|NT)_", col)
    return m.group(1), m.group(2)


info = pd.DataFrame([parse(c) for c in counts.columns], index=counts.columns, columns=["donor", "arm"])
keep = info.index[info.arm.isin(["HIF1a_KO", "HIF2a_KO", "NT"])]
info = info.loc[keep]
assert len(keep) == 12 and info.donor.nunique() == 4
c = counts[keep]
cpm = c / c.sum() * 1e6
c = c[(cpm >= 1).sum(axis=1) >= 4]
lc = np.log2(c / c.sum() * 1e6 + 1)
z = lc.sub(lc.mean(axis=1), axis=0).div(lc.std(axis=1), axis=0).T


def sc(genes):
    cols = sorted(set(genes) & set(z.columns))
    return z[cols].mean(axis=1), len(cols)


def deltas(s):
    d1, d2 = {}, {}
    for dn in sorted(info.donor.unique()):
        sel = info[info.donor == dn]
        nt = s[sel.index[sel.arm == "NT"]].iloc[0]
        d1[dn] = float(s[sel.index[sel.arm == "HIF1a_KO"]].iloc[0] - nt)
        d2[dn] = float(s[sel.index[sel.arm == "HIF2a_KO"]].iloc[0] - nt)
    return pd.Series(d1), pd.Series(d2)


res = {"n_genes_after_filter": int(z.shape[1])}
hif, n_hif = sc(HIF2)
hyp, n_hyp = sc(HYP)
res["coverage"] = {"hif2dep_expressed": n_hif, "hif2dep_genes": sorted(set(HIF2) & set(z.columns)), "hyp_all_expressed": n_hyp}
res["gate_pass"] = n_hif >= 30
ko = {}
for g in ["HIF1A", "EPAS1", "FLT1"]:
    if g in z.columns:
        a, b = deltas(z[g]); ko[g] = {"delta_HIF1KO": float(a.mean()), "delta_HIF2KO": float(b.mean())}
    else:
        ko[g] = "not expressed after filter"
res["single_genes"] = ko
h1, h2 = deltas(hyp)
res["positive_control_hyp_all"] = {"delta_HIF1KO": float(h1.mean()), "delta_HIF2KO": float(h2.mean()),
                                   "D_hyp_all": float((h2 - h1).mean()), "pass": bool(h1.mean() < 0)}
a, b = deltas(hif)
D = b - a
signs = [np.mean(D.values * np.array(s)) for s in itertools.product([1, -1], repeat=len(D))]
res["Q1"] = {"delta_HIF1KO": float(a.mean()), "delta_HIF2KO": float(b.mean()), "D": float(D.mean()),
             "per_donor_D": D.to_dict(), "n_donors_D_negative": int((D < 0).sum()),
             "p_signflip_one_sided": float(np.mean(np.array(signs) <= D.mean() + 1e-12))}
mean_expr = lc.mean(axis=1)
pool = sorted(set(z.columns) - HYP - HIF2 - {"FLT1"})
hin = sorted(HIF2 & set(z.columns))
dec = pd.qcut(mean_expr.reindex(sorted(z.columns)), 10, labels=False)
need = dec.reindex(hin).value_counts()
by = {k: np.array([g for g in pool if dec.get(g) == k]) for k in range(10)}
rng = np.random.default_rng(0)
null = []
for _ in range(1000):
    pick = np.concatenate([rng.choice(by[k], n, replace=False) for k, n in need.items()])
    x1, x2 = deltas(z[list(pick)].mean(axis=1)); null.append(float((x2 - x1).mean()))
null = np.array(null)
res["Q1_calibration"] = {"null_median": float(np.median(null)), "null_5th": float(np.percentile(null, 5)),
                         "percentile_of_observed": float((null < D.mean()).mean() * 100)}
res["reading"] = ("uninterpretable (coverage gate failed)" if not res["gate_pass"] else
                  ("expressed part of the program is HIF-2-preferring in human T cells"
                   if D.mean() < 0 and D.mean() < res["Q1_calibration"]["null_5th"] else "not shown"))
(HERE / "Q_results.json").write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps({k: v for k, v in res.items() if k != "coverage"} | {"coverage_n": n_hif}, indent=2))
