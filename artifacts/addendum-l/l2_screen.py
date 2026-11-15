"""Addendum L step 2: screen H1 (HIF2dep-FLT1 partial Spearman | hyp_all) in the three used cohorts.

Usage: python l2_screen.py GSE75010_series_matrix.txt GPL6244.txt
"""

import json
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from scipy import stats
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXT = ROOT / "data_external"
DATA = ROOT / "openbio-datathon-2026/track-01/data"
gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(EXT / "h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
HIF2 = set((HERE / "L1_HIF2dep_genes.txt").read_text().split())


def vst_of(counts, groups):
    counts = counts.loc[(counts >= 10).sum(axis=1) >= 5]
    md = pd.DataFrame({"g": groups}, index=counts.columns)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=md, design="~g",
                       inference=DefaultInference(n_cpus=4), quiet=True)
    dds.fit_size_factors()
    dds.vst(use_design=False)
    return pd.DataFrame(dds.layers["vst_counts"], index=counts.columns, columns=counts.index)


def score_frame(z):
    return pd.DataFrame({"hyp_all": z[sorted(g for g in HYP if g in z.columns)].mean(axis=1),
                         "hif2dep": z[sorted(g for g in HIF2 if g in z.columns)].mean(axis=1),
                         "FLT1": z["FLT1"]}), int(sum(g in z.columns for g in HIF2))


def partial_spearman(frame):
    X = np.column_stack([np.ones(len(frame)), frame["hyp_all"].values])
    resid = lambda y: y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
    rho, p = stats.spearmanr(resid(frame["hif2dep"].values), resid(frame["FLT1"].values))
    return float(rho), float(p)


rows = []
# GSE203507 villous: disease = EOPE, EOPE+FGR, FGR
counts = pd.read_csv(DATA / "track01_counts_vt.csv", index_col=0)
meta = pd.read_csv(DATA / "track01_metadata_vt.csv", index_col=0).loc[counts.columns]
vst = vst_of(counts, meta.pe_status.values)
S, n = score_frame((vst - vst.mean()) / vst.std())
S = S[meta.group.isin(["EOPE", "EOPE+FGR", "FGR"]).values]
rows.append({"cohort": "GSE203507", "n": len(S), "hif2dep_genes": n, **dict(zip(["rho", "p"], partial_spearman(S)))})

# GSE114691: disease = PE, PE+IUGR, IUGR
tx2gene = pd.read_csv(EXT / "grch37_tx2gene.tsv", sep="\t").set_index("Transcript stable ID")["Gene name"]
parts = []
for g, f in {"control": "ControlONLY", "IUGR": "IUGRONLY", "PE": "PEONLY", "PE+IUGR": "PEIUGR"}.items():
    part = pd.read_csv(EXT / f"GSE114691_MasterCount_{f}.txt.gz", sep="\t", index_col=0)
    part.columns = [f"{g}|{c}" for c in part.columns]
    parts.append(part)
tx = pd.concat(parts, axis=1)
m = tx.index.isin(tx2gene.index)
gc = tx[m].groupby(tx2gene.reindex(tx.index[m]).values).sum().round().astype(int)
grp = gc.columns.str.split("|").str[0]
vst = vst_of(gc, np.where(grp.str.startswith("PE"), "PE", "no_PE"))
S, n = score_frame((vst - vst.mean()) / vst.std())
S = S[(grp != "control")]
rows.append({"cohort": "GSE114691", "n": len(S), "hif2dep_genes": n, **dict(zip(["rho", "p"], partial_spearman(S)))})

# GSE75010 (Addendum K processing): disease = all PE + non-PE IUGR
src = (ROOT / "artifacts/addendum-k/k_confirm_gse75010.py").read_text()
prefix = src[: src.index("gmt = {")]
ns = {"__file__": str(ROOT / "artifacts/addendum-k/k_confirm_gse75010.py"), "__name__": "k"}
sys.argv = ["k", sys.argv[1], sys.argv[2]]
exec(compile(prefix, "k_confirm_prefix", "exec"), ns)
z, meta75 = ns["z"], ns["meta"]
S, n = score_frame(z)
disease = (meta75["diagnosis"] == "PE") | (meta75["iugr diagnosis"] == "Yes")
S = S[disease.values]
rows.append({"cohort": "GSE75010", "n": len(S), "hif2dep_genes": n, **dict(zip(["rho", "p"], partial_spearman(S)))})

out = pd.DataFrame(rows)
zf = np.arctanh(out.rho)
w = out.n - 3
pooled_z = float((w * zf).sum() / w.sum())
se = float(1 / np.sqrt(w.sum()))
p_one = float(1 - stats.norm.cdf(pooled_z / se))
result = {"per_cohort": rows, "pooled_rho": float(np.tanh(pooled_z)), "pooled_one_sided_p": p_one,
          "all_positive": bool((out.rho > 0).all())}
result["screen_pass"] = bool(result["all_positive"] and p_one < 0.05)
(HERE / "L2_screen.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
