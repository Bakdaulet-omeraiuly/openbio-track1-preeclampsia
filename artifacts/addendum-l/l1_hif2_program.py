"""Addendum L step 1: HIF2dep program from GSE216799 (PT2385 vs DMSO, primary CTB at 2.5% O2, 3 placentas).

Input: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE216nnn/GSE216799/suppl/GSE216799_raw_counts.tsv.gz
Usage: python l1_hif2_program.py GSE216799_raw_counts.tsv.gz GSE216799_series_matrix.txt.gz
"""

import gzip
import json
from pathlib import Path
import sys
import warnings

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

counts = pd.read_csv(sys.argv[1], sep="\t", index_col=0)
meta = {}
with gzip.open(sys.argv[2], "rt") as handle:
    for line in handle:
        values = [v.strip('"') for v in line.rstrip("\n").split("\t")]
        if values[0] in ("!Sample_geo_accession", "!Sample_description", "!Sample_title"):
            meta[values[0]] = values[1:]
        if values[0] == "!Sample_characteristics_ch1" and values[1].startswith("treatment"):
            meta["treatment"] = [v.split(": ", 1)[1] for v in values[1:]]
md = pd.DataFrame({"gsm": meta["!Sample_geo_accession"], "desc": meta["!Sample_description"],
                   "treatment": meta["treatment"]})
md["donor"] = md["desc"].str.split("-").str[0]
md["treatment"] = np.where(md["treatment"].str.startswith("PT2385"), "PT2385", "DMSO")
print("count columns:", list(counts.columns))
print(md.to_string(index=False))
# match count columns to samples by GSM, description or title
column_map = {}
for col in counts.columns:
    hit = md[(md.gsm == col) | (md.desc == col) | (md.desc.str.replace("-", "_") == col)]
    assert len(hit) == 1, f"cannot map count column {col}"
    column_map[col] = hit.index[0]
counts = counts.rename(columns=lambda c: md.loc[column_map[c], "desc"])
md = md.set_index("desc").loc[counts.columns]
assert sorted(md.groupby("donor").size().tolist()) == [2, 2, 2]

# gene symbols: counts index may be Ensembl or symbols
index_is_ensembl = counts.index.astype(str).str.startswith("ENSG").mean() > 0.5
if index_is_ensembl:
    # Ensembl gene id -> gene name from the Ensembl GRCh37.87 GTF (data_external/ensembl_gene_id_to_name_GRCh37.87.tsv);
    # ids without a name in that release are dropped; counts of ids sharing a name are summed.
    names = pd.read_csv(ROOT / "data_external/ensembl_gene_id_to_name_GRCh37.87.tsv", sep="\t", header=None,
                        index_col=0).iloc[:, 0]
    mapped = counts.index.str.split(".").str[0].map(names)
    print("Ensembl ids mapped to names:", int(mapped.notna().sum()), "of", len(mapped))
    counts = counts[mapped.notna()]
    counts.index = mapped[mapped.notna()]
counts = counts.groupby(counts.index).sum()
keep = (counts >= 10).sum(axis=1) >= 3
counts = counts[keep].round().astype(int)

dds = DeseqDataSet(counts=counts.T, metadata=md[["donor", "treatment"]], design="~ donor + treatment",
                   inference=DefaultInference(n_cpus=4), quiet=True)
dds.deseq2()
ds = DeseqStats(dds, contrast=["treatment", "PT2385", "DMSO"], inference=DefaultInference(n_cpus=4), quiet=True)
ds.summary()
res = ds.results_df.copy()
res.to_csv(HERE / "L1_GSE216799_PT2385_vs_DMSO.csv")

gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(ROOT / "data_external/h.all.v2024.1.Hs.symbols.gmt")}
excluded = set(gmt["HALLMARK_HYPOXIA"]) | {"FLT1"}
ranked = res[~res.index.isin(excluded)].dropna(subset=["stat"]).sort_values("stat")
hif2dep = ranked.index[:100].tolist()
(HERE / "L1_HIF2dep_genes.txt").write_text("\n".join(hif2dep) + "\n")

f1 = pd.read_csv(ROOT / "notebooks/results/F1_module_genes.csv")["gene"].tolist()
vst = np.log2(counts / counts.sum() * 1e6 + 1).T
z = (vst - vst.mean()) / vst.std()
f1_present = [g for g in f1 if g in z.columns]
f1_score = z[f1_present].mean(axis=1)
paired = {d: float(f1_score[(md.donor == d) & (md.treatment == "PT2385")].iloc[0]
                   - f1_score[(md.donor == d) & (md.treatment == "DMSO")].iloc[0]) for d in sorted(md.donor.unique())}
f1_res = res.reindex(f1_present).dropna(subset=["log2FoldChange"])
summary = {
    "genes_tested": int(len(res)),
    "FLT1": {k: (None if pd.isna(v) else float(v)) for k, v in res.loc["FLT1", ["log2FoldChange", "stat", "pvalue", "padj"]].items()}
    if "FLT1" in res.index else None,
    "HIF2dep_size": len(hif2dep),
    "HIF2dep_top10": hif2dep[:10],
    "HIF2dep_max_stat": float(ranked["stat"].iloc[99]),
    "F1_genes_present": len(f1_present),
    "F1_score_paired_difference_PT2385_minus_DMSO": paired,
    "F1_fraction_negative_LFC": float((f1_res["log2FoldChange"] < 0).mean()),
    "all_genes_fraction_negative_LFC": float((res["log2FoldChange"].dropna() < 0).mean()),
    "note": "n = 3 placentas; descriptive only",
}
(HERE / "L1_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary, indent=2))
