"""Addendum K, step 2: screen seven pre-listed Hallmark programs for PE vs FGR-only in the two used cohorts.

Processing follows notebooks 01/05/08: DESeq2 VST, per-gene z-score, mean over set genes.
Writes screen ranking and the winner. No GSE75010 data are read here.
"""

import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXT = ROOT / "data_external"
DATA = ROOT / "openbio-datathon-2026/track-01/data"
CANDIDATES = ["HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
              "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
              "HALLMARK_MTORC1_SIGNALING", "HALLMARK_P53_PATHWAY", "HALLMARK_COMPLEMENT"]

gmt = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2:] for l in open(EXT / "h.all.v2024.1.Hs.symbols.gmt")}
HYP = set(gmt["HALLMARK_HYPOXIA"])
SETS = {name: set(gmt[name]) - HYP for name in CANDIDATES}


def vst_of(counts, groups):
    counts = counts.loc[(counts >= 10).sum(axis=1) >= 5]
    md = pd.DataFrame({"g": groups}, index=counts.columns)
    dds = DeseqDataSet(counts=counts.T.astype(int), metadata=md, design="~g",
                       inference=DefaultInference(n_cpus=4), quiet=True)
    dds.fit_size_factors()
    dds.vst(use_design=False)
    return pd.DataFrame(dds.layers["vst_counts"], index=counts.columns, columns=counts.index)


def scores(vst):
    z = (vst - vst.mean()) / vst.std()
    return pd.DataFrame({k: z[sorted(g for g in v if g in z.columns)].mean(axis=1) for k, v in SETS.items()}), \
        {k: sum(g in z.columns for g in v) for k, v in SETS.items()}


# GSE203507 villous (discovery)
counts = pd.read_csv(DATA / "track01_counts_vt.csv", index_col=0)
meta = pd.read_csv(DATA / "track01_metadata_vt.csv", index_col=0).loc[counts.columns]
S_d, n_d = scores(vst_of(counts, meta.pe_status.values))
d = S_d.join(meta[["group", "ga_delivery_weeks", "fetal_sex"]])
d = d[d.group.isin(["EOPE", "EOPE+FGR", "FGR"])].copy()
d["PE"] = (d.group != "FGR").astype(int)

# GSE114691 (replication cohort), same processing as notebook 05/08
tx2gene = pd.read_csv(EXT / "grch37_tx2gene.tsv", sep="\t").set_index("Transcript stable ID")["Gene name"]
parts = []
for g, f in {"control": "ControlONLY", "IUGR": "IUGRONLY", "PE": "PEONLY", "PE+IUGR": "PEIUGR"}.items():
    part = pd.read_csv(EXT / f"GSE114691_MasterCount_{f}.txt.gz", sep="\t", index_col=0)
    part.columns = [f"{g}|{c}" for c in part.columns]
    parts.append(part)
tx = pd.concat(parts, axis=1)
m = tx.index.isin(tx2gene.index)
gc_r = tx[m].groupby(tx2gene.reindex(tx.index[m]).values).sum().round().astype(int)
grp_r = gc_r.columns.str.split("|").str[0]
S_r, n_r = scores(vst_of(gc_r, np.where(grp_r.str.startswith("PE"), "PE", "no_PE")))
r = S_r.copy()
r["group"] = grp_r
r = r[r.group != "control"].copy()
r["PE"] = r.group.str.startswith("PE").astype(int)

rows = []
for name in CANDIDATES:
    fd = smf.ols(f"{name} ~ PE + ga_delivery_weeks + C(fetal_sex)", data=d).fit()
    fr = smf.ols(f"{name} ~ PE", data=r).fit()
    td, tr = float(fd.tvalues["PE"]), float(fr.tvalues["PE"])
    rows.append({"program": name, "genes_GSE203507": n_d[name], "genes_GSE114691": n_r[name],
                 "est_GSE203507": float(fd.params["PE"]), "t_GSE203507": td,
                 "est_GSE114691": float(fr.params["PE"]), "t_GSE114691": tr,
                 "same_sign": bool(np.sign(td) == np.sign(tr)), "min_abs_t": min(abs(td), abs(tr))})
screen = pd.DataFrame(rows).sort_values(["same_sign", "min_abs_t"], ascending=False)
screen.to_csv(HERE / "K_screen_ranking.csv", index=False)
eligible = screen[screen.same_sign]
winner = None if eligible.empty else {
    "program": eligible.iloc[0]["program"],
    "direction": "PE higher" if eligible.iloc[0]["t_GSE114691"] > 0 else "PE lower",
    "sign": int(np.sign(eligible.iloc[0]["t_GSE114691"])),
    "min_abs_t": float(eligible.iloc[0]["min_abs_t"])}
(HERE / "K_screen_winner.json").write_text(json.dumps({"winner": winner}, indent=2) + "\n")
pd.set_option("display.width", 250)
print(screen.round(3).to_string(index=False))
print(json.dumps(winner, indent=2))
