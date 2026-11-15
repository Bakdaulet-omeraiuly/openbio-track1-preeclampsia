"""Addendum J2: GSE203507 villous FLT1 isoform gate and pre-specified disease-vs-iPTB contrast.

Written before any GSE203507 counts existed.
  python j2_gse203507.py --work WORK --gate       # label-free gate (labels only for criterion 5) + join check
  python j2_gse203507.py --work WORK --contrast   # refuses unless the gate and the join check passed
"""

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, spearmanr

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SFLT = "ENST00000541932_unique_junction_1"
CAN1 = "ENST00000282397_unique_junction_1"
CAN2 = "ENST00000282397_unique_junction_2"
SFLT_BODY = "ENST00000541932_terminal_body_1"
CAN_BODY = "ENST00000282397_terminal_body_1"
SFLT_BODY_EFF = 567 - 49
CAN_BODY_EFF = 2782 - 49
EXPLORATORY = {
    "ENST00000539099": ["ENST00000539099_unique_junction_1"],
    "ENST00000540678": ["ENST00000540678_unique_junction_1", "ENST00000540678_unique_junction_2"],
    "ENST00000543394": ["ENST00000543394_unique_junction_1"],
}
DISEASE = {"EOPE", "EOPE+FGR", "FGR"}

parser = argparse.ArgumentParser()
parser.add_argument("--work", required=True, type=Path)
parser.add_argument("--gate", action="store_true")
parser.add_argument("--contrast", action="store_true")
args = parser.parse_args()

manifest = pd.read_csv(HERE / "gse203507-vt-manifest.csv")
rows = []
for run in manifest["run"]:
    record = json.loads((args.work / "counts" / f"{run}.json").read_text())
    expected = int(manifest.loc[manifest["run"] == run, "read_count"].iloc[0])
    assert record["total_reads"] == expected, run
    row = {"run": run, "total_reads": record["total_reads"], "candidate_reads": record["candidate_reads"],
           "unique_reads_in_flt1_span": record["alignment_stats"]["unique_reads_in_flt1_span"],
           "flt1_span_forward_strand": record["alignment_stats"]["flt1_span_forward_strand"],
           "flt1_span_reverse_strand": record["alignment_stats"]["flt1_span_reverse_strand"]}
    row.update(record["hisat2_unique_counts_flank10"])
    rows.append(row)
counts = pd.DataFrame(rows)
assert len(counts) == 30
counts["junction_ratio"] = np.log2((counts[SFLT] + 0.5) / ((counts[CAN1] + counts[CAN2]) / 2 + 0.5))
counts["body_ratio"] = np.log2((counts[SFLT_BODY] / SFLT_BODY_EFF + 0.5 / SFLT_BODY_EFF)
                               / (counts[CAN_BODY] / CAN_BODY_EFF + 0.5 / CAN_BODY_EFF))

if args.gate:
    counts.drop(columns=["junction_ratio", "body_ratio"]).to_csv(HERE / "gse203507-marker-counts.csv", index=False)
    each = (counts[SFLT] >= 5) & (counts[CAN1] >= 5) & (counts[CAN2] >= 5)
    depth = (counts[SFLT] + counts[CAN1] + counts[CAN2]) >= 20
    usable = each & depth
    rho = float(spearmanr(counts["junction_ratio"], counts["body_ratio"]).statistic)
    labels = counts[["run"]].merge(manifest[["run", "group"]], on="run")["group"]
    arm = labels.map(lambda g: "disease" if g in DISEASE else g)
    per_group = {a: {"n": int((arm == a).sum()), "usable": int(usable[arm == a].sum()),
                     "fraction": round(float(usable[arm == a].mean()), 4)} for a in ["disease", "iPTB"]}

    processed = pd.read_csv(ROOT / "openbio-datathon-2026/track-01/data/track01_counts_vt.csv", index_col=0)
    cpm = processed.loc["FLT1"] / processed.sum() * 1e6
    sample = counts[["run"]].merge(manifest[["run", "sample_id"]], on="run")["sample_id"]
    join_rho = float(spearmanr(counts["unique_reads_in_flt1_span"] / counts["total_reads"],
                               cpm.reindex(sample).values).statistic)

    result = {
        "population": "30 GSE203507 villous samples (labels used only for criterion 5)",
        "criterion_2_each_primary_junction_ge5": {"fraction": round(float(each.mean()), 4), "pass": bool(each.mean() >= 0.70)},
        "criterion_3_depth_ge20": {"fraction": round(float(depth.mean()), 4), "pass": bool(depth.mean() >= 0.70)},
        "criterion_4_body_vs_junction_rho": {"rho": round(rho, 4), "pass": bool(rho >= 0.70)},
        "criterion_5_usable_per_group": {"groups": per_group,
                                         "pass": bool(all(g["fraction"] >= 0.60 for g in per_group.values()))},
        "sample_join_check": {"spearman_rho_raw_flt1_fraction_vs_processed_flt1_cpm": round(join_rho, 4),
                              "pass": bool(join_rho >= 0.90)},
        "coverage_median": {m: float(counts[m].median()) for m in [SFLT, CAN1, CAN2, SFLT_BODY, CAN_BODY]},
        "coverage_min": {m: float(counts[m].min()) for m in [SFLT, CAN1, CAN2]},
        "strand_reverse_fraction": round(float(counts["flt1_span_reverse_strand"].sum()
                                               / (counts["flt1_span_forward_strand"].sum()
                                                  + counts["flt1_span_reverse_strand"].sum())), 4),
    }
    result["gate_pass"] = all(v["pass"] for k, v in result.items() if k.startswith("criterion_"))
    (HERE / "gse203507-coverage-gate.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))

if args.contrast:
    gate = json.loads((HERE / "gse203507-coverage-gate.json").read_text())
    if not (gate["gate_pass"] and gate["sample_join_check"]["pass"]):
        sys.exit("J2 gate or sample-join check failed: contrast not run (pre-specified).")
    meta = pd.read_csv(ROOT / "openbio-datathon-2026/track-01/data/track01_metadata_vt.csv")
    processed = pd.read_csv(ROOT / "openbio-datathon-2026/track-01/data/track01_counts_vt.csv", index_col=0)
    scores = pd.read_csv(ROOT / "notebooks/results/pathway_scores.csv", index_col=0)
    data = counts.merge(manifest[["run", "sample_id"]], on="run").merge(
        meta[["sample_id", "group", "ga_delivery_weeks", "fetal_sex"]], on="sample_id", validate="one_to_one")
    data["log_flt1"] = np.log2((processed.loc["FLT1"] / processed.sum() * 1e6).reindex(data["sample_id"]).values + 1)
    data["hypoxia"] = scores["HALLMARK_HYPOXIA"].reindex(data["sample_id"]).values
    canonical = (data[CAN1] + data[CAN2]) / 2
    for t, markers in EXPLORATORY.items():
        data[f"expl_{t}"] = np.log2((data[markers].mean(axis=1) + 0.5) / (canonical + 0.5))
    data["disease"] = data["group"].isin(DISEASE).astype(int)
    data["PE"] = data["group"].isin({"EOPE", "EOPE+FGR"}).astype(int)
    base = "ga_delivery_weeks + C(fetal_sex)"

    def fit(frame, term, outcome="junction_ratio", extra="", cov="nonrobust"):
        model = smf.ols(f"{outcome} ~ {term} + {base}{extra}", data=frame).fit(cov_type=cov)
        low, high = model.conf_int().loc[term]
        return {"n": int(model.nobs), "n_term1": int(frame[term].sum()), "coef_log2": float(model.params[term]),
                "CI_low": float(low), "CI_high": float(high), "p": float(model.pvalues[term])}

    primary_frame = data[data["group"] != "TB"].copy()
    rows = [{"analysis": "PRIMARY: disease (EOPE, EOPE+FGR, FGR) vs iPTB, +GA +sex, OLS",
             **fit(primary_frame, "disease")}]
    rows.append({"analysis": "Sensitivity (a): primary model, HC3", **fit(primary_frame, "disease", cov="HC3")})
    observed = rows[0]["coef_log2"]
    rng = np.random.default_rng(0)
    perm = primary_frame.copy()
    hits = 0
    for _ in range(10000):
        perm["disease"] = rng.permutation(primary_frame["disease"].values)
        coef = smf.ols(f"junction_ratio ~ disease + {base}", data=perm).fit().params["disease"]
        hits += abs(coef) >= abs(observed)
    rows.append({"analysis": "Sensitivity (b): label-permutation p, 10,000 perms, seed 0",
                 "n": len(primary_frame), "coef_log2": observed, "p": (hits + 1) / 10001})
    mw = mannwhitneyu(primary_frame.loc[primary_frame.disease == 1, "junction_ratio"],
                      primary_frame.loc[primary_frame.disease == 0, "junction_ratio"], alternative="two-sided")
    rows.append({"analysis": "Sensitivity (c): unadjusted Mann-Whitney, disease vs iPTB",
                 "n": len(primary_frame), "p": float(mw.pvalue)})
    pe_frame = data[data["group"].isin({"EOPE", "EOPE+FGR", "iPTB"})].copy()
    rows.append({"analysis": "Secondary 1: PE vs iPTB", **fit(pe_frame, "PE")})
    fgr_frame = data[data["group"].isin({"FGR", "iPTB"})].copy()
    fgr_frame["FGR"] = (fgr_frame["group"] == "FGR").astype(int)
    rows.append({"analysis": "Secondary 2: FGR-only vs iPTB", **fit(fgr_frame, "FGR")})
    pf_frame = data[data["group"].isin(DISEASE)].copy()
    rows.append({"analysis": "Secondary 3: PE vs FGR-only (J1 replication)", **fit(pf_frame, "PE")})
    rows.append({"analysis": "Secondary 4: disease vs iPTB, + log2 FLT1 CPM + hypoxia",
                 **fit(primary_frame, "disease", extra=" + log_flt1 + hypoxia")})
    for t, markers in EXPLORATORY.items():
        detected = float((data[markers] >= 5).all(axis=1).mean())
        if detected >= 0.70:
            rows.append({"analysis": f"Exploratory: {t}", "detected_fraction": detected,
                         **fit(primary_frame, "disease", outcome=f"expl_{t}")})
        else:
            rows.append({"analysis": f"Exploratory: {t} not detectable, not modelled", "detected_fraction": detected})
    result = pd.DataFrame(rows)
    result.to_csv(HERE / "gse203507-contrast-results.csv", index=False)
    summary = data.groupby("group")["junction_ratio"].agg(["count", "median", "mean", "std"])
    summary.to_csv(HERE / "gse203507-ratio-by-group-descriptive.csv")
    data[["run", "sample_id", "group", "ga_delivery_weeks", "fetal_sex", "junction_ratio", "body_ratio",
          "log_flt1", "hypoxia"]].to_csv(HERE / "gse203507-ratio-per-sample.csv", index=False)
    pd.set_option("display.width", 250)
    print(result.round(4).to_string(index=False))
    print(summary.round(3).to_string())
