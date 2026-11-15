"""Group-blinded donor-coverage gate for the FLT1 isoform analysis (raw-reanalysis-plan.md + Addendum J1).

Reads per-run marker counts and evaluates gate criteria 1-5. Clinical labels are used ONLY to count
usable donors per group (criterion 5). No group-wise ratio, mean, test or plot is produced here.

Usage: python evaluate_gate.py --work WORK
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
SFLT = "ENST00000541932_unique_junction_1"
CAN1 = "ENST00000282397_unique_junction_1"
CAN2 = "ENST00000282397_unique_junction_2"
SFLT_BODY = "ENST00000541932_terminal_body_1"
CAN_BODY = "ENST00000282397_terminal_body_1"
SFLT_BODY_EFF = 567 - 49
CAN_BODY_EFF = 2782 - 49

parser = argparse.ArgumentParser()
parser.add_argument("--work", required=True, type=Path)
args = parser.parse_args()

manifest = pd.read_csv(HERE / "sra-sample-manifest.csv")
rows = []
for run in manifest["run"]:
    path = args.work / "counts" / f"{run}.json"
    record = json.loads(path.read_text())
    assert record["total_reads"] == int(manifest.loc[manifest["run"] == run, "spots"].iloc[0]), run
    row = {"run": run, "total_reads": record["total_reads"], "candidate_reads": record["candidate_reads"]}
    row.update({f"{m}": c for m, c in record["hisat2_unique_counts_flank10"].items()})
    row.update({f"{m}__flank5": c for m, c in record["hisat2_unique_junction_counts_flank5"].items()})
    row.update({f"{m}__exact50": c for m, c in record["exact_50mer_counts"].items()})
    row["unique_reads_in_flt1_span"] = record["alignment_stats"]["unique_reads_in_flt1_span"]
    row["flt1_span_forward_strand"] = record["alignment_stats"]["flt1_span_forward_strand"]
    row["flt1_span_reverse_strand"] = record["alignment_stats"]["flt1_span_reverse_strand"]
    rows.append(row)
counts = pd.DataFrame(rows)
assert len(counts) == 79
counts.to_csv(HERE / "cohort-marker-counts.csv", index=False)  # no group labels in this file

table = counts.merge(manifest[["run", "gate_population", "primary_contrast_group"]], on="run")
disease = table[table["gate_population"] != "control"].copy()
assert len(disease) == 58

each_junction = (disease[SFLT] >= 5) & (disease[CAN1] >= 5) & (disease[CAN2] >= 5)
depth = disease[SFLT] + disease[CAN1] + disease[CAN2]
enough_depth = depth >= 20
usable = each_junction & enough_depth

junction_ratio = np.log2((disease[SFLT] + 0.5) / ((disease[CAN1] + disease[CAN2]) / 2 + 0.5))
body_ratio = np.log2((disease[SFLT_BODY] / SFLT_BODY_EFF + 0.5 / SFLT_BODY_EFF)
                     / (disease[CAN_BODY] / CAN_BODY_EFF + 0.5 / CAN_BODY_EFF))
rho = float(spearmanr(junction_ratio, body_ratio).statistic)

mappability = json.loads((HERE / "hisat2-marker-mappability.json").read_text())["markers"]
primary_markers = [SFLT, CAN1, CAN2]
criterion1 = all(mappability[m]["fraction"] >= 0.90 for m in primary_markers)

group_usable = {
    group: {"donors": int(len(index)), "usable": int(usable.loc[index].sum()),
            "fraction": round(float(usable.loc[index].mean()), 4)}
    for group, index in disease.groupby("primary_contrast_group").groups.items()
}
criterion5 = all(g["fraction"] >= 0.60 for g in group_usable.values())

forward = int(table["flt1_span_forward_strand"].sum())
reverse = int(table["flt1_span_reverse_strand"].sum())


def quantiles(series):
    return {q: float(series.quantile(p)) for q, p in
            [("min", 0), ("q25", 0.25), ("median", 0.5), ("q75", 0.75), ("max", 1)]}


result = {
    "as_of": "2026-11-15",
    "population": "58 disease donors (PE/IUGR labels not used for criteria 1-4)",
    "criterion_1_primary_marker_unique_mappability": {
        "fractions": {m: mappability[m]["fraction"] for m in primary_markers}, "pass": bool(criterion1)},
    "criterion_2_each_primary_junction_ge5": {
        "donors_meeting": int(each_junction.sum()), "fraction": round(float(each_junction.mean()), 4),
        "threshold": 0.70, "pass": bool(each_junction.mean() >= 0.70),
        "per_junction_fraction_ge5": {m: round(float((disease[m] >= 5).mean()), 4) for m in primary_markers}},
    "criterion_3_total_primary_junction_depth_ge20": {
        "donors_meeting": int(enough_depth.sum()), "fraction": round(float(enough_depth.mean()), 4),
        "threshold": 0.70, "pass": bool(enough_depth.mean() >= 0.70)},
    "criterion_4_body_vs_junction_ratio_spearman": {
        "rho": round(rho, 4), "threshold": 0.70, "pass": bool(rho >= 0.70)},
    "criterion_5_usable_fraction_per_group": {"groups": group_usable, "threshold": 0.60, "pass": bool(criterion5)},
    "coverage_summary_disease_donors": {m: quantiles(disease[m]) for m in
                                        [SFLT, CAN1, CAN2, SFLT_BODY, CAN_BODY]},
    "strand_check_all_runs_unique_reads_in_flt1_span": {
        "forward": forward, "reverse": reverse,
        "reverse_fraction": round(reverse / (forward + reverse), 4)},
}
result["gate_pass"] = all(result[k]["pass"] for k in result if k.startswith("criterion_"))
(HERE / "cohort-coverage-gate.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
