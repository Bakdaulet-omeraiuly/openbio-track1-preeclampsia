"""Pre-specified PE-versus-FGR isoform-composition contrast (raw-reanalysis-plan.md + Addendum J1).

This script was written before any cohort counts existed. It refuses to run unless the group-blinded
coverage gate (`cohort-coverage-gate.json`) passed.

Usage: python primary_contrast.py
"""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

HERE = Path(__file__).resolve().parent
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

gate = json.loads((HERE / "cohort-coverage-gate.json").read_text())
if not gate["gate_pass"]:
    sys.exit("Coverage gate did not pass: the clinical contrast is not run (pre-specified).")

counts = pd.read_csv(HERE / "cohort-marker-counts.csv")
manifest = pd.read_csv(HERE / "sra-sample-manifest.csv")[["run", "disease_group", "gate_population",
                                                          "primary_contrast_group"]]
covariates = pd.read_csv(HERE / "run-covariates.csv")[["run", "total_FLT1_CPM", "hypoxia_score"]]
data = counts.merge(manifest, on="run", validate="one_to_one").merge(covariates, on="run", validate="one_to_one")
canonical = (data[CAN1] + data[CAN2]) / 2
data["junction_ratio"] = np.log2((data[SFLT] + 0.5) / (canonical + 0.5))
data["body_ratio"] = np.log2((data[SFLT_BODY] / SFLT_BODY_EFF + 0.5 / SFLT_BODY_EFF)
                             / (data[CAN_BODY] / CAN_BODY_EFF + 0.5 / CAN_BODY_EFF))
data["log_flt1"] = np.log2(data["total_FLT1_CPM"] + 1)
for transcript, markers in EXPLORATORY.items():
    data[f"expl_{transcript}"] = np.log2((data[markers].mean(axis=1) + 0.5) / (canonical + 0.5))

disease = data[data["gate_population"] != "control"].copy()
disease["PE"] = (disease["primary_contrast_group"] == "PE").astype(int)
assert len(disease) == 58 and disease["PE"].sum() == 40


def fit(frame, outcome, formula_rhs="PE + log_flt1 + hypoxia_score"):
    model = smf.ols(f"{outcome} ~ {formula_rhs}", data=frame).fit(cov_type="HC3")
    low, high = model.conf_int().loc["PE"]
    return {"outcome": outcome, "model": formula_rhs, "n": int(model.nobs), "n_PE": int(frame["PE"].sum()),
            "PE_coefficient_log2": float(model.params["PE"]), "CI_low": float(low), "CI_high": float(high),
            "p": float(model.pvalues["PE"])}


rows = [{"analysis": "PRIMARY: junction ratio, PE (PE only + PE with IUGR) vs IUGR only", **fit(disease, "junction_ratio")}]
rows.append({"analysis": "Secondary: body-read ratio", **fit(disease, "body_ratio")})
pe_only = disease[disease["disease_group"].isin(["Preeclampsia only"]) | (disease["PE"] == 0)]
rows.append({"analysis": "Secondary: PE only vs IUGR only, junction ratio", **fit(pe_only, "junction_ratio")})
rows.append({"analysis": "Descriptive: unadjusted junction ratio", **fit(disease, "junction_ratio", "PE")})
exploratory = [{"analysis": f"Exploratory: {t} vs canonical", **fit(disease, f"expl_{t}")} for t in EXPLORATORY]
adjusted = multipletests([row["p"] for row in exploratory], method="fdr_bh")[1]
for row, q in zip(exploratory, adjusted):
    row["BH_q"] = float(q)
result = pd.DataFrame(rows + exploratory)
result.to_csv(HERE / "primary-contrast-results.csv", index=False)

summary = data.groupby("disease_group")[["junction_ratio", "body_ratio"]].agg(["count", "median", "mean", "std"])
summary.to_csv(HERE / "isoform-ratio-by-group-descriptive.csv")
data[["run", "disease_group", "junction_ratio", "body_ratio", "log_flt1", "hypoxia_score"]].to_csv(
    HERE / "isoform-ratio-per-donor.csv", index=False)
pd.set_option("display.width", 250)
print(result.round(4).to_string(index=False))
print(summary.round(3).to_string())
