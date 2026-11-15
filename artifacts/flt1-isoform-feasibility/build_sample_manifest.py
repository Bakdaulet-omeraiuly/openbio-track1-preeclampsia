"""Join official GEO series metadata to the saved SRA run table."""

import argparse
import csv
import gzip
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent


def read_series_rows(path):
    wanted = {}
    disease = None
    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if not row:
                continue
            if row[0] in {"!Sample_title", "!Sample_geo_accession"}:
                wanted[row[0]] = row[1:]
            elif row[0] == "!Sample_characteristics_ch1" and row[1].startswith("disease group:"):
                disease = [value.removeprefix("disease group: ") for value in row[1:]]
    assert set(wanted) == {"!Sample_title", "!Sample_geo_accession"}
    assert disease is not None
    assert len(wanted["!Sample_title"]) == len(wanted["!Sample_geo_accession"]) == len(disease) == 79
    return pd.DataFrame({
        "geo_title": wanted["!Sample_title"],
        "gsm": wanted["!Sample_geo_accession"],
        "disease_group": disease,
    })


parser = argparse.ArgumentParser()
parser.add_argument("--series-matrix-gz", required=True, type=Path)
args = parser.parse_args()

geo = read_series_rows(args.series_matrix_gz)
runinfo = pd.read_csv(HERE / "sra-runinfo.csv", usecols=["Run", "SampleName", "spots", "bases", "avgLength"])
manifest = runinfo.merge(geo, left_on="SampleName", right_on="gsm", validate="one_to_one")
assert len(manifest) == 79
manifest["gate_population"] = manifest.disease_group.map(
    lambda value: "control" if value == "Control" else "disease"
)
manifest["primary_contrast_group"] = manifest.disease_group.map({
    "Control": "descriptive_control",
    "Preeclampsia only": "PE",
    "Preeclampsia and Intrauterine Growth Restriction": "PE",
    "Intrauterine Growth Restriction only": "IUGR_only",
})
assert manifest.gate_population.value_counts().to_dict() == {"disease": 58, "control": 21}
assert manifest.primary_contrast_group.value_counts().to_dict() == {
    "PE": 40,
    "descriptive_control": 21,
    "IUGR_only": 18,
}
manifest.rename(columns={"Run": "run", "SampleName": "gsm_from_sra"}).to_csv(
    HERE / "sra-sample-manifest.csv", index=False
)
print(manifest.groupby(["gate_population", "primary_contrast_group"]).size().to_string())
