"""Link SRA runs to processed GSE114691 columns and build the pre-specified covariates (Addendum J1, item 9).

GEO `Sample_description` holds the donor code used as the processed-matrix column name
(e.g. `267-05` -> `X267.05`); the second description row names the matrix file, which must agree. Outputs one row per run with gene-level FLT1 CPM and the existing
hypoxia score (`hyp_all` from notebooks/results/C1_scores_GSE114691.csv). No isoform data is read here.

Usage: python build_covariates.py
"""

import gzip
import io
from pathlib import Path
import re
import urllib.request

import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
EXT = ROOT / "data_external"
URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114691/matrix/GSE114691_series_matrix.txt.gz"

text = gzip.decompress(urllib.request.urlopen(URL).read()).decode()
fields = {}
for line in text.splitlines():
    if line.startswith("!Sample_geo_accession"):
        fields["gsm"] = [v.strip('"') for v in line.split("\t")[1:]]
    elif line.startswith("!Sample_description"):
        values = [v.strip('"') for v in line.split("\t")[1:]]
        fields["matrix_file" if values[0].startswith("MasterCount_") else "description"] = values
geo = pd.DataFrame(fields)


def column_name(description):
    # GEO writes e.g. `267-05`, `PE 12 -03`, `396-10_a`; the matrices use `X267.05`, `PE12.03`, `X396.10`.
    name = re.sub(r"_a$", "", re.sub(r"\s+", "", description)).replace("-", ".")
    return "X" + name if name[0].isdigit() else name


FILE_GROUP = {"MasterCount_ControlONLY.txt": "control", "MasterCount_IUGRONLY.txt": "IUGR",
              "MasterCount_PEONLY.txt": "PE", "MasterCount_PEIUGR.txt": "PE+IUGR"}
geo["sample"] = geo["matrix_file"].map(FILE_GROUP) + "|" + geo["description"].map(column_name)

tx2gene = pd.read_csv(EXT / "grch37_tx2gene.tsv", sep="\t").set_index("Transcript stable ID")["Gene name"]
parts = []
for group, name in {"control": "ControlONLY", "IUGR": "IUGRONLY", "PE": "PEONLY", "PE+IUGR": "PEIUGR"}.items():
    part = pd.read_csv(EXT / f"GSE114691_MasterCount_{name}.txt.gz", sep="\t", index_col=0)
    part.columns = [f"{group}|{c}" for c in part.columns]
    parts.append(part)
tx = pd.concat(parts, axis=1)
flt1_rows = tx.index[tx.index.map(lambda t: tx2gene.get(t) == "FLT1")]
cpm = tx.loc[flt1_rows].sum() / tx.sum() * 1e6
scores = pd.read_csv(ROOT / "notebooks/results/C1_scores_GSE114691.csv", index_col=0)

processed = pd.DataFrame({"sample": tx.columns})
processed["total_FLT1_CPM"] = cpm.reindex(processed["sample"]).values
processed["hypoxia_score"] = scores["hyp_all"].reindex(processed["sample"]).values
assert processed["sample"].is_unique

manifest = pd.read_csv(HERE / "sra-sample-manifest.csv")[["run", "gsm"]]
out = manifest.merge(geo, on="gsm", how="left", validate="one_to_one").merge(
    processed, on="sample", how="left", validate="one_to_one", indicator=True)
print("FLT1 transcript rows:", len(flt1_rows), "| unmatched runs:", int((out["_merge"] != "both").sum()))
assert len(out) == 79 and (out["_merge"] == "both").all() and out[["total_FLT1_CPM", "hypoxia_score"]].notna().all().all()
out[["run", "gsm", "description", "sample", "total_FLT1_CPM", "hypoxia_score"]].to_csv(
    HERE / "run-covariates.csv", index=False)
print(out.head(3).to_string())
