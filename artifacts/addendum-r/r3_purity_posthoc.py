"""Addendum R, POST HOC purity check (written after R1 was seen, and after a marker-only diagnostic showed that
about 33% of the selected barcodes also carry a non-trophoblast lineage marker and only about 65% are KRT7+).

Rule (fixed before running):
  * Barcodes: the R1 STB selection AND KRT7 detected AND none of PTPRC, CD14, PECAM1, COL1A1, HBB detected.
  * Covariates: as in R1 (log10 UMI, log10 genes, mito fraction), plus a library indicator for donors with
    two libraries.
  * Same gate (>= 30 barcodes in >= 8 donors), sign-flip test and 1,000 random-set null as R1.
Reading (pre-stated):
  * Mean rho > 0 with sign-flip p < 0.05 AND above the random-set 95th percentile: the within-STB link stands
    for marker-pure trophoblast barcodes.
  * Otherwise the R1 result may reflect contamination or cell mix, and we say so.
Usage: python r3_purity_posthoc.py DIR
"""

import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r_within_donor as R  # noqa: E402

LINEAGE = ["PTPRC", "CD14", "PECAM1", "COL1A1", "HBB"]
orig_stb = R.stb_cells
library_of = {}  # donor -> list of per-cell library index, filled while loading


def pure_stb(m, genes, min_markers=3):
    sel, total, ngenes, mito = orig_stb(m, genes, min_markers)
    det = lambda g: np.asarray(m[np.where(genes == g)[0][0]].todense()).ravel() > 0
    lineage = np.any([det(g) for g in LINEAGE], axis=0)
    return sel & det("KRT7") & ~lineage, total, ngenes, mito


R.stb_cells = pure_stb
d = Path(sys.argv[1])
tables = R.donor_tables(d)
# library indicator: rebuild per donor by re-selecting per library (same rule, same order as donor_tables)
for donor, rows in R.LIBS.groupby("donor"):
    lab = []
    for i, (_, r) in enumerate(rows.iterrows()):
        m, genes = R.load_library(d, r)
        lab.append(np.full(int(pure_stb(m, genes)[0].sum()), float(i)))
    library_of[donor] = np.concatenate(lab)
    assert len(library_of[donor]) == tables[donor]["m"].shape[1]

counts = {x: int(t["m"].shape[1]) for x, t in tables.items()}
eligible = sorted(x for x, n in counts.items() if n >= R.MIN_STB_CELLS)
out = {"barcodes_per_donor": counts, "eligible_donors": eligible, "gate_passed": len(eligible) >= R.MIN_DONORS}
if out["gate_passed"]:
    orig_cs = np.column_stack

    def run(donors):
        order = iter(donors)

        def patched(arrs):
            res = orig_cs(arrs)
            if res.ndim == 2 and res.shape[1] == 4 and np.allclose(res[:, 0], 1):
                x = next(order)
                lib = library_of[x]
                return orig_cs([res, lib]) if np.ptp(lib) > 0 else res
            return res

        R.np.column_stack = patched
        try:
            return R.analyse(tables, donors, "pure trophoblast barcodes + library")
        finally:
            R.np.column_stack = orig_cs

    res = run(eligible)
    res["R1"]["supported"] = bool(res["R1"]["sign_flip_p"] < 0.05 and res["R1"]["mean_rho"] > res["R1"]["random_null_p95"])
    out["R1_pure"] = res["R1"]
    out["R2_pure_descriptive"] = res["R2_descriptive"]
    pe = [x for x in eligible if x in R.PE]; ct = [x for x in eligible if x in R.CONTROL]
    out["FLT1_detection_rate_in_pure_barcodes"] = {
        x: float((np.asarray(tables[x]["m"][int(np.where(tables[x]["genes"] == "FLT1")[0][0])].todense()).ravel() > 0).mean())
        for x in eligible}
(Path(__file__).resolve().parent / "R_purity_posthoc.json").write_text(json.dumps(out, indent=2, default=float) + "\n")
print("gate", out["gate_passed"], "eligible", len(eligible), "barcodes", counts)
for k in ("R1_pure", "R2_pure_descriptive"):
    if k in out:
        v = out[k]
        print(k, "mean %.3f positive %d/%d sign-flip p %.4f null p95 %.3f percentile %.1f supported=%s" % (
            v["mean_rho"], v["n_positive"], len(v["per_donor_rho"]), v["sign_flip_p"], v["random_null_p95"],
            v["random_percentile"], v.get("supported")))
