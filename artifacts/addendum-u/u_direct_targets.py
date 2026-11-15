"""Addendum U: S2-type within-donor test in intact tissue (Addendum S code) for ARNT-bound vs unbound program genes.

Usage: python u_direct_targets.py SPATIAL_DIR   (the Addendum S data directory)
"""

import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "addendum-s"))
import s_spatial as S  # noqa: E402

d = Path(sys.argv[1])
tables = S.donor_tables(d)
spots = {x: int(t["m"].shape[1]) for x, t in tables.items()}
eligible = sorted(x for x, v in spots.items() if v >= S.MIN_SPOTS)
out = {"eligible": eligible}
for name in ("bound", "unbound"):
    genes = set((HERE / f"HIF2dep_ARNT_{name}.txt").read_text().split())
    S.HIF2 = genes  # analyse() reads the module-level HIF2 set
    r = S.analyse(tables, eligible)["S2_descriptive"]
    if name == "bound":
        r["supported"] = bool(r["sign_flip_p"] < 0.05 and r["mean_rho"] > r["random_null_p95"])
    out[f"U_{name}"] = r
b, u = out["U_bound"]["per_donor_rho"], out["U_unbound"]["per_donor_rho"]
out["bound_minus_unbound"] = {"per_donor": {x: b[x] - u[x] for x in eligible},
                              "mean": sum(b[x] - u[x] for x in eligible) / len(eligible),
                              "n_bound_higher": sum(b[x] > u[x] for x in eligible)}
(HERE / "U_results.json").write_text(json.dumps(out, indent=2) + "\n")
for k in ("U_bound", "U_unbound"):
    v = out[k]
    print(k, "n_genes_used", "mean %.3f pos %d/%d signflip %.4f p95 %.3f pct %.1f" % (
        v["mean_rho"], v["n_positive"], len(v["per_donor_rho"]), v["sign_flip_p"], v["random_null_p95"],
        v["random_percentile"]), "supported" if v.get("supported") else "")
print("bound-unbound mean %.3f, bound higher in %d/%d" % (out["bound_minus_unbound"]["mean"],
      out["bound_minus_unbound"]["n_bound_higher"], len(eligible)))
