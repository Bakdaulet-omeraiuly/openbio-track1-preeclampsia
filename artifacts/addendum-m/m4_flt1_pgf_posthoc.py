"""Addendum M, POST HOC (written after M1-M3 were seen): placental FLT1:PGF log-ratio, descriptive only.

Not confirmatory: the M expression data had already been opened when this was proposed.
Usage: python m4_flt1_pgf_posthoc.py DIR   (same inputs as m_confirm.py)
"""

import ast
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

HERE = Path(__file__).resolve().parent
DIR = Path(sys.argv[1])
# reuse the reader functions from m_confirm.py without running its analysis
src = (HERE / "m_confirm.py").read_text()
tree = ast.parse(src)
ns = {}
exec("import gzip\nimport numpy as np\nimport pandas as pd\n", ns)
for node in tree.body:
    if isinstance(node, ast.FunctionDef) and node.name in ("read_series", "read_gpl", "gene_matrix", "maybe_log"):
        exec(compile(ast.Module([node], []), "m_confirm", "exec"), ns)
read_series, read_gpl, gene_matrix, maybe_log = (ns[k] for k in ("read_series", "read_gpl", "gene_matrix", "maybe_log"))

spec = {
    "GSE44711": ("GPL10558.txt", "Symbol", True,
                 lambda m: m["ch1:condition"].map({"EOPET": "PE", "Control": "control"}), " + GA",
                 lambda m: m["ch1:gestational age (weeks)"].astype(float)),
    "GSE74341": ("GPL16699.txt", "GENE_SYMBOL", False,
                 lambda m: m["ch2:subject status"].map(lambda v: "PE" if v.startswith("patients with early-onset") else
                                                       ("control" if v.startswith("controls who delivered preterm") else "other")),
                 "", None),
    "GSE66273": ("GPL4133.txt", "GENE_SYMBOL", True,
                 lambda m: m["ch1:condition"].map(lambda v: "control" if v.startswith("Preterm control") else "PE"),
                 " + GA", lambda m: m["ch1:gestational weeks"].astype(float)),
}
out = {}
for c, (gpl, col, uselog, grpf, cov, gaf) in spec.items():
    m, e = read_series(DIR / f"{c}_series_matrix.txt.gz")
    if uselog:
        e, _ = maybe_log(e)
    g = gene_matrix(e, read_gpl(DIR / gpl, col))
    d = pd.DataFrame({"ratio": g.loc["FLT1"] - g.loc["PGF"], "group": grpf(m)})
    if gaf is not None:
        d["GA"] = gaf(m)
    f = d[d.group.isin(["PE", "control"])].assign(PE=lambda x: (x.group == "PE").astype(int))
    fit = smf.ols(f"ratio ~ PE{cov}", data=f).fit()
    out[c] = {"diff_log2": float(fit.params["PE"]), "se": float(fit.bse["PE"]),
              "p_one_sided_normal": float(1 - stats.norm.cdf(fit.params["PE"] / fit.bse["PE"])),
              "median_PE": float(f.ratio[f.PE == 1].median()), "median_control": float(f.ratio[f.PE == 0].median()),
              "n": [int(f.PE.sum()), int((f.PE == 0).sum())]}
villous = [out[c] for c in ("GSE44711", "GSE74341")]
w = np.array([1 / v["se"] ** 2 for v in villous]); est = np.array([v["diff_log2"] for v in villous])
fe = float((w * est).sum() / w.sum()); se = float(1 / np.sqrt(w.sum()))
out["pooled_villous_descriptive"] = {"diff_log2": fe, "CI_low": fe - 1.96 * se, "CI_high": fe + 1.96 * se}
out["note"] = ("Post hoc, descriptive. GSE74341 values are log2 ratios to a common pooled reference, so FLT1-PGF is "
               "a relative ratio; the PE-vs-control difference is still interpretable on the log2 scale.")
(HERE / "M4_flt1_pgf_posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
