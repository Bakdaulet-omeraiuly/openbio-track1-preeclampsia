"""Addendum R, POST HOC (written after R1/R2 were seen; it can only limit the claim).

Question: is the within-STB hypoxia–FLT1 link explained by STB maturity/identity, i.e. cells that are "more STB"
having more of both? We add the STB marker score (mean z of the 6 selection markers: CGA, CYP19A1, PSG1, PSG4,
KISS1, LGALS13) as a covariate, and also the number of markers detected (3–6). The same donors, gene universe,
sign-flip test and random-set null are used as in r_within_donor.py.

Pre-stated reading (written before running):
  - R1 still above the random-set 95th percentile and sign-flip p < 0.05 means the link is not explained by
    this marker score.
  - If it drops to the null, the R1 link may reflect the STB state, and we say so.
Usage: python r2_stb_state_posthoc.py DIR
"""

import json
from pathlib import Path
import sys

import numpy as np
from scipy import sparse

sys.path.insert(0, str(Path(__file__).resolve().parent))
import r_within_donor as R  # noqa: E402

d = Path(sys.argv[1])
tables = R.donor_tables(d)
eligible = sorted(x for x, t in tables.items() if t["m"].shape[1] >= R.MIN_STB_CELLS)
genes = tables[eligible[0]]["genes"]
marker_idx = [int(np.where(genes == g)[0][0]) for g in R.STB_MARKERS]
# analyse() builds each donor's covariates from total UMI, genes and mito; the wrapper below appends the
# two STB-state columns to that matrix without changing the committed primary script.
state = {}
for x in eligible:
    t = tables[x]
    lg = (sparse.diags(1e4 / np.maximum(t["total"], 1)) @ t["m"].T).tocsc()[:, marker_idx].log1p().toarray()
    z = (lg - lg.mean(0)) / np.where(lg.std(0) > 0, lg.std(0), np.nan)
    state[x] = np.column_stack([np.nanmean(z, axis=1), (lg > 0).sum(axis=1).astype(float)])

orig_column_stack = np.column_stack


def analyse_with_state(tables, donors, label):
    # Patch np.column_stack inside R only while building each donor's covariate matrix (4 columns: 1, UMI,
    # genes, mito). That call happens once per donor, in donor order.
    order = iter(donors)

    def patched(arrs):
        out = orig_column_stack(arrs)
        if out.ndim == 2 and out.shape[1] == 4 and np.allclose(out[:, 0], 1):
            x = next(order)
            assert len(out) == len(state[x])
            return orig_column_stack([out, state[x]])
        return out

    R.np.column_stack = patched
    try:
        return R.analyse(tables, donors, label)
    finally:
        R.np.column_stack = orig_column_stack


res = analyse_with_state(tables, eligible, "R1/R2 + STB marker score + n markers")
out = {"R1_adjusted_for_STB_state": res["R1"], "R2_adjusted_for_STB_state": res["R2_descriptive"]}
(Path(__file__).resolve().parent / "R_stb_state_posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
for k, v in out.items():
    print(k, "mean %.3f positive %d/%d sign-flip p %.4f null p95 %.3f percentile %.1f" % (
        v["mean_rho"], v["n_positive"], len(v["per_donor_rho"]), v["sign_flip_p"], v["random_null_p95"],
        v["random_percentile"]))
