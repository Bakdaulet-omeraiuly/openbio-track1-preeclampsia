"""Addendum W, POST HOC (written after W2 was seen; descriptive only): random-set calibration of the W2 MVM effects.

  * HIF2dep: MVM coefficient vs 1,000 random gene sets of the same size, decile-matched on mean expression,
    excluding hypoxia, HIF2dep and FLT1 genes (seed 0). Same scoring and model as W1.
  * FLT1: MVM coefficient vs 1,000 random single genes from FLT1's mean-expression decile (seed 0), plus FLT1's
    rank among all expressed genes.
Pre-stated reading: "calibrated" only if the observed value is above the random 95th percentile. Because this was
chosen after seeing W2, any result is a lead, not a confirmation.
Usage: python w2_calibration_posthoc.py DIR
"""

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import w_mvm as W  # noqa: E402

d = Path(sys.argv[1])
counts, m = W.load(d)
cpm = counts / counts.sum() * 1e6
lg = np.log2(cpm[(cpm >= 1).sum(axis=1) >= 0.1 * cpm.shape[1]] + 1)
z = lg.sub(lg.mean(axis=1), axis=0).div(lg.std(axis=1), axis=0)
dec = pd.qcut(lg.mean(axis=1).rank(method="first"), 10, labels=False)
pool = dec[~dec.index.isin(W.HYP | W.HIF2 | {"FLT1"})]
out = {}

obs = W.fit(W.score(z, W.HIF2), m)[0]
need = dec[z.index.intersection(sorted(W.HIF2))].value_counts()
rng = np.random.default_rng(0); null = []
for _ in range(1000):
    g = np.concatenate([rng.choice(pool.index[pool == k], size=n, replace=False) for k, n in need.items()])
    null.append(W.fit(W.score(z, set(g)), m)[0])
null = np.array(null)
out["HIF2dep"] = {"observed": obs, "random_median": float(np.median(null)), "random_p95": float(np.percentile(null, 95)),
                  "percentile": float((null < obs).mean() * 100), "calibrated": bool(obs > np.percentile(null, 95))}

zf = lambda g: (z.loc[g] - z.loc[g].mean()) / z.loc[g].std()
obs_f = W.fit(zf("FLT1"), m)[0]
same = pool.index[pool == dec["FLT1"]]
rng = np.random.default_rng(0)
null_f = np.array([W.fit(zf(g), m)[0] for g in rng.choice(same, size=min(1000, len(same)), replace=False)])
out["FLT1"] = {"observed": obs_f, "random_median": float(np.median(null_f)), "random_p95": float(np.percentile(null_f, 95)),
               "percentile_in_decile": float((null_f < obs_f).mean() * 100), "calibrated": bool(obs_f > np.percentile(null_f, 95))}
allc = pd.Series({g: W.fit(zf(g), m)[0] for g in z.index})
out["FLT1"]["rank_among_all_genes"] = int((allc > obs_f).sum() + 1)
out["FLT1"]["n_genes"] = int(len(allc))
out["FLT1"]["genome_percentile"] = float((allc < obs_f).mean() * 100)
out["all_genes_median_MVM_coef"] = float(allc.median())
(Path(__file__).resolve().parent / "W2_calibration_posthoc.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=1))
