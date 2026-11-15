"""POST HOC (after P was seen; can only limit): random-set null for P1 under the P1-cs adjustment (hyp_all + STB + EVT)
and within PE only. Same draws rule as p_confirm.py (decile-matched, seed 0). Usage: same args as p_confirm.py."""
from pathlib import Path
import json, sys
import numpy as np
HERE = Path(__file__).resolve().parent
src = (HERE / "p_confirm.py").read_text()
ns = {"__file__": str(HERE / "p_confirm.py"), "__name__": "p"}
exec(compile(src[: src.index("rx_cov = ")], "p_prefix", "exec"), ns)
z, dis, pe, by, need, resid, stats, res = ns["z"], ns["dis"], ns["pe"], ns["by"], ns["need"], ns["resid"], ns["stats"], ns["res"]
rng = np.random.default_rng(0)
out = {}
for name, frame, covs, obs in [("cs_hyp_STB_EVT", dis, ["hyp_all", "stb", "evt"], res["P1_cs_STB_EVT"]["rho"]),
                               ("PE_only_hyp", pe, ["hyp_all"], res["P1_PE_only"]["rho"])]:
    null = []
    for _ in range(1000):
        pick = np.concatenate([rng.choice(by[k], n, replace=False) for k, n in need.items()])
        tmp = frame.assign(rs=z.loc[frame.index, list(pick)].mean(axis=1))
        null.append(stats.spearmanr(resid(tmp, "rs", covs), resid(tmp, "FLT1", covs)).statistic)
    null = np.array(null)
    out[name] = {"observed": obs, "null_median": float(np.median(null)), "null_p95": float(np.percentile(null, 95)),
                 "percentile_of_observed": float((null < obs).mean() * 100)}
(HERE / "P_posthoc_null_cs.json").write_text(json.dumps(out, indent=2) + "\n")
print(json.dumps(out, indent=2))
