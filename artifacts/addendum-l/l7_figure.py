"""Figure for Addendum L (descriptive; reads committed JSON outputs and recomputes the L4 null draws).

Left: C4 (HIF2dep vs FLT1 | hypoxia) pooled and per cohort, without and with composition adjustment.
Right: the L4 random-gene-set null (same seed and draws as l4_calibration_posthoc.py) vs observed.
Usage: python l7_figure.py DIR   (same DIR as l3_confirm.py)
"""

import json
from pathlib import Path
import runpy
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
L3 = json.loads((HERE / "L3_confirmation.json").read_text())["main"]
L5 = json.loads((HERE / "L5_stb_adjustment_posthoc.json").read_text())
L6 = json.loads((HERE / "L6_deconvolution_posthoc.json").read_text())

# Re-run L4 in memory to get the null draws (it rewrites L4_calibration_posthoc.json identically).
sys.argv = [str(HERE / "l4_calibration_posthoc.py"), sys.argv[1]]
l4 = runpy.run_path(str(HERE / "l4_calibration_posthoc.py"), run_name="__main__")
null, observed = l4["null"], l4["observed"]

COH = ["GSE35574", "GSE24129", "GSE147776"]
COL = {"GSE35574": "#2a78d6", "GSE24129": "#eb6834", "GSE147776": "#1baf7a"}
MRK = {"GSE35574": "o", "GSE24129": "s", "GSE147776": "^"}
rows = [
    ("Hypoxia only\n(pre-specified C4)", {"pooled_rho": L3["C4"]["pooled_rho"], "CI_low": L3["C4"]["CI_low"],
                                         "CI_high": L3["C4"]["CI_high"]},
     {c: L3["per_cohort"][c]["C4_partial_rho_n"][0] for c in COH}),
    ("+ 4 STB marker genes\n(post hoc)", L5["primary_hyp_plus_STB"],
     {c: L5["primary_hyp_plus_STB"]["per_cohort"][c]["rho"] for c in COH}),
    ("+ deconvolved STB, CTB\n(post hoc)", L6["primary_hyp_STB_CTB"],
     {c: L6["primary_hyp_STB_CTB"]["per_cohort"][c]["rho"] for c in COH}),
    ("+ deconvolved STB, CTB, EVT\n(post hoc)", L6["secondary_hyp_STB_CTB_EVT"],
     {c: L6["secondary_hyp_STB_CTB_EVT"]["per_cohort"][c]["rho"] for c in COH}),
]

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e"})
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.2), gridspec_kw={"width_ratios": [1.25, 1]})
fig.patch.set_facecolor("#fcfcfb")
for ax in (a, b):
    ax.set_facecolor("#fcfcfb")

y = np.arange(len(rows))[::-1]
a.axvline(0, color="#52514e", lw=0.8)
for yi, (lab, pooled, per) in zip(y, rows):
    a.plot([pooled["CI_low"], pooled["CI_high"]], [yi, yi], color="#0b0b0b", lw=2, solid_capstyle="round")
    a.plot(pooled["pooled_rho"], yi, "D", color="#0b0b0b", ms=8)
    a.text(pooled["CI_high"] + 0.03, yi, f"{pooled['pooled_rho']:.2f} [{pooled['CI_low']:.2f}, {pooled['CI_high']:.2f}]",
           va="center", fontsize=8, color="#0b0b0b")
    for j, c in enumerate(COH):
        a.plot(per[c], yi - 0.22 + 0.11 * j, MRK[c], color=COL[c], ms=8, mec="#fcfcfb", mew=1.5,
               label=c if yi == y[0] else None)
a.set_yticks(y); a.set_yticklabels([r[0] for r in rows], fontsize=8)
a.set_xlim(-0.3, 1.05); a.set_xlabel("Partial Spearman rho, HIF-2 program vs FLT1 (PE + FGR-only)")
a.set_title("Link holds when pooled; one small cohort turns negative\nwith full composition adjustment", fontsize=9.5, loc="left")
a.legend(title="Cohort (black diamond = pooled, line = 95% CI)", fontsize=7.5, title_fontsize=7.5, frameon=False,
         loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3)

b.hist(null, bins=35, color="#86b6ef", edgecolor="#fcfcfb", linewidth=1)
b.axvline(observed, color="#0b0b0b", lw=2)
b.text(observed - 0.02, b.get_ylim()[1] * 0.92, f"HIF-2 program\n{observed:.2f}\n(79 genes on all\nthree platforms)", ha="right", va="top", fontsize=8)
b.set_xlabel("Pooled partial rho vs FLT1 (| hypoxia)"); b.set_ylabel("Random 100-gene sets")
b.set_title(f"Above all 1,000 expression-matched random sets\n(post hoc; random median {np.median(null):.2f})",
            fontsize=9.5, loc="left")
fig.suptitle("Addendum L: genes lowered by HIF-2 blockade track placental FLT1 in three untouched cohorts (138 placentas)",
             x=0.01, ha="left", fontsize=10.5)
fig.tight_layout()
out = ROOT / "figures/fig16_hif2_program_flt1.png"
fig.savefig(out, dpi=300, bbox_inches="tight", facecolor=fig.get_facecolor())
print(out)
