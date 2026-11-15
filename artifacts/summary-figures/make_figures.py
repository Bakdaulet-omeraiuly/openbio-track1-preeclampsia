"""Summary figures for the article (fig17-fig20). Reads only committed result files; no raw data needed.

Usage: python artifacts/summary-figures/make_figures.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures"
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e6e5e1"
plt.rcParams.update({"font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})


def fisher_ci(r, n):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    return np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


# ---- fig17: hypoxia-FLT1 coupling within PE across 10 bulk cohorts
cohorts = [  # (label, rho, n, pre-registered?)  values from the committed results (see hypothesis card)
    ("GSE203507 (discovery)", 0.78, 15, False), ("GSE114691", 0.59, 40, True), ("GSE75010", 0.58, 80, True),
    ("GSE35574", 0.65, 19, True), ("GSE24129", 0.62, 8, True), ("GSE147776", 0.56, 13, True),
    ("GSE44711", 0.61, 8, True), ("GSE74341", 0.68, 7, True), ("GSE66273 (chorionic plate)", 0.26, 12, True),
    ("GSE186257", 0.69, 26, True)]
L = json.load(open(ROOT / "artifacts/addendum-l/L3_confirmation.json"))["main"]["per_cohort"]
assert abs(L["GSE35574"]["C3_rho_n"][0] - 0.65) < 0.01
fig, ax = plt.subplots(figsize=(8.6, 5.4))
for i, (lab, r, n, pre) in enumerate(cohorts[::-1]):
    lo, hi = fisher_ci(r, n)
    ax.plot([lo, hi], [i, i], color=BLUE if pre else MUTED, lw=2, solid_capstyle="round")
    ax.plot(r, i, "o", ms=8, color=BLUE if pre else "white", mec=BLUE if pre else MUTED, mew=2)
    ax.text(1.04, i, f"n = {n}", va="center", fontsize=9, color=MUTED)
ax.axvline(0, color=MUTED, lw=1)
ax.set_yticks(range(len(cohorts))); ax.set_yticklabels([c[0] for c in cohorts[::-1]])
ax.set_xlim(-0.6, 1.0)
ax.set_xlabel("How closely the low-oxygen alarm and FLT1 rise together in preeclampsia placentas\n(0 = not at all, 1 = perfectly; bars show the range of plausible values)")
ax.set_ylabel("Patient group (public dataset ID)")
fig.suptitle("The low-oxygen alarm rises with FLT1 in all 10 patient groups", x=0.02, ha="left", fontsize=12.5, color=INK)
fig.text(0.02, 0.01, "Filled dots: prediction written down before the data were opened. Open: where we first saw it. Small groups have wide bars.", fontsize=8.5, color=MUTED)
ax.grid(axis="x", color=GRID, lw=0.8); ax.set_axisbelow(True)
fig.tight_layout(rect=(0, 0.04, 1, 0.95)); fig.savefig(FIG / "fig17_coupling_across_cohorts.png", dpi=200)
fig.savefig(FIG / "headline.png", dpi=200); plt.close(fig)  # headline = the single most important finding

# ---- fig18: within-placenta links (R pure single cells; S intact tissue)
Rp = json.load(open(ROOT / "artifacts/addendum-r/R_purity_posthoc.json"))["R1_pure"]
S = json.load(open(ROOT / "artifacts/addendum-s/S_results.json"))["primary"]["S1"]
fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2), sharey=True)
for ax, res, title, color in ((axes[0], Rp, "Single trophoblast cells\n(19 placentas, GSE319236)", BLUE),
                              (axes[1], S, "Intact tissue slices\n(19 placentas, same study)", ORANGE)):
    vals = np.sort(np.array(list(res["per_donor_rho"].values())))
    x = np.random.default_rng(0).uniform(-0.12, 0.12, len(vals))
    ax.axhspan(-1, res["random_null_p95"], color=GRID, zorder=0)
    ax.axhline(0, color=MUTED, lw=1)
    ax.scatter(x, vals, s=60, color=color, edgecolor="#fcfcfb", linewidth=1.5, zorder=3)
    ax.plot([-0.25, 0.25], [res["mean_rho"]] * 2, color=INK, lw=2)
    ax.text(0.28, res["mean_rho"], f"mean {res['mean_rho']:.2f}\n{res['n_positive']}/{len(vals)} above 0", va="center", fontsize=9, color=INK)
    ax.set_xlim(-0.5, 0.75); ax.set_xticks([]); ax.set_title(title, fontsize=10.5, color=INK)
axes[0].set_ylim(-0.2, 0.45)
axes[0].set_ylabel("Within-placenta correlation,\nlow-oxygen score vs FLT1")
axes[1].text(-0.48, -0.17, "Grey band: up to the 95th percentile of 1,000 random gene lists", fontsize=8, color=MUTED)
fig.suptitle("Inside one placenta, cells and tissue spots with more alarm carry more FLT1", x=0.02, ha="left", fontsize=12.5, color=INK)
fig.tight_layout(); fig.savefig(FIG / "fig18_within_placenta.png", dpi=200); plt.close(fig)

# ---- fig19: HIF-2 vs HIF-1 knockouts at 3% oxygen (V)
V = json.load(open(ROOT / "artifacts/addendum-v/V_results.json"))
rows = []
for study, gene in (("GSE288288", "HIF-2α (EPAS1) removed"), ("GSE302997", "HIF-1α (HIF1A) removed")):
    for s in V[study]["per_sample"]:
        if s["condition"] == "hypoxia":
            rows.append({"panel": gene, "group": "removed" if s["null"] else "normal", "FLT1": s["FLT1"],
                         "program": s["program"], "hyp_all": s["hyp_all"]})
df = pd.DataFrame(rows)
outs = [("FLT1", "FLT1 (log2 CPM)"), ("program", "HIF-2 program score"), ("hyp_all", "Classic low-oxygen genes")]
fig, axes = plt.subplots(2, 3, figsize=(9, 5.6), sharey="col")
for i, panel in enumerate(["HIF-2α (EPAS1) removed", "HIF-1α (HIF1A) removed"]):
    for j, (col, lab) in enumerate(outs):
        ax = axes[i, j]; d = df[df.panel == panel]
        for k, (g, color) in enumerate((("normal", MUTED), ("removed", BLUE if i == 0 else ORANGE))):
            v = d[d.group == g][col].values
            x = k + np.random.default_rng(k).uniform(-0.1, 0.1, len(v))
            ax.scatter(x, v, s=45, color=color, edgecolor="#fcfcfb", linewidth=1.2, zorder=3)
            ax.plot([k - 0.22, k + 0.22], [v.mean()] * 2, color=INK, lw=2)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["normal", "removed"]); ax.set_xlim(-0.6, 1.6)
        if i == 0: ax.set_title(lab, fontsize=10.5, color=INK)
        if j == 0: ax.set_ylabel(panel, fontsize=10, color=INK)
        ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
fig.suptitle("At low oxygen, removing HIF-2 erases FLT1; removing HIF-1 does not", x=0.02, ha="left", fontsize=12.5, color=INK)
fig.text(0.02, 0.005, "Stem-cell-derived placental cells (one cell line), 3% oxygen; each dot is one culture. Data: GSE288288, GSE302997.\nScores are standardised within each study, so compare within a row only.", fontsize=8, color=MUTED)
fig.tight_layout(rect=(0, 0.05, 1, 0.95)); fig.savefig(FIG / "fig19_hif2_vs_hif1.png", dpi=200); plt.close(fig)

# ---- fig20: the chain, link by link
links = [("1", "Alarm higher with pathology signs of poor blood supply?", "No association found (2 pathology cohorts; Doppler)", "no"),
         ("2", "Alarm higher in early preeclampsia vs same-stage births?", "Yes, but fragile (5 of 6 matched cohorts)", "part"),
         ("3", "Does HIF drive the program genes in placental cells?", "Yes: drug, knockdown, ARNT binding, knockout", "yes"),
         ("4", "Is it HIF-2 (not HIF-1) that FLT1 needs?", "Yes, in a stem-cell model (one line); confirms Sasagawa 2018", "yes"),
         ("5", "Do the alarm and FLT1 rise together in placentas?", "Yes: 10/10 cohorts; within cells and intact tissue", "yes"),
         ("6", "Does placental FLT1 RNA track the mother's illness?", "No clear association in our data", "no")]
fig, ax = plt.subplots(figsize=(9.5, 6.6)); ax.axis("off")
mark = {"yes": ("✔", BLUE), "part": ("~", "#eda100"), "no": ("✘", ORANGE)}
for i, (k, link, ev, st) in enumerate(links):
    y = 0.95 - i * 0.15
    sym, col = mark[st]
    ax.text(0.0, y, k, fontsize=13, color=MUTED, va="center", fontweight="bold", transform=ax.transAxes)
    ax.text(0.05, y, link, fontsize=11.5, color=INK, va="center", transform=ax.transAxes)
    ax.text(0.05, y - 0.055, sym, fontsize=14, color=col, va="center", fontweight="bold", transform=ax.transAxes)
    ax.text(0.09, y - 0.055, ev, fontsize=10, color=INK, va="center", transform=ax.transAxes)
ax.set_title("Six questions we tested", loc="left", fontsize=12.5, color=INK)
ax.text(0.0, -0.02, "✔ supported   ~ supported but fragile   ✘ not supported (no association found, not proof of none)\n"
        "Questions 3–4: cell experiments; 5: patient tissue; 1, 2, 6: patient groups. FLT1 is measured as RNA, not sFlt-1 protein.",
        fontsize=8.5, color=MUTED, transform=ax.transAxes)
fig.tight_layout(); fig.savefig(FIG / "fig20_chain_evidence.png", dpi=200); plt.close(fig)
print("ok")
