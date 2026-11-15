# Addendum O1 results: ARNT knockdown and the PT2385-defined program (GSE80995)

The plan and script were committed before any expression value was read (`a823f57`). The script ran once.

- **Plan time correction:** `plan.md` says "13:05 EST". The plan was actually written around 12:45 and
  committed at 12:46 EST. The order (commit, then run) is unchanged.
- **Disclosure:** the plan assumed the series values were already log2. The shared reader applied log2
  because the maximum value was > 100. This is reported, not changed.

**Source.** Wakeland et al. 2017, Am J Pathol, [10.1016/j.ajpath.2016.11.018](https://doi.org/10.1016/j.ajpath.2016.11.018).
The study used primary first-trimester villous cytotrophoblasts from 5 placentas, with shARNT or
scrambled shRNA, at 2% or 20% O2. Δ means shARNT minus scramble within a donor, in gene-z units.

## Pre-specified

| Test | Result | Verdict |
|---|---|---|
| Gate: ARNT Δ at 2% O2 | −1.98 (5/5 donors lower) | passed |
| Positive control: hyp_all Δ at 2% O2 | −0.46 (5/5 lower, sign-flip p = 0.031) | passed |
| **O1 primary: HIF2dep Δ at 2% O2** | **−1.15 (5/5 donors lower, −0.82 to −1.30); sign-flip p = 0.031, the minimum possible with 5 donors** | **supported** |
| O1 calibration: 1,000 matched random sets | random median +0.05, 5th percentile −0.09; HIF2dep below all 1,000 | supported |
| O1 beyond hyp_all (descriptive) | HIF2dep Δ − hyp_all Δ = −0.69 | the program falls more than the generic hypoxia genes |
| FLT1 Δ at 2% O2 (descriptive) | −0.35 (2 donors up, 2 down, 1 flat) | **not consistent** |

## Amendment, after the run: a sample-label conflict

- **The conflict.** In donor 630P, the two normoxia arrays (GSM2139930 and GSM2139932; an earlier version of this file mistyped them as GSM2139927/GSM2139930. GSM2139927 is the 630P hypoxia scramble array, whose labels agree) have titles
  (Scramble / shARNT) that contradict their `treatment` field. The script used `treatment`.
- **ARNT expression supports the titles, but does not prove them.** We have no confirmation from the authors, so both versions are kept and the title-label numbers are a sensitivity analysis under the alternative labelling. Under the `treatment` labels, 630P's "shARNT" array had
  *higher* ARNT (Δ +1.82), while every other donor had a lower one (−1.3 to −2.1).
- **The rerun** (`o1_amendment_title_labels.py` → `O1_results_title_labels.json`) uses the titles.
- **Effect.** The hypoxia pairs are identical under both rules, so the primary result, the gate, the
  positive control and the calibration are unchanged. Only the normoxia numbers change:

| | `treatment` labels (as run, primary) | title labels (sensitivity; supported by ARNT expression) |
|---|---|---|
| HIF2dep Δ at 20% O2 | −0.52 | **−0.85** |
| Specificity (Δ at 2% − Δ at 20%) | −0.63 | **−0.30** |

**Reading:** the program falls with ARNT knockdown at both oxygen levels, and only a little more at 2% O2.
Much of its ARNT dependence is therefore not oxygen-driven. Possible sources are basal HIF activity, AhR
(another ARNT partner) or differentiation.

## Post hoc, descriptive (`o1_posthoc_hypoxia_induction.py`)

In the scramble arm, 2% vs 20% O2 changed the scores as follows:

| Score | Change at 2% vs 20% O2 | Per donor |
|---|---|---|
| HIF2dep | +0.88 | rises in 5/5 donors |
| hyp_all | +0.22 | rises in 5/5 donors |
| FLT1 | **−1.31** | falls in 4/5 donors |

In these first-trimester cells, low oxygen does not raise FLT1. FLT1 itself was measured after ARNT knockdown, and the change was inconsistent across donors. So this system gave **no support for a consistent FLT1 fall after ARNT knockdown**. That limits the conclusion; the low-oxygen pattern above is a possible reason.

## What O1 shows and does not show

- **Shows.** The genes lowered by the drug PT2385 in term cytotrophoblasts were also lowered when ARNT
  (HIF-1β) expression was reduced by shRNA in first-trimester cytotrophoblasts. This is a different lab, a
  different perturbation and a different gestational stage. The drop was in 5/5 donors and beyond all random
  gene sets. This is **agreement between two independent perturbations**.
- **Not shown.** That every drug off-target effect is excluded. ARNT has other partners (HIF-1α, AhR), and
  both perturbations also change differentiation. The pre-written reading label in `o1_arnt.py` ("depends on
  the HIF complex genetically") is stronger than this, and we do not use it as the conclusion.
- **Does not show:**
  - HIF-2α specificity: EPAS1 mRNA rose with the knockdown (+1.25).
  - That the dependence is oxygen-specific (see the sensitivity analysis).
  - A consistent FLT1 fall after ARNT knockdown. FLT1 was measured, but the change was inconsistent across
    donors, so this result limits the conclusion.
  - Independence from cell state: Wakeland et al. link HIF to trophoblast differentiation.
- **Method note.** The knockdown is shRNA-mediated, not a gene knockout.

## Companion test O2, by the second session (`artifacts/addendum-o2/`)

- **Primary, not supported.** Within syncytiotrophoblast alone (26 donors), the partial rho of HIF2dep with
  FLT1 | hyp_all was 0.26 [−0.14, 0.59], with permutation p = 0.11.
- **Reading.** The direction is the same, but we cannot claim that L-C4 is free of cell-mix effects.
