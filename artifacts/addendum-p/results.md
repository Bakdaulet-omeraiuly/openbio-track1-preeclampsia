# Addendum P results: GSE186257 (fourth untouched cohort for the HIF2dep–FLT1 link)

The plan and script were committed before any expression value was read (`d2be780`). The script ran once.
Full output: `P_results.json`. Source: Horii et al. 2023, [10.1016/j.modpat.2022.100035](https://doi.org/10.1016/j.modpat.2022.100035).

## Pre-specified

| Test | Result | Verdict |
|---|---|---|
| **P1:** HIF2dep–FLT1 partial rho \| hyp_all, PE + SGA-only (n = 34) | **0.40 [0.06, 0.65]**, permutation p = 0.011, Holm p = 0.011 | **replicates L-C4** |
| **P2:** hyp_all–FLT1 within PE (n = 26) | **0.69 [0.42, 0.85]**, Holm p = 0.0002 | **replicates** |
| P1-cs: P1 additionally adjusted for STB and EVT marker scores | **0.63 [0.35, 0.80]** | not explained by these marker scores (coarse proxies) |
| P1 within PE only (n = 26) | 0.26 [−0.15, 0.59], p = 0.10 | same direction, not clear |
| Descriptive: 1,000 matched random sets for P1 (unadjusted) | random median **+0.23**, 95th percentile 0.48; HIF2dep at the **83rd percentile** | **specificity not shown** in this cohort |
| Descriptive: hyp_only, PE vs non-SGA controls | +0.06 [−0.16, 0.28] | no difference (no GA; controls probably later) |
| Descriptive: hyp_only, PE vs SGA-only | +0.12 [−0.11, 0.36] | no clear difference |

- **Descriptive correlations.** In the disease group, HIF2dep correlates with the EVT marker score (rho
  0.54) more than with the STB score (0.18).
- **Random sets behave differently here.** Unlike the three array cohorts of Addendum L, where the median
  was −0.43, random gene sets in this RNA-seq cohort tend to move *with* FLT1 after hypoxia adjustment.
  The unadjusted P1 therefore does not show that HIF2dep is special here.

## Post hoc (`p_posthoc_null_cs.py`, written after P was seen; descriptive, it does not upgrade the claim)

| Adjustment | HIF2dep | Random median | Random 95th percentile | Percentile |
|---|---|---|---|---|
| hyp_all + STB + EVT (as P1-cs) | 0.63 | 0.00 | 0.32 | **100** |
| PE only, \| hyp_all | 0.26 | 0.41 | 0.67 | 27 |

**Reading.** Once the marker scores are removed, random sets no longer track FLT1, but HIF2dep does. That
fits the idea that the unadjusted random-set co-movement in this cohort reflects cell mix. This is post hoc
and was not planned. Within PE alone, HIF2dep is not above random sets.

## Summary

**Supported by the pre-specified tests:**

- The HIF2dep–FLT1 link (L-C4) replicates in a fourth untouched cohort, on another platform (RNA-seq) and
  from another lab.
- The link stays positive after adjusting for STB and EVT marker scores.
- The hypoxia–FLT1 link within PE replicates.

**Not shown:**

- Specificity against random gene sets without composition adjustment (83rd percentile).
- A clear link within PE alone.
- Any PE vs control or PE vs SGA difference in hypoxia score in this cohort, which lacks GA.
