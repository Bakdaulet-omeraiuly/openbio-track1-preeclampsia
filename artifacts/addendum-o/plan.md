# Addendum O1: does genetic loss of the HIF complex lower the PT2385-defined program? (written before running)

Written on 2026-11-15, 13:05 EST, before any expression value of GSE80995 was read. Only the metadata lines
(sample titles, treatment, oxygen, platform, processing note) were viewed. The plan is committed together
with `o1_arnt.py`.

## Why

Addendum L defined "HIF2dep" as the genes lowered by a drug (PT2385) in term cytotrophoblasts. A drug can act
off-target, and a 4-day treatment changes differentiation. This test asks whether the same genes also fall
when the HIF complex is removed genetically, in another lab and another system. That is a step toward
showing that the program depends on a specific regulator, not only that it correlates with FLT1.

## Data

[GSE80995](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE80995) (Wakeland et al. 2017, PMID 28167044):

- **Cells:** primary first-trimester villous cytotrophoblasts from 5 placentas.
- **Treatment:** shARNT (HIF-1β, the obligatory partner of HIF-1α and HIF-2α) or scrambled shRNA.
- **Oxygen:** cultured at 2% or 20% O2.
- **Platform:** Illumina HT-12 v4 (GPL10558), lumi RSN-normalised (log2 assumed; the script checks the
  range).
- **Design:** 20 cultured arrays (5 donors × 2 shRNA × 2 O2). The 4 freshly isolated arrays are not used.

## Scores

- **Gene values:** probes are mapped to `Symbol` (GPL10558) and averaged per gene. Each gene is then
  z-scored across the 20 cultured arrays.
- **Score sets:**
  - HIF2dep: `artifacts/addendum-l/L1_HIF2dep_genes.txt`, frozen.
  - hyp_all: HALLMARK_HYPOXIA.
  - FLT1, ARNT and EPAS1: single genes.

## Tests

The knockdown effect is Δ = shARNT − scramble, computed within donor.

| Test | Contrast | Prediction | Statistic |
|---|---|---|---|
| **Gate (must pass for O1 to be read)** | ARNT gene Δ at 2% O2 | < 0 | knockdown happened |
| **Positive control** | hyp_all Δ at 2% O2 | < 0 | the HIF response is reduced |
| **O1 primary** | HIF2dep Δ at 2% O2 | < 0 | exact sign-flip permutation over the 5 donors, one-sided; the minimum possible p is 1/32 = 0.031 |
| O1 calibration | the same Δ for 1,000 random 100-gene sets, decile-matched on mean expression, excluding HYP, HIF2dep and FLT1 | HIF2dep below the random 5th percentile | seed 0 |
| O1 specificity (descriptive) | HIF2dep Δ at 2% O2 minus Δ at 20% O2 | < 0 | — |
| O1 beyond hypoxia genes (descriptive) | HIF2dep Δ minus hyp_all Δ | — | reported |
| FLT1 (descriptive) | FLT1 Δ at 2% O2 | < 0 | — |

## Reading rule

- **Gate or positive control fails:** O1 is uninterpretable and is reported as such.
- **O1 Δ < 0, permutation p ≤ 0.05 and below the random 5th percentile:** "the PT2385-defined program
  also depends on the HIF complex genetically".
- **Otherwise:** "not shown".

## Limits, stated in advance

- **ARNT is not HIF-2-specific.** It partners HIF-1α and AhR too, so O1 tests HIF-complex dependence, not
  HIF-2α specifically.
- **The ARNT effect may run through differentiation.** Wakeland et al. report that ARNT knockdown attenuates
  hypoxia-driven EVT differentiation, so a fall may partly reflect a change in differentiation.
- **Different system:** first-trimester CTB rather than term CTB.
- **Small n:** only 5 donors.
