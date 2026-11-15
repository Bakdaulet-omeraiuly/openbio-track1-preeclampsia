# Addendum O2 results (2026-11-15)

The plan and script were committed (`8abfcc1`) and run **once** with no changes. Output: `O2_results.json`.

## Primary (confirmatory)

| Test | n | Partial ρ (approx. 95% CI) | Permutation p (one-sided) | Verdict |
|---|---:|---|---:|---|
| STB, all donors: ρ(HIF2dep, FLT1 \| hyp_all) | 26 | **0.26 [−0.14, 0.59]** | **0.11** | **Not supported** |

The 26 STB donors are 10 early PE, 7 late PE, 3 early controls and 6 late controls. 67 HIF2dep genes and 146
hypoxia genes passed the expression filter.

## Pre-specified secondaries (descriptive, no verdict)

| Analysis | n | Partial ρ (approx. 95% CI) |
|---|---:|---|
| S1 STB + GA, PE status, onset group | 26 | 0.32 [−0.11, 0.65] |
| S2 villous cytotrophoblast | 26 | 0.57 [0.23, 0.78] |
| S2 extravillous trophoblast | 15 | **−0.54 [−0.83, −0.04]** |
| S3 STB random-set null | — | median of random sets −0.35, 95th percentile −0.14; observed 0.26 is above all 1,000 (empirical p = 0.001) |
| S4 STB, PE donors only | 17 | 0.21 [−0.30, 0.63] |

## Reading (by the rules written before the run)

- **Primary: positive but not supported.** The direction is the same as in L-C4, but the result is not
  significant with 26 donors. This is not evidence against L-C4, and it is not a confirmation either.
- **We cannot claim that L-C4 is free of cell composition.** The within-STB estimate (0.26) is lower than the
  bulk estimate (0.39) and its CI includes zero.
- **S3.** HIF2dep ranks above all random gene sets, so the STB estimate is not an artefact of gene-set size or
  expression level.
  - The random-set null is **centred well below zero** (median −0.35), as in L4. A random gene-set score that is
    residualised on the hypoxia score tends to move against FLT1.
  - So "above all random sets" and "about zero against a permutation null" are both true. The permutation p is the
    pre-specified decision.
- **S2 differs by cell type.** Villous cytotrophoblast is positive (0.57), while EVT is **negative** (−0.54,
  n = 15).
  - The HIF2dep program was defined in primary cytotrophoblasts (GSE216799), so a positive VCT result is
    expected. That makes it descriptive context, not a new finding.
  - The EVT sign is the opposite, and we do not interpret it further at this n.

## What this means for the claim

L-C4 remains a **bulk-tissue co-expression result in three untouched cohorts**. This test could not show that it
holds inside a single cell type with confidence. The direction is positive in STB and VCT and negative in EVT.
A larger single-cell cohort, or a pre-registered deconvolution-adjusted bulk test, is still needed.

## Disclosures

- These pseudobulk data were used before in Addenda A/D/F and L6, but never for HIF2dep per donor (see plan.md).
- The approximate CIs use Fisher z with n − 3 − k. They are descriptive; the permutation p decides.
- S4 re-applies the gene filter and z-scoring within the 17 PE donors.
