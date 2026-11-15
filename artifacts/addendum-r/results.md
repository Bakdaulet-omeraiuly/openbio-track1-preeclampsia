# Addendum R results (2026-11-15)

Plan and script were committed and pushed (`df0298c`) before any GSE319236 expression file was downloaded. The
script was run **once** with no changes. Outputs: `R_results.json`, input checksums `inputs.sha256`.

## Feasibility gate: passed

All 20 donors had ≥ 30 STB cells (≥ 3 of 6 STB markers detected): range 113–4,445 cells per donor. Gene
universe: 7,146 genes detected in ≥ 5% of STB cells. Of these, 105 are Hallmark hypoxia genes and 28 are HIF2dep
genes.

## R1 (confirmatory): within-donor, within-STB partial Spearman ρ(hypoxia score, FLT1 | UMI, genes, mito)

| Statistic | Value |
|---|---|
| Mean ρ across 20 donors | **0.137** |
| Donors with ρ > 0 | **19 / 20** (the one negative: FLJ, −0.01) |
| Sign-flip p (one-sided, 10,000 draws) | **0.0001** |
| 1,000 random gene sets (decile-matched) | median −0.036, 95th percentile 0.038; observed above **all 1,000** |
| **Verdict (both pre-specified conditions)** | **Supported** |

## Pre-specified sensitivity analyses (R1)

| Analysis | Mean ρ | ρ > 0 | Sign-flip p | Random-set 95th pct | Meets both conditions |
|---|---|---|---|---|---|
| (a) without FND | 0.125 | 18/19 | 0.0001 | 0.039 | yes |
| (b) PE donors only | **0.234** | 10/10 | 0.001 | 0.053 | yes |
| (b) control donors only | 0.045 | 8/10 | 0.011 | 0.054 | **no** (92nd percentile) |
| (c) STB = ≥ 4 of 6 markers | 0.155 | 20/20 | 0.0001 | 0.021 | yes |

## R2 (descriptive): ρ(HIF2dep, FLT1 | hypoxia score + technical covariates)

The mean ρ is 0.240, with **20/20** donors positive. Sign-flip p is 0.0001, and the observed value is above all
1,000 random sets (95th percentile 0.064).

## R3 (descriptive): PE vs control, donor means in STB

- FLT1 is higher in PE STB: median 0.86 vs 0.28 log-normalised, Mann–Whitney p = 0.0006. This is a sanity check
  that the data carry the known PE signal.
- The mean hypoxia-gene level does not differ: 0.47 vs 0.44, p = 0.19.

## Post hoc (written and committed after R1/R2 were seen, before running; `r2_stb_state_posthoc.py`, `f948085`)

This check adds the STB marker score (mean z of the 6 selection markers) and the number of markers detected as
covariates, to test whether "more STB-like" cells simply carry more of both.

| | Mean ρ | ρ > 0 | Sign-flip p | Random-set 95th pct |
|---|---|---|---|---|
| R1 + STB state | 0.135 | 18/20 | 0.0001 | 0.019 (observed above all) |
| R2 + STB state | 0.212 | 20/20 | 0.0001 | 0.039 (observed above all) |

The link is **not explained by this STB-identity score.**

## Reading (by the rules written before the run)

- **R1 is supported.** Within individual placentas, syncytiotrophoblast cells with a stronger hypoxia program
  carry more FLT1, beyond technical covariation and beyond random gene sets. Comparing cells of one type inside
  one person removes cell mix. This answers link 5a, which had not been established in O2's pseudobulk test,
  with a different design: cell level and within donor.
- **The coupling is stronger in PE than in controls.**
  - PE donors: 0.23, 10/10.
  - Control donors: 0.045, not beyond random sets.
  - This is a descriptive split, not a pre-specified contrast.
- **R2:** the HIF2dep program tracks FLT1 beyond the generic hypoxia score within cells (20/20 donors). This
  matches bulk L-C4, but it is descriptive.

## Limits

- This is association across cells, not cause.
- **Dissociation.** These are dissociated scRNA-seq cells, and dissociation stress can induce hypoxia-like genes.
  That would add noise, but it is not expected to track FLT1 specifically.
- **STB in scRNA-seq.** STB is a syncytium. "STB cells" here are barcodes with ≥ 3 STB markers, which may be
  fragments or STB-attached nuclei.
- **Effect size.** Per-cell ρ values are small (about 0.14). Single-cell data are sparse, so the
  donor-consistency (19/20) and the random-set margin carry the evidence, not the size.
- **Missing metadata.** Onset group and GA are not in the GEO metadata. The donor code FND was taken as "FD" (PE),
  and the result holds without it.
- **One cohort.** Replication in a second single-nucleus or single-cell PE dataset would be the next test.
