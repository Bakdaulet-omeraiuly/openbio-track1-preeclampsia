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
| Sign-flip p (one-sided, 10,000 draws) | **≤ 1e-4** (the floor of 10,000 draws) |
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
  carry more FLT1, beyond technical covariation and beyond random gene sets.
  - The result **survives the purity check below**, on KRT7+ barcodes with no lineage marker.
  - This **supports link 5a in one cohort**, with a different design from O2's pseudobulk test: cell level and
    within donor.
- **PE vs control split (0.23 vs 0.045): not interpretable.**
  - FLT1 is much lower in control STB: in pure barcodes it is detected in a median of 16% of cells vs 37% in
    PE.
  - More zeros attenuate a Spearman correlation, so this split is not a disease difference. It is not carried
    into the hypothesis card.
- **R2:** the HIF2dep program tracks FLT1 beyond the generic hypoxia score within cells (20/20 donors). This
  matches bulk L-C4, but it is descriptive.

## Limits

- This is association across cells, not cause.
- **Dissociation stress is a real alternative, not just noise.** Dissociation can induce hypoxia and stress
  genes, and **FLT1 is itself hypoxia-inducible**. Cells stressed more during processing could therefore raise
  both, giving a within-cell link that does not exist in vivo. This dataset cannot rule that out. Single-nucleus
  data or spatial or in situ data from intact tissue would.
- **STB in scRNA-seq.** STB is a syncytium. "STB cells" here are barcodes with ≥ 3 STB markers, which may be
  fragments or STB-attached nuclei.
- **Effect size.** Per-cell ρ values are small (about 0.14). Single-cell data are sparse, so the
  donor-consistency (19/20) and the random-set margin carry the evidence, not the size.
- **Missing metadata.** Onset group and GA are not in the GEO metadata. The donor code FND was taken as "FD" (PE),
  and the result holds without it.
- **One cohort.** Replication in a second single-nucleus or single-cell PE dataset would be the next test.

## Post hoc purity check (`r3_purity_posthoc.py`, committed `d4d2c18` before running)

**Why.** A marker-only diagnostic, run before this check and without hypoxia or FLT1 values, showed the problem.
Of the R1-selected barcodes, a median of **33% also carry a non-trophoblast lineage marker** (PTPRC, CD14,
PECAM1, COL1A1, HBB), and only **65% are KRT7+**. The six STB markers are abundant ambient transcripts.

**Rule.**

- Barcodes: the R1 selection, **and** KRT7 detected, **and** none of the 5 lineage markers.
- Covariates: as in R1, plus a library indicator for donors with two libraries.
- The same gate, sign-flip test and random-set null as R1.

| | Donors (≥ 30 barcodes) | Mean ρ | ρ > 0 | Sign-flip p | Random-set 95th pct | Meets both |
|---|---|---|---|---|---|---|
| **R1, pure trophoblast barcodes** | 19 (FLR has 23 barcodes) | **0.099** | **17/19** | **0.0006** | −0.025 (median −0.084; observed above all) | **yes** |
| R2, pure (descriptive) | 19 | 0.217 | 18/19 | ≤ 1e-4 | 0.017 | — |

**Reading (pre-stated).** The within-donor link **stands for marker-pure trophoblast barcodes**. It is weaker
than in the marker-selected set (0.099 vs 0.137), so contamination inflated R1 somewhat, but it did not create it.

Pure barcodes range from 49 to 2,536 per donor. The two negative donors are FVB and FVT, both controls.
