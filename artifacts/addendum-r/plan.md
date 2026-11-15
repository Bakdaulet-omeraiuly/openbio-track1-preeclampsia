# Addendum R: is the hypoxia–FLT1 link present *within* syncytiotrophoblast cells of one donor? (GSE319236)

Written on 2026-11-15, before any expression value of GSE319236 was downloaded or read. The plan,
`r_within_donor.py` and `libraries.tsv` are committed together. The data are then downloaded and the script is
run **once**.

## Why

Link 5 in the hypothesis card, hypoxia ↔ FLT1 in placental tissue, is positive in 10 of 10 bulk cohorts. Link 5a
asks whether this holds within one cell type, rather than reflecting the cell mix. It is not established: within
syncytiotrophoblast (STB) pseudobulk across 26 donors (O2), ρ = 0.26 [−0.14, 0.59].

A **within-donor, cell-level** test removes between-donor and between-cell-type composition entirely. Within one
person, do STB cells with a higher hypoxia score also carry more FLT1?

## Data (viewed so far: GEO sample metadata only)

- **GSE319236**, placental single-cell RNA-seq, 10x 3′ v3, CellRanger filtered matrices.
- 25 libraries ("PVBP", not "Spatial") from 20 donors, listed in `libraries.tsv`.
- The GEO series design lists the groups:
  - **PE:** FCM FAM FGS FEP FVS FFM FRH FSH FD FMM;
  - **control:** FVQ FLR FJJ FVB FLJ FJD FVT FYI FRK FNS.
  - The library donor code "FND" is not in either list. We take it to be "FD" (PE), and an analysis without it is
    pre-specified as a sensitivity check.
- Early/late onset and GA are not in the GEO metadata.
- This dataset has not been used anywhere in this repository.

## Cell selection (fixed, no expression-driven tuning)

- **QC:** ≥ 500 detected genes and mitochondrial fraction < 20% (genes starting with "MT-").
- **STB marker set:** CGA, CYP19A1, PSG1, PSG4, KISS1, LGALS13.
- **STB cell:** at least **3 of the 6** markers detected (count > 0).
- Libraries of the same donor are pooled.

## Feasibility gate (label-free, checked before any hypoxia or FLT1 value is computed)

The test runs only if **≥ 8 donors have ≥ 30 STB cells**. Otherwise the result is "not testable in this dataset",
and nothing else is computed. Dissociation-based scRNA-seq is known to capture STB poorly.

## Scores, per donor and within its STB cells

- **Normalisation:** log1p(count / total × 10⁴).
- **`hyp_all`:** mean z of Hallmark hypoxia genes (v2024.1). Genes are z-scored across the donor's STB cells, and
  zero-variance genes are dropped. FLT1 is not in the set (asserted).
- **`hif2dep`:** the same, using `artifacts/addendum-l/L1_HIF2dep_genes.txt`.
- **FLT1:** log-normalised value.
- **Technical covariates:** log10 total UMI, log10 detected genes, mitochondrial fraction.

## Tests (the donor is the unit)

**R1 (confirmatory).** Within-donor partial Spearman ρ(`hyp_all`, FLT1 | technical covariates).

- The covariates are removed by linear residuals of both variables; then Spearman is computed.
- Statistic: mean ρ across eligible donors.
- **Supported only if both:**
  - (a) the one-sided sign-flip permutation p across donors is < 0.05 (exact if ≤ 2^16 sign patterns, otherwise
    10,000 draws, seed 42); and
  - (b) the observed mean ρ is above the 95th percentile of **1,000 random gene sets**.
    - The random sets have the same size as the measured hypoxia set, drawn within pooled-STB mean-expression
      deciles matched to the hypoxia set.
    - Excluded from the draws: hypoxia genes, HIF2dep genes and FLT1. Seed 0.
    - Condition (b) protects against technical co-variation: in scRNA-seq, any gene set tends to correlate with
      any gene.

**R2 (descriptive).** Within-donor partial ρ(`hif2dep`, FLT1 | `hyp_all` + technical covariates), as in L-C4.
Reported with its sign-flip p and its random-set percentile, without a verdict.

**R3 (descriptive).** PE vs control in the donor-level mean STB `hyp_all` and mean FLT1 (Mann–Whitney,
two-sided).

**Sensitivity analyses (pre-specified):**

- (a) R1 without donor FND;
- (b) R1 in PE donors only and in control donors only;
- (c) R1 with an STB threshold of ≥ 4 of 6 markers.

## Pre-written readings

- **R1 supported.** "Within individual placentas, STB cells with a stronger hypoxia program carry more FLT1,
  beyond technical covariation and beyond random gene sets. This does not depend on cell mix."
  - This is still association within cells, not cause.
  - The cells are dissociated, so stress during dissociation could induce hypoxia-like genes. That would add
    noise, though it is not obviously expected to track FLT1.
- **R1 fails (a) or (b).** "No within-STB, within-donor link beyond technical covariation was found in this
  dataset." This weakens link 5a further.
- **Gate fails.** "STB is too sparse in this dissociated scRNA-seq dataset to test link 5a." Not evidence either
  way.
