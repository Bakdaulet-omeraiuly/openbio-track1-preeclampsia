# Addendum S: does the within-donor hypoxia–FLT1 link exist in intact, undissociated tissue? (GSE319236 Visium)

Written on 2026-11-15, before any spatial expression file was downloaded or read. The plan, `s_spatial.py` and
`sections.tsv` are committed together. The data are then downloaded and the script is run **once**.

## Why

Addendum R found a within-donor, within-trophoblast link between the hypoxia score and FLT1 in dissociated
scRNA-seq: ρ = 0.10 in pure barcodes, 17/19 donors. Its main open alternative is **dissociation stress**.
Processing can induce hypoxia and stress genes, and FLT1 is itself hypoxia-inducible.

Visium spatial transcriptomics is taken from frozen tissue sections that were **never dissociated**. If the link
also exists across spots within a donor, after adjusting for each spot's cell composition, dissociation cannot
explain it.

## Data (viewed so far: GEO sample metadata only)

- GSE319236 placenta "PVBP-Spatial": 24 sections from 19 donor codes (`sections.tsv`).
  - Two donors have two sections each.
  - FJJ has no section.
  - The code "ND" is taken as donor "FND" ("FD", PE in the series design). An analysis without it is
    pre-specified as a sensitivity check.
- Groups are as in Addendum R.
- No spatial file has been used anywhere in this repository.

## Spots (fixed)

- In-tissue spots (`tissue_positions` in_tissue = 1) with ≥ 500 UMI.
- Sections of the same donor are pooled, with a section indicator as a covariate.

**Gate:** the test runs only if ≥ 8 donors have ≥ 200 spots.

## Scores, per donor and across its spots

- **Normalisation:** log1p(count / total × 10⁴).
- **Gene universe:** genes detected in ≥ 5% of pooled spots, with FLT1 excluded.
- **`hyp_all`:** mean z of the Hallmark hypoxia genes in the universe, z-scored across the donor's spots.
- **`hif2dep`:** the same, using `artifacts/addendum-l/L1_HIF2dep_genes.txt`.
- **FLT1:** log-normalised value.
- **Composition covariates:** mean z of each marker set below. No marker is a hypoxia, HIF2dep or FLT1 gene; the
  script asserts this.
  - **STB:** CGA, CYP19A1, PSG1, PSG4, KISS1, LGALS13.
  - **Endothelium:** PECAM1, VWF, CDH5. FLT1 is also endothelial, so this covariate is essential.
  - **Stroma:** COL1A1, COL3A1, LUM. DCN was replaced before the data were opened, because it is a Hallmark hypoxia gene.
  - **Immune:** PTPRC, CD14, CD68.
  - **Erythrocytes:** HBB, HBA1.
- **Technical covariates:** log10 UMI, log10 detected genes, mitochondrial fraction, section indicator.

## Tests (the donor is the unit)

**S1 (confirmatory).** Within-donor partial Spearman ρ(`hyp_all`, FLT1 | composition + technical covariates).

- The covariates are removed by linear residuals of both variables; then Spearman is computed.
- **Supported only if both:**
  - (a) the one-sided sign-flip p across donors is < 0.05; and
  - (b) the mean ρ is above the 95th percentile of 1,000 decile-matched random gene sets.
    - The sets have the same size and are drawn without hypoxia, HIF2dep, FLT1 or marker genes. Seed 0.
    - With spots, (b) matters even more: neighbouring spots share technical and tissue-level variation.

**S2 (descriptive).** ρ(`hif2dep`, FLT1 | `hyp_all` + all covariates), with its sign-flip p and random-set
percentile.

**Sensitivity analyses (pre-specified):**

- (a) without "ND";
- (b) without the composition covariates (technical covariates only), to show how much they matter;
- (c) restricted to STB-rich spots, i.e. the top half of the donor's STB score.

## Pre-written readings

- **S1 supported.** "In intact, undissociated placental tissue, spots with a stronger hypoxia program carry more
  FLT1 within the same placenta, beyond spot composition, technical covariation and random gene sets.
  Dissociation stress cannot explain this."
  - It is still an association.
  - A spot is not a single cell, so composition is adjusted for, not removed.
- **S1 fails.** "No composition-adjusted within-donor link was found in intact tissue." The R result may then
  partly reflect dissociation, or the Visium resolution and sensitivity may be too low. We cannot tell which.
- **Gate fails.** "Not testable."
