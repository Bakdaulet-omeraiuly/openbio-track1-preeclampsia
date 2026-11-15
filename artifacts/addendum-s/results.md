# Addendum S results (2026-11-15)

The plan and script were committed and pushed (`5d385f3`) before any spatial file was downloaded.

- One change was made before the download: DCN, a Hallmark hypoxia gene, was replaced by LUM in the stroma
  markers. The script's own assert caught it, before any data existed.
- The script was run **once**. Outputs: `S_results.json`; input checksums: `inputs.sha256`.

## Gate: passed

All 19 donor codes have ≥ 200 in-tissue spots with ≥ 500 UMI: 1,005–4,999 spots per donor, 46,453 in total.

- Gene universe: 7,675 genes. Of these, 139 are hypoxia genes and 62 are HIF2dep genes.
- No marker gene is a hypoxia, HIF2dep or FLT1 gene (asserted).

## S1 (confirmatory): within-donor partial Spearman ρ(hypoxia score, FLT1 | spot composition + technical covariates)

**Covariates.**

- Composition: STB, endothelium, stroma, immune, erythrocyte marker scores.
- Technical: UMI, genes, mito, section.

| Statistic | Value |
|---|---|
| Mean ρ across 19 donors | **0.118** |
| Donors with ρ > 0 | **18 / 19** (negative: FGS) |
| Sign-flip p (one-sided) | **≤ 1e-4** (the floor of 10,000 draws) |
| 1,000 random gene sets | median 0.013, 95th percentile 0.058; observed at the **99.9th percentile** |
| **Verdict (both conditions)** | **Supported** |

## Pre-specified sensitivity analyses (S1)

| Analysis | Mean ρ | ρ > 0 | Sign-flip p | Random-set pct | Meets both |
|---|---|---|---|---|---|
| (a) without "ND" | 0.116 | 17/18 | 0.0002 | 100 | yes |
| (b) technical covariates only | 0.140 | 17/19 | ≤ 1e-4 | 99.9 | yes |
| (c) STB-rich half of spots | 0.078 | 18/19 | ≤ 1e-4 | 100 | yes |

Adjusting for composition lowers the estimate only modestly (0.140 to 0.118).

## S2 (descriptive): ρ(HIF2dep, FLT1 | hypoxia score + all covariates)

The mean ρ is 0.170, with **19/19** donors positive. Sign-flip p ≤ 1e-4, and the observed value is above all
1,000 random sets (95th percentile 0.041).

## Reading (by the rule written before the run)

**S1 is supported.** In intact placental tissue that was never dissociated, spots with a stronger hypoxia program
carry more FLT1 within the same placenta. This holds beyond spot composition, technical covariation and random
gene sets. **Dissociation stress cannot explain this**, which was the main open alternative to Addendum R.

Together:

- **R** (single cells, pure trophoblast barcodes) removes cell mix;
- **S** (intact tissue) removes dissociation.

They point the same way, with similar sizes: ρ ≈ 0.10 and 0.12.

## Limits

- This is an association, not cause.
- A Visium spot (55 µm) holds several cells. Composition is **adjusted for** with marker scores, not removed.
  The STB-rich half still shows the link (0.078, 18/19).
- R and S come from the **same 19–20 donors**, so they are two designs in one cohort, not two cohorts.
- Spatial autocorrelation makes neighbouring spots non-independent. Inference is therefore at the donor level
  (sign-flip), and the random-set null shares the same spot structure.
- Onset group and GA are not in the GEO metadata. "ND" is taken as FND; the result holds without it.
