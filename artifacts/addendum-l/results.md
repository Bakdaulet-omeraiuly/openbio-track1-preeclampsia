# Addendum L results (2026-11-15)

The plan, the screen and the confirmation script were each committed before the data they use were
opened (`fe4ee07`, `25892db`, `a9aa9d8`). The untouched cohorts were opened once (`ba57966`).

## Step 1: perturbation-defined gene set (GSE216799, descriptive, n = 3 placentas)

GSE216799 is primary term cytotrophoblasts at 2.5% O2, treated for 4 days with PT2385 (an HIF-2
inhibitor) or DMSO. The GEO record shows no linked paper. Its design matches the PT2385 study in
Hypertension 2023 ([10.1161/HYPERTENSIONAHA.122.20739](https://doi.org/10.1161/HYPERTENSIONAHA.122.20739)),
but we did not verify that link.

- **FLT1 was lower with PT2385:** log2FC −1.93, DESeq2 padj ≈ 10⁻⁷¹. This is a re-analysis of that
  experiment, in line with Sasagawa et al. 2018; it is not independent confirmation.
- **Gene set.** We defined it as the 100 genes most lowered by PT2385, after excluding Hallmark hypoxia
  genes and FLT1, and call it **"genes lower with PT2385"** (HIF2dep in the code). An inhibitor given
  during 4 days of differentiation can shift cell state, so these are not proven direct HIF-2 targets.
  Several are syncytiotrophoblast-associated genes (e.g. ERVH48-1, SIGLEC6).
- **The F1 (Blue-type) module score was lower with PT2385 in 3/3 donors** (−0.28, −0.50, −0.49 z).
  60% of F1 genes had a negative log2FC, against 51% of all genes.

## Step 2: screen in the used cohorts

The statistic is the partial Spearman correlation between the HIF2dep score and FLT1, adjusting for
`hyp_all`, within disease placentas.

| Cohort | rho |
|---|---:|
| GSE203507 | 0.15 |
| GSE114691 | 0.18 |
| GSE75010 | 0.50 |
| Pooled | 0.37 (p = 5 × 10⁻⁷) |

The screen passed. These cohorts have been heavily viewed, so this is not a finding.

## Step 3: confirmation in untouched cohorts (opened once)

The cohorts were GSE35574 (86 placentas after averaging 8 technical repeats), GSE24129 (24) and
GSE147776 (28). Across them there were 40 PE, 42 FGR-only and 56 control placentas.

| Test | Pooled result | Holm p | Outcome |
|---|---|---:|---|
| **L-C4** HIF2dep–FLT1 partial rho, PE + FGR-only | **0.39 [0.19, 0.57]**; per cohort 0.34 / 0.44 / 0.49 | **5 × 10⁻⁴** | **Supported** |
| **L-C3** hypoxia–FLT1 rho within PE | **0.62 [0.36, 0.79]**; per cohort 0.65 / 0.62 / 0.56 | **1 × 10⁻⁴** | **Supported** |
| L-C1 hypoxia-only, PE vs FGR-only (two-sided) | −0.19 SD [−0.56, 0.19] | 0.33 | No clear difference |
| L-C2 hypoxia-only, PE vs control | fixed 0.34 [−0.06, 0.75], τ² = 0.94; per cohort −0.61 / +0.93 / +1.14 | 0.09 | **Not supported** |

**Sensitivity.** Leave-one-cohort-out pooled C4 stayed between 0.36 and 0.47, with every p < 0.003.
C3 stayed between 0.58 and 0.65. C1 stayed non-significant. Dropping the labour and batch covariates
in GSE35574 did not change C1, C3 or C4.

## Post hoc checks (written after C4 was seen; they can only limit the claim)

**Random-gene-set calibration (`L4_calibration_posthoc.json`).** We drew 1,000 random 100-gene sets,
matched by expression and excluding hypoxia genes, FLT1 and HIF2dep, and applied the same pooled
partial correlation.

- Observed: 0.42 (79 genes measured in all three cohorts).
- Random sets: median −0.43, 95th percentile −0.14.
- HIF2dep ranked above all 1,000, so the empirical p was 0.001.
- In these arrays a random gene set tends to move *against* FLT1 once hypoxia is removed, so HIF2dep's
  positive link is clearly not generic co-expression.

**Overlap with the known modules.** HIF2dep shares:

- 8 of 100 genes with our F1 module;
- 19 with Ackerman VT Blue (1,001 genes);
- 13 with Ackerman VT Black (391 genes).

So C4 is not simply the Blue-type module again.

**Why C2 failed (post hoc reading, not tested).** GSE35574 controls are almost all term deliveries
(median 38 weeks; 1 of 40 delivered before 34 weeks), while its PE placentas have a median of 34 weeks.
The GA adjustment therefore extrapolates. The early-onset hypoxia signal may simply not be measurable
against term controls. This is an untested explanation.

**Processing disclosures.**

- The Ensembl-id-to-symbol mapping for GSE216799, from the Ensembl GRCh37.87 GTF, was added after the
  plan.
- All 8 GSE35574 technical repeats are IUGR arrays in batch B whose originals are in batch A. Each
  averaged placenta keeps the original's batch label (A).

## Conclusion (answering the brief)

- **What was tested.** A perturbation-defined trophoblast gene program, a pre-registered screen and a
  one-time confirmation in three untouched cohorts.
- **Supported:**
  - The genes lowered by HIF-2 inhibition in hypoxic trophoblasts track placental FLT1 beyond the generic
    hypoxia score, in three new cohorts and more strongly than matched random gene sets.
  - Hypoxia tracks FLT1 within PE in three new cohorts.
- **Not supported:**
  - The hypoxia-only score being higher in PE than in controls (C2) did not hold in these mostly term
    cohorts.
  - No PE-vs-FGR difference was found, now with 42 FGR-only placentas: −0.19 SD [−0.56, 0.19].
- **What this does not show.**
  - That HIF-2 directly causes FLT1 variation in vivo.
  - That the program is free of cell-composition effects; for example, it contains
    syncytiotrophoblast genes.
- **Next decisive test.**
  - Genetic HIF-2 (EPAS1) knockdown in trophoblasts, with RNA-seq, to separate direct HIF-2 targets from
    differentiation effects.
  - Repeating C4 with a syncytiotrophoblast-content covariate, pre-registered, in a new cohort.
