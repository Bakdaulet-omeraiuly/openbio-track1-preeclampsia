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

**STB-marker adjustment (`l5_stb_adjustment_posthoc.py`, `L5_stb_adjustment_posthoc.json`; committed
before running, `ba1135a`).** The covariate is the datathon marker matrix's syncytiotrophoblast class
(CGA, CGB, ERVW-1, PSG1). None of these are HIF2dep genes. CGB uses a fixed alias rule
(CGB/CGB3/CGB5/CGB8), and GSE35574 has no ERVW-1 probe, so only 3 to 4 genes are available. The
baseline reproduced L3 exactly (0.34 / 0.44 / 0.49, pooled 0.39). Adjusted pooled intervals use
n − 3 − k weights, where k counts covariates beyond hyp_all; the baseline keeps n − 3, as in L3. The EVT
score has 3 markers on GPL6244, which has no ITGA1.

| Adjustment (PE + FGR-only) | Pooled partial rho [95% CI] | Per cohort |
|---|---|---|
| hyp_all (L3 C4) | 0.39 [0.19, 0.57] | 0.34 / 0.44 / 0.49 |
| hyp_all + STB (primary) | 0.39 [0.18, 0.57] | 0.36 / 0.31 / 0.52 |
| hyp_all + STB + EVT (secondary) | 0.42 [0.20, 0.60] | 0.46 / 0.24 / 0.44 |

- **Calibration under the same adjustment:** observed 0.43 (79 common genes); random-set median −0.38;
  95th percentile −0.15; empirical p = 0.001.
  The random pool here also excludes the STB/EVT markers and CGB aliases. The same seed therefore
  gives slightly different draws from L4, despite the script docstring's "same 1,000 draws".
- **The STB score correlates with FLT1** (rho 0.30 to 0.43) **but inconsistently with HIF2dep** (−0.05,
  0.36, 0.35).
- **Where the check has leverage.** In GSE35574, 46 of the 82 disease placentas, the STB score is
  unrelated to HIF2dep (rho −0.05), so adjusting for it could not move C4 there. That may be a poor
  proxy on GPL6102, which has only CGA, CGB5/CGB8 and PSG1, or HIF2dep may not be STB-driven in bulk
  tissue; the data cannot tell these apart. In the two small cohorts where the STB score does track
  HIF2dep, C4 moved 0.44 → 0.31 (GSE24129, n = 16) and 0.49 → 0.52 (GSE147776, n = 20).
  "Unchanged" therefore describes the pooled value, not every cohort.
- **Reading, by the rule fixed before the run:** C4 is *not explained by this STB marker score*. That
  is not the same as "free of composition", because 3 to 4 markers are a coarse proxy. STB content
  could also be a mediator rather than a confounder, since HIF-2 is reported to impair
  syncytialisation. This check limits the claim less than feared, but it does not upgrade it.

**Single-cell-reference deconvolution (`l6_deconvolution_posthoc.py`, `L6_deconvolution_posthoc.json`,
`L6_estimated_proportions_*.csv`; committed before running, `abdbe7b`).**

- **Reference:** the datathon single-cell pseudobulk, 14 cell classes.
- **Method:** the existing Addendum B5 NNLS rule, giving 700 signature genes. Signature genes were
  limited to genes on all three platforms, with HIF2dep, FLT1 and hypoxia genes removed.
- **Not the pre-registered new-cohort test:** these are the same three cohorts as C4.
- **Baseline:** reproduced L3.

| Adjustment (PE + FGR-only) | Pooled partial rho [95% CI] | Per cohort (GSE35574 / GSE24129 / GSE147776) |
|---|---|---|
| hyp_all (L3 C4) | 0.39 [0.19, 0.57] | 0.34 / 0.44 / 0.49 |
| hyp_all + STB + CTB proportions (primary) | **0.43 [0.21, 0.60]** | 0.40 / 0.37 / 0.52 |
| hyp_all + STB + CTB + EVT proportions (secondary) | 0.34 [0.11, 0.53] | 0.40 / **−0.15** / 0.45 |

- **Calibration under the primary adjustment:** observed 0.45 (79 common genes); random-set median
  −0.33; 95th percentile −0.11; empirical p = 0.001.
- **Reading, by the rule fixed before the run:** C4 is *not explained by deconvolution-estimated
  STB/CTB proportions*.
- **The STB axis is plausible.** The estimated STB proportion agrees with the STB marker score
  (rho 0.68 / 0.78 / 0.97).
- **The absolute proportions are not credible.** For example, the GSE147776 median Hofbauer fraction
  is 0.35 and its CTB fraction is near 0, while the GSE35574 CTB fraction is 0.58. NNLS on log-scale
  arrays turned back to linear values is rough, so only the rank of the STB estimate should be
  trusted.
- **A real caution.** GSE24129 is the cohort where estimated STB tracks HIF2dep (rho 0.55). There,
  adding the EVT proportion as well turned the per-cohort link negative (−0.15). That cohort has only
  16 placentas and 4 covariates, so the estimate is unstable, but it is exactly where composition
  had leverage. Where composition varies with the program, we cannot rule it out.
- **Overall.** Both composition checks (L5, L6) leave the pooled C4 link positive and above random
  sets. Neither is decisive. The decisive test is still a new, pre-registered cohort with a better
  composition estimate, or genetic EPAS1 knockdown.

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
  - That the program is free of cell-composition effects. It contains syncytiotrophoblast genes. A
    post hoc 3–4-marker STB adjustment (0.39 [0.18, 0.57]) and a single-cell-reference
    deconvolution adjustment (0.43 [0.21, 0.60]) did not weaken pooled C4. Both are rough, and in
    GSE24129, where estimated STB tracks the program, the full-composition model turned that
    cohort's estimate negative (n = 16).
- **Next decisive test.**
  - Genetic HIF-2 (EPAS1) knockdown in trophoblasts, with RNA-seq, to separate direct HIF-2 targets from
    differentiation effects.
  - Repeating C4 in a new cohort, pre-registered, with a better composition covariate, e.g.
    deconvolution against a placental single-cell reference. The post hoc marker check above is not
    that test.


## Post hoc l8: does diagnosis explain L-C4? (`l8_diagnosis_sensitivity_posthoc.py`, committed before running; reviewer request)

| Sensitivity | Pooled partial ρ [95% CI] | p (one-sided) | Per cohort (GSE35574 / GSE24129 / GSE147776) |
|---|---|---|---|
| Baseline L-C4 (PE + FGR-only) | 0.39 [0.19, 0.57] | 0.0002 | 0.34 / 0.44 / 0.49 |
| (a) + PE indicator (diagnosis-adjusted) | **0.35 [0.13, 0.53]** | 0.001 | 0.32 / 0.14 / 0.53 |
| (b) PE only (n = 19 / 8 / 13) | **0.36 [0.03, 0.62]** | 0.018 | 0.36 / 0.10 / 0.48 |
| (c) FGR-only only (n = 27 / 8 / 7) | 0.31 [−0.02, 0.58] | 0.031 | — |
| (d) GSE35574 + PE, GA, sex, labour, batch | 0.20 (from 0.34 on the same set) | — | — |

**Reading.** The link is not explained by diagnosis-group differences: it holds when diagnosis is adjusted for,
and within PE alone. It is attenuated, though still positive, when GA, sex, labour and batch are added in the one
cohort that records them. This is post hoc and can only limit the claim.
