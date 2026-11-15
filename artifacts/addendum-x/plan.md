# Addendum X: does placental infarction track the HIF-2/FLT1 axis in an independent cohort? (ELGAN, GSE154829)

Written on 2026-11-15, before any expression value of GSE154829 was downloaded. The plan, `x_infarct.py` and
`samples.tsv` (from GEO metadata only) are committed together, then run **once**.

## Why

In Addendum W (GSE303463), histological maternal vascular malperfusion (MVM) did **not** raise the generic hypoxia
score (W1, pre-registered). A post hoc calibration found, however, that MVM shifted the **HIF2dep program**
(above all 1,000 random sets) and **FLT1** (98.6th percentile), while generic hypoxia stayed at the 7.9th
percentile.

That post hoc lead is tested here, pre-registered, in an independent cohort with a malperfusion lesion recorded:
**placental infarction**, a core MVM lesion.

## Data

- **GSE154829** (ELGAN; PMIDs 32721520, 33308293): 785 placentas from births at 23–27 weeks, RNA-seq counts,
  14 US centres. Independent of GSE303463 and of every cohort used so far.
- **Exposure:** infarct (1 = present, 127; 0 = absent, 616; NA = 42, excluded).
- **Covariates:** fetal sex, chorioamnionitis, fetal vessel neutrophilic infiltration, smoking exposure, maternal
  age, maternal BMI.
  - **Not adjusted:** birthweight and placental weight, which may be consequences of infarction.
  - GA is not available per sample; the cohort's GA range is 23–27 weeks.
  - Preeclampsia is not recorded.
- **Primary analysis set:** singletons (`fused2` = 0): 92 infarct vs 359 without. Twins share a mother and are
  not independent.

## Processing (fixed)

- Counts → CPM. Keep genes with CPM ≥ 1 in ≥ 10% of samples; log2(CPM + 1). Gene z-scores across the analysis
  set.
- Scores (mean z, standardised):
  - **`hif2dep`:** `artifacts/addendum-l/L1_HIF2dep_genes.txt`;
  - **`hyp_only`:** Hallmark hypoxia minus glycolysis;
  - **FLT1:** gene z.

## Tests (OLS, HC3; `y ~ infarct + sex + chorio + fvni + smoke + mage + bmi`)

**X1 (confirmatory).** `hif2dep` is higher with infarction.

- (a) One-sided p < 0.05 (Holm with X2); **and**
- (b) the coefficient is above the 95th percentile of 1,000 random gene sets of the same size, decile-matched on
  mean expression, without hypoxia, HIF2dep or FLT1 genes. Seed 0.

**X2 (confirmatory).** FLT1 is higher with infarction.

- (a) One-sided p < 0.05 (Holm with X1); **and**
- (b) the coefficient is above the 95th percentile of 1,000 random single genes from FLT1's expression decile.
  Seed 0.

**X3 (descriptive; the dissociation).** The same model for `hyp_only`, with its random-set percentile. W predicts
that it ranks lower than X1.

**Sensitivity analyses (descriptive):**

- (a) all births, with cluster-robust SE by mother (title without its last digit, which groups twin and triplet
  sets);
- (b) without the inflammation covariates (chorio, fvni).

## Pre-written readings

- **X1 and X2 supported.** "In an independent cohort of extremely preterm placentas, infarction, a malperfusion
  lesion, goes with higher HIF-2 program and FLT1 expression beyond random genes. This replicates the post hoc
  lead from Addendum W."
  - It is association at delivery.
  - Infarcted tissue may be sampled differently.
  - Preeclampsia is unrecorded and may mediate the effect.
- **Only one supported.** Reported as partial.
- **Neither.** "The Addendum W lead does not replicate." Link 1 then stays unsupported.
