# Addendum W: does maternal vascular malperfusion (MVM) raise the placental hypoxia program? (link 1, GSE303463)

Written on 2026-11-15, before any expression value of GSE303463 was downloaded. The plan, `w_mvm.py` and
`samples.tsv` (from GEO metadata only) are committed together, then run **once**.

## Why

Link 1 is "malperfusion → placental hypoxia". So far it rested on the literature only. Addendum N found no link
between the hypoxia score and uterine-artery Doppler in GSE75010. MVM is the histopathological diagnosis of
inadequate maternal perfusion, so it is the direct tissue-level measure of link 1.

## Data

- **GSE303463** (Horii, Parast et al.; PMID 41467353): 168 term or preterm placentas, bulk RNA-seq (raw counts).
- Per-sample metadata: placental pathology MVM (49) vs nonMVM (119); disease (Normal 68, siPE 54, PE 22,
  cHTN 13, gHTN 11); GA (weeks); sex; batch (5).
- **Overlap.** This series re-uses GSE186257, which Addendum P used for hypoxia–FLT1 coupling. **The MVM labels
  have never been used in this repository.**
- **Prior work.** The depositors analysed placental histopathology and molecular data together. MVM-associated
  pathways may already have been reported, so this is a pre-registered, confounder-adjusted test of our specific
  scores, not a claim of first observation.

**MVM is strongly confounded with disease.** It is present in 4/68 Normal vs 45/100 hypertensive placentas. The
primary model therefore adjusts for disease group.

## Processing (fixed)

- Raw counts → CPM. Keep genes with CPM ≥ 1 in ≥ 10% of samples; log2(CPM + 1).
- Scores, as in our bulk addenda (mean z of set genes, then standardised):
  - **`hyp_only`:** Hallmark hypoxia minus glycolysis;
  - **`hyp_all`:** Hallmark hypoxia;
  - **`hif2dep`:** the HIF2dep program;
  - **FLT1:** z of log2 CPM.

## Tests

**W1 (confirmatory).** `hyp_only ~ MVM + GA + sex + batch + disease` (OLS). Prediction: the MVM coefficient is
> 0.

**Supported** if both hold:

- (a) the one-sided p is < 0.05 (HC3 robust SE); and
- (b) the MVM coefficient is above the 95th percentile of 1,000 random gene sets.
  - The sets have the same size as `hyp_only` and are drawn within mean-expression deciles, without hypoxia or
    HIF2dep genes. Seed 0.
  - Each random set is scored and modelled identically.

**Pre-specified sensitivity analyses (descriptive):**

- (a) within hypertensive disorders only (cHTN + gHTN + PE + siPE; 45 MVM vs 55 nonMVM), same covariates minus
  Normal;
- (b) without the disease covariate;
- (c) `hyp_all` instead of `hyp_only`.

**W2 (descriptive).** The same model for `hif2dep` and for FLT1. Prediction: both MVM coefficients > 0.

**W3 (descriptive).** Does the hypoxia–FLT1 coupling differ by MVM? Spearman of `hyp_all` vs FLT1 within MVM and
within nonMVM.

## Pre-written readings

- **W1 supported.** "Placentas with histological maternal vascular malperfusion carry a stronger hypoxia program
  than placentas without it, beyond disease group, gestational age, sex and batch, and beyond random gene sets.
  This is tissue-level support for link 1 (malperfusion → hypoxia)."
  - It is still association: MVM and expression come from the same placenta at delivery.
- **W1 not supported.** "No MVM effect on the hypoxia program beyond disease and GA in this cohort." Link 1 then
  remains literature-only and is weakened. This would be consistent with the null Doppler result (N).

## Amendment 1 (before any count value was read; only the header row and the ID column were viewed)

Gene IDs are written as `ENSG…$SYMBOL`. The symbol is taken after `$`, and counts of duplicated symbols are
summed. Sample columns are the GEO titles (e.g. P1007).
