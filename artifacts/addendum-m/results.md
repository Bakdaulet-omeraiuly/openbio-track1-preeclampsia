# Addendum M results (2026-11-15)

The plan and script were committed and pushed (`0b12ad5`) before any series matrix was downloaded. The three
cohorts were then downloaded and the script was run **once**, with no changes to it. Output:
`M_confirmation.json`.

## Confirmatory tests (Holm across M1 and M2; decision p = stratified permutation, 10,000 permutations, seed 42)

| Test | Pooled estimate [95% CI] | Permutation p (one-sided) | Holm p | Verdict |
|---|---|---:|---:|---|
| **M1** hyp_only, PE vs GA-matched or preterm controls | **+0.67 SD [0.06, 1.28]**, τ² = 0 | 0.033 | **0.033** | **Supported** |
| **M2** hyp_all–FLT1 Spearman correlation within PE | **0.47 [0.05, 0.75]** | 0.016 | **0.031** | **Supported** |
| M3 (descriptive only) HIF2dep–FLT1 partial ρ given hyp_all, within PE | 0.45 [0.02, 0.74] | — | — | no verdict (pre-specified) |

Per cohort:

| Cohort (tissue) | n PE / control | M1 estimate (SE) | M2 ρ | M3 partial ρ |
|---|---|---|---|---|
| GSE44711 (chorionic villi) | 8 / 8 | +1.05 (0.44) | 0.61 | 0.33 |
| GSE74341 (placental tissue) | 7 / 5 | +0.61 (0.69) | 0.68 | 0.86 |
| GSE66273 (chorionic plate) | 12 / 5 | +0.07 (0.57) | 0.26 | 0.20 |

## Pre-specified sensitivity analyses (normal-theory p; they do not change the verdict)

| Analysis | M1 [95% CI], p | M2 [95% CI] |
|---|---|---|
| without GSE66273 (villous or placental tissue only) | +0.92 [0.19, 1.65], p = 0.007 | 0.64 [0.11, 0.89] |
| without GSE74341 | +0.68 [0.00, 1.37], p = 0.025 | 0.40 [−0.10, 0.74] |
| **without GSE44711** | **+0.29 [−0.58, 1.16], p = 0.26** | 0.41 [−0.11, 0.75] |
| GSE66273 without the HELLP cases | +0.29 [−0.16, 0.74], p = 0.11 | 0.54 [0.04, 0.83] |
| GSE44711 controls < 34 wk only (6 controls) | +0.58 [−0.09, 1.24], p = 0.045 | — |
| M1 without covariates | +0.67 [0.07, 1.27], p = 0.014 | — |

## Reading (by the rules written before the run)

- **M1 is supported.** Early-onset or preterm PE placentas have a higher hypoxia-only score than GA-matched or
  preterm controls in this pre-registered test.
  - Counting the earlier cohorts with preterm controls (GSE203507, GSE114691, GSE75010), the direction is
    positive in all six such cohorts. It is significant pooled, but not in every single small cohort.
  - The result **leans on GSE44711.** Without it, the pooled estimate is +0.29 with a CI that includes zero.
  - The chorionic-plate cohort (GSE66273) shows almost no difference.
  - With 27 PE vs 18 controls, the CI is wide (0.06 to 1.28 SD).
  - This does **not** test why L-C2 failed against term controls. That explanation stays untested.
- **M2 is supported.** The hypoxia score tracks FLT1 within PE in three more independent cohorts. Together with
  Addenda C, K and L, the coupling has now been tested in **nine bulk cohorts** and is positive in every one.
- **M3 is descriptive.** Its estimate (0.45) is close to L-C4 (0.39), but PE-only and small. It cannot confirm L-C4.

## Post hoc, descriptive (proposed after the data were opened): placental FLT1:PGF log2 ratio

Script: `m4_flt1_pgf_posthoc.py`. Output: `M4_flt1_pgf_posthoc.json`. Covariates are as in M1, except fetal
sex in GSE66273.

| Cohort | PE − control, log2 (SE) |
|---|---|
| GSE44711 | +0.59 (0.90) |
| GSE74341 | +2.48 (0.78) |
| GSE66273 | +3.44 (0.68) |
| pooled, villous cohorts only | +1.67 [0.51, 2.83] |

**Reading.** This is a sanity check that the cohorts carry the known PE placental signature: FLT1 up and PGF down.
It is *not* a confirmatory test.

- The ratio is strongly raised in the chorionic-plate cohort, where the hypoxia score is not.
- So the null M1 there is not simply "the tissue shows nothing". In that cohort the clinical-marker signal is
  present without the hypoxia program, which fits our earlier finding that only a subset of PE placentas carry
  it. This is an untested reading.

## Processing disclosures

- `scores()` in `m_confirm.py` drops genes with any missing value before z-scoring; L3 did not. This was fixed
  in the committed script before the run.
- In GSE44711 (Illumina, maximum > 100), log2(x + 1) was applied. GSE74341 and GSE66273 were used as provided
  (already on a log scale).
- GSE74341 values are log2(sample/pooled reference). Every array was asserted to have Cy3 = sample, so there are
  no dye swaps.
- File checksums for the downloaded matrices and GPL tables were recorded locally at download. The URLs are in
  the script header.
