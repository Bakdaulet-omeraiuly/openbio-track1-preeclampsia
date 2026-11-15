# Addendum Y: does placental FLT1 track maternal disease? (link 6, GSE75010 + GSE303463)

Written on 2026-11-15, before FLT1 was related to any maternal-severity variable in either cohort. The plan and
`y_maternal.py` are committed together, then run **once**.

## Why, and what is tested

Link 6 is "sFlt-1 → maternal endothelial disease". The causal step, sFlt-1 damaging maternal endothelium, rests
on animal and clinical literature (Maynard et al. 2003; Levine et al. 2004), and expression data cannot test it.

What can be tested is the placental end of that step: **is placental FLT1 expression higher when the mother has
more endothelial disease?** This is association only.

## Prior use (disclosure)

- **GSE75010** was used in Addenda K, L (screen) and N. N related the *hypoxia score* (not FLT1) to systolic BP,
  HELLP, newborn and placental weight. FLT1 was related only to uterine-artery PI (N2s). **FLT1 vs maternal BP,
  proteinuria or HELLP has not been computed.**
- **GSE303463** was used in Addendum W for MVM (FLT1 vs MVM, descriptive and post hoc). **FLT1 vs disease
  subtype has not been computed.**

## Tests

**Y1 (GSE75010, within the 80 PE placentas; Holm across the 4 outcomes).** Partial Spearman of FLT1 (gene z, as
in Addendum K) with each outcome, adjusting for GA:

| Outcome | Predicted sign |
|---|---|
| maximum systolic BP | positive |
| maximum diastolic BP | positive |
| mode proteinuria (ordinal 0–4; NA excluded) | positive |
| HELLP (yes = 1) | positive |

**Y2 (GSE303463, independent; confirmatory).** FLT1 is higher in preeclampsia with end-organ (endothelial)
disease than in hypertension without it.

- Groups: **PE + siPE** (n = 76) vs **gHTN + cHTN** (n = 24).
- Model: `FLT1_z ~ PE + GA + sex + batch` (HC3), one-sided.
- **Supported** if p < 0.05 **and** the PE coefficient is above the 95th percentile of 1,000 random genes from
  FLT1's expression decile, under the same model. Seed 0.

**Y3 (descriptive).** In GSE303463, the same model for Normal vs gHTN + cHTN. This shows whether FLT1 is already
up with hypertension alone.

## Pre-written readings

- **Y2 supported.** "Placental FLT1 is higher when hypertension comes with the endothelial disease that defines
  preeclampsia than with hypertension alone, beyond random genes."
  - This is consistent with link 6 at the placental end; the causal step itself is from the literature.
- **Y1.** Each outcome is reported with its Holm p. Any supported outcome adds within-PE severity support.
- **Neither.** "Placental FLT1 does not track maternal disease in these data." Link 6 then rests on the
  literature only.
