# Addendum W results (2026-11-15): MVM does not raise the generic hypoxia program in GSE303463 (link 1 not supported)

The plan (`7cd1559`) and amendment 1 (`5f9ab27`, ID parsing only) were committed before any count value was
read. The script was run **once**. Outputs: `W_results.json`; input checksum: `inputs.sha256`.

Data: 168 placentas, 49 MVM and 119 nonMVM. Genes used: hyp_only 122, hyp_all 179, HIF2dep 83.

## W1 (confirmatory): `hyp_only ~ MVM + GA + sex + batch + disease`

| | Estimate |
|---|---|
| MVM coefficient | **+0.12 SD** (HC3 SE 0.18), one-sided **p = 0.25** |
| 1,000 random gene sets, same model | median **+0.20**, 95th percentile +0.29; observed at the **7.9th percentile** |
| **Verdict** | **Not supported** (fails both conditions) |

## Pre-specified sensitivity analyses

| Analysis | MVM coefficient | p (one-sided) |
|---|---|---|
| (a) hypertensive disorders only (45 vs 55) | −0.01 | 0.53 |
| (b) without the disease covariate | +0.18 | 0.13 |
| (c) `hyp_all` instead of `hyp_only` | +0.13 | 0.22 |

## W2 / W3 (descriptive, no verdict)

| Outcome | MVM coefficient (same model) | p (one-sided) |
|---|---|---|
| HIF2dep program | +0.46 SD | 0.012 |
| FLT1 | +0.61 SD | 0.0007 |

**Hypoxia–FLT1 coupling (Spearman):** 0.41 within MVM vs 0.17 within nonMVM.

The random gene sets themselves shift by about +0.2 with MVM (W1 null), so the W2 coefficients should be read
against that background. No random-set null was run for W2.

## Post hoc calibration of W2 (`w2_calibration_posthoc.py`, committed `28e8ad3` before running; a lead, not a confirmation)

| Outcome | MVM coefficient | Random comparison | Position |
|---|---|---|---|
| HIF2dep program | +0.46 SD | 1,000 random sets of 83 genes (decile-matched): median +0.16, 95th pct +0.27 | **above all 1,000** |
| FLT1 | +0.61 SD | 1,000 random genes from FLT1's expression decile: median +0.11, 95th pct +0.48 | **98.6th percentile** |
| FLT1 genome-wide | — | all 17,128 expressed genes (median +0.05) | **rank 121** (99.3rd percentile) |
| (for contrast, W1) generic `hyp_only` | +0.12 SD | random sets: median +0.20 | 7.9th percentile |

**Post hoc reading.** With MVM, the **HIF-2 program and FLT1** shift beyond random expectation, while the
**generic Hallmark hypoxia score** does not. This dissociation in tissue mirrors the genetic one in Addendum V,
where HIF-2 loss removes FLT1 and the program and HIF-1 loss removes canonical hypoxia genes.

- It was chosen after seeing W2, so it needs a pre-registered test in another cohort with histopathology before
  it can count.
- That FLT1 or sFlt-1 is higher with MVM is consistent with prior clinical-pathology reports. The program-level
  contrast with generic hypoxia is our observation.

## Reading (by the rule written before the run)

- **W1 is not supported.** Placentas with histological MVM do not carry a stronger *generic* hypoxia program than
  those without it, once disease group, GA, sex and batch are accounted for. They do not differ from random gene
  sets either. Within hypertensive disorders the difference is zero.
  - **Link 1 therefore stays literature-only, and it is weakened in our data.** This agrees with the null Doppler
    result in Addendum N.
- **Descriptive pattern, not a finding.** FLT1 and the HIF2dep program are higher in MVM placentas, while the
    generic Hallmark hypoxia score is not.
  - This fits the HIF-2 / FLT1 axis of Addendum V more than a canonical HIF-1 response.
  - It was not pre-specified as confirmatory and has no random-set calibration, so it is a lead only.

## Limits

- MVM is a binary pathology call made at delivery. Expression is measured at delivery too, so the timing of
  malperfusion relative to sampling is unknown.
- MVM is strongly tied to disease (4/68 Normal vs 45/100 hypertensive), so the disease-adjusted contrast rests
  mostly on the hypertensive placentas.
- Part of the series (GSE186257) was used in Addendum P for a different question.
