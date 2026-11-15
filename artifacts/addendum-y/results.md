# Addendum Y results (2026-11-15): placental FLT1 does not clearly track maternal endothelial disease (link 6 not supported)

The plan and script were committed (`37916b0`) before running. The script was run **once**. Output:
`Y_results.json`.

## Y1: GSE75010, within 80 PE placentas, partial Spearman with GA

| Outcome | ρ | n | p (one-sided) | Holm p | Verdict |
|---|---|---|---|---|---|
| maximum systolic BP | +0.11 | 80 | 0.16 | 0.32 | not supported |
| maximum diastolic BP | +0.16 | 80 | 0.074 | 0.22 | not supported |
| proteinuria (mode, 0–4) | +0.19 | 78 | 0.044 | 0.17 | not supported |
| HELLP | +0.05 | 80 | 0.33 | 0.33 | not supported |

All four point in the predicted direction; none survives multiplicity.

## Y2: GSE303463, PE + siPE (76) vs gHTN + cHTN (24), `FLT1 ~ PE + GA + sex + batch`

| Coefficient | p (one-sided) | Random same-decile genes: median / 95th pct | Position | Verdict |
|---|---|---|---|---|
| +0.32 SD (SE 0.25) | 0.093 | −0.06 / +0.33 | 94.8th percentile | **Not supported** (fails both conditions, narrowly) |

## Y3 (descriptive): gHTN + cHTN (24) vs Normal (68)

FLT1 is **+0.68 SD** (p = 0.0007) in hypertension without PE.

## Reading (by the rules written before the run)

- **Link 6 is not supported at the placental end in these data.**
  - Within PE, placental FLT1 is not clearly higher with worse maternal blood pressure, proteinuria or HELLP.
  - In an independent cohort it is not clearly higher in PE than in hypertension without end-organ disease.
- **Placental FLT1 is already raised in hypertension alone (Y3).** In these data FLT1 marks hypertensive
  pregnancy in general, and does not specifically mark the endothelial disease that defines preeclampsia.
- The causal step (sFlt-1 → maternal endothelium) remains supported **only by the literature**. Placental mRNA
  at delivery is a weak proxy for circulating sFlt-1 protein.

## Limits

- Placental FLT1 mRNA, not serum sFlt-1 protein.
- Severity variables are maxima recorded clinically.
- Y2 has only 24 hypertension-alone placentas.
- GSE75010 has been used several times before (disclosed in the plan).
