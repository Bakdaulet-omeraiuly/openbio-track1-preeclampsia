# Addendum X results (2026-11-15): infarction goes with higher FLT1, but not with the HIF2dep program (ELGAN, GSE154829)

The plan (`dd8b1cf`) and amendment 1 (`20ed690`) were committed before any analysis.

- Amendment 1 restricts the analysis to mRNA libraries. The first run stopped while loading, because the series
  also lists miRNA libraries. That also corrects the sample counts quoted in the plan.
- The script was then run **once**. Outputs: `X_results.json`; input checksum: `inputs.sha256`.

Primary set: singleton births with complete covariates: **40 with infarction vs 163 without**. The model is
`y ~ infarct + sex + chorioamnionitis + fetal vessel neutrophilic infiltration + smoking + maternal age + BMI`,
with HC3 SE.

| Test | Infarct coefficient | One-sided p (Holm) | Random comparison: median / 95th pct; position | Verdict |
|---|---|---|---|---|
| **X2** FLT1 | **+0.42 SD** (SE 0.16) | **0.005 (0.011)** | random same-decile genes: −0.10 / +0.13; **above all 1,000** | **Supported** |
| **X1** HIF2dep program | +0.03 SD (SE 0.18) | 0.43 (0.43) | random sets: −0.06 / −0.01; 99.8th percentile | **Not supported** (fails the p condition) |
| X3 generic `hyp_only` (descriptive) | −0.03 SD | 0.56 | random sets: −0.07 / −0.03; 96.7th percentile | — |

In this cohort the random-set null for gene-set scores is centred below zero, so a near-zero program coefficient
ranks high. That is why X1 needs both conditions.

## Sensitivity (descriptive)

| Analysis | FLT1 | HIF2dep | `hyp_only` |
|---|---|---|---|
| (a) all births, SE clustered by mother (55 vs 285) | +0.31, p = 0.023 | −0.07 | −0.15 |
| (b) without inflammation covariates | +0.46, p = 0.004 | +0.06 | −0.01 |

## Reading (by the rules written before the run)

- **Partial replication of the Addendum W lead.**
  - **FLT1 replicates.** In an independent cohort of extremely preterm placentas, infarction goes with higher
    FLT1 beyond covariates and beyond expression-matched random genes. This is consistent with prior clinical
    links between malperfusion lesions and sFlt-1.
  - **The HIF2dep program does not** replicate: about zero.
  - **The generic hypoxia score does not move**, again, as in W1.
- **For link 1, generic hypoxia is not raised by malperfusion lesions in two cohorts** (W, X), nor tracked by
  Doppler (N). FLT1 is raised by them in two cohorts (W post hoc, X pre-registered).
  - So in tissue, malperfusion links to FLT1 more robustly than to any hypoxia program we measure.
  - This could mean that malperfusion raises FLT1 by routes not captured by these scores, or that tissue
    hypoxia at delivery is not what bulk scores measure.

## Limits

- Placentas from extremely preterm births (23–27 weeks).
- Infarction is one MVM lesion and is binary.
- Preeclampsia is unrecorded and may mediate the effect.
- Expression is measured at delivery.
- The infarcted region may differ in sampling.
