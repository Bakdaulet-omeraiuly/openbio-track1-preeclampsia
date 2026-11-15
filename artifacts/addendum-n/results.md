# Addendum N results (exploratory; GSE75010 previously opened in Addendum K)

The plan and script were committed before running (`494bde1`). Full output is in `N_results.json`.
Partial Spearman values adjust for GA. Holm covers N2 plus the four N3 tests; p-values are one-sided in the
predicted direction.

| Question | Result | Reading |
|---|---|---|
| N1 early PE vs <34-week controls (= K-A1) | +0.23 [0.06, 0.41] (49 vs 27) | as in K |
| N1 late PE vs ≥34-week controls | +0.16 [0.02, 0.31] (31 vs 36) | late-onset PE is also raised |
| N1 interaction (early minus late) | +0.17 [−0.03, 0.37], p = 0.09 | direction as predicted, not clear |
| **N2 hypoxia vs uterine-artery PI** | rho −0.01, n = 52, Holm p = 1.0 | **no association** |
| N2s FLT1 vs uterine PI | rho 0.19, n = 52, p = 0.09 | weak, not clear |
| N2s hypoxia vs umbilical PI | rho 0.17, n = 92, p = 0.05 | weak, descriptive |
| N3 within PE: max systolic BP | rho 0.12, Holm p = 0.57 | no clear association |
| N3 within PE: newborn weight z | rho 0.01, Holm p = 1.0 | none |
| N3 within PE: placental weight z | rho 0.05, Holm p = 1.0 | none |
| N3 within PE: HELLP | rho 0.16, Holm p = 0.38 | no clear association |

## Reading

- **The placental hypoxia score did not track maternal or fetal clinical measures in this cohort.** That
  includes uterine-artery Doppler, the most direct perfusion measure available.
- **Uterine PI was recorded mostly in PE pregnancies** (37 PE vs 15 non-PE), so its range among controls is
  narrow.
- **N1 is consistent with a stronger early-onset signal but does not establish it.** Late-onset PE was also
  raised over GA-matched controls.
- **What this means for the hypothesis.** The score behaves as a marker of placental cell state, tied to FLT1,
  rather than as a readout of clinical severity or measured perfusion. By the plan's own rule, the N2 null
  weakens the reading of the score as a malperfusion marker.
- **Limits.** These are exploratory analyses in previously opened data. Each N3 test is within 80 PE
  placentas, which gives ~80% power only for rho ≈ 0.28 or larger.
