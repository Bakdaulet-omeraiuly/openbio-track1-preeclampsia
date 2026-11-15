# Addendum N: clinical anchors of the placental hypoxia score in GSE75010 (exploratory; written before running)

Written on 2026-11-15 before `n_explore.py` was run. It is committed together with the script.

**Status: exploratory.** GSE75010 expression was opened once in Addendum K, which computed the hypoxia scores
and tested A1–A4 and B1. None of the clinical associations below were computed or viewed. Because the data
have been opened, nothing here counts as confirmation. Results are reported whatever they are.

## Questions and predictions

Scores, gene sets and preprocessing are exactly as in `artifacts/addendum-k/k_confirm_gse75010.py`. Scores
are mean z of `hyp_only` (HALLMARK_HYPOXIA minus HALLMARK_GLYCOLYSIS) and of `hyp_all`; FLT1 is gene z.

### N1: early- vs late-onset PE

The "two diseases" view holds that placental malperfusion drives early-onset PE more than late-onset PE.

- **Early contrast:** PE delivered <34 weeks vs non-PE, non-IUGR deliveries <34 weeks. This is the K-A1
  contrast, recomputed for comparison and not new.
- **Late contrast:** PE delivered ≥34 weeks vs non-PE, non-IUGR deliveries ≥34 weeks.
- **Model:** OLS of `hyp_only ~ PE + GA + sex + attempted labour`, within each GA stratum.
- **Interaction:** an interaction model on both strata together, `hyp_only ~ PE * early + GA + sex + labour`.
- **Prediction:** the late contrast is smaller than the early contrast, i.e. the interaction coefficient is
  > 0. This is reported, not gated.

### N2: maternal uterine-artery Doppler (primary exploratory question)

Raised uterine-artery pulsatility index (PI) measures impaired maternal perfusion of the placenta.

- **Statistic:** partial Spearman of `hyp_only` vs `mean uterine pi`, adjusting for GA. It uses all
  placentas with a recorded PI; the n is reported.
- **Prediction:** positive.
- **Secondary:** the same for FLT1 (prediction: positive) and for `mean umbilical pi` (descriptive).

### N3: severity within PE

Within the 80 PE placentas, partial Spearman of `hyp_only` vs each of the following, adjusting for GA:

| Variable | Predicted sign |
|---|---|
| maximum systolic BP | positive |
| newborn weight z-score | negative |
| placental weight z-score | negative |
| HELLP (yes = 1) | positive |

**Multiplicity.** Holm is applied across one family: N2 primary plus the four N3 tests. All p-values are
one-sided in the predicted direction. N1 and the N2 secondaries are descriptive.

## What would count against the hypothesis

- No association with uterine PI (N2).
- A late-onset effect as large as the early-onset one (N1).

Either would weaken the reading of the score as a marker of placental malperfusion.
