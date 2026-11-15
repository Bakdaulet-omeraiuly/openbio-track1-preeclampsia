# Addendum M: hypoxia program in early-onset PE against gestational-age-matched preterm controls, three untouched cohorts

Written on 2026-11-15 before any expression value from these three cohorts was downloaded or read. The plan
and `m_confirm.py` are committed together. Each cohort is then opened **once**, and the result is reported
whatever it is.

## What was viewed before this plan

- The GEO sample metadata (characteristics, titles, labels, value definition) of the cohorts below and of the
  candidates we did not use.
- The GPL annotation tables (probe to gene symbol), for the coverage gate.
- No expression values.

## Why these cohorts (selection rule)

The rule had three parts:

- placental bulk expression;
- early-onset or preterm PE, plus **controls delivered preterm** or matched by gestational age (GA) by design;
- not used anywhere in this repository.

Screened candidates:

| GEO | Decision | Reason |
|---|---|---|
| GSE44711 (Blair et al. 2013, Robinson lab, Illumina HT-12 v4) | **used** | 8 EOPE vs 8 GA-matched controls, chorionic villi, GA per sample |
| GSE74341 (Liang et al. 2016, Agilent G3 v2, two-channel) | **used** | 7 EOPE (<34 wk) vs 5 controls delivered <34 wk |
| GSE66273 (Than/Romero/Tarca, Agilent 4x44K) | **used** (with a tissue caveat) | 12 preterm PE (6 with HELLP) vs 5 preterm controls, GA and sex per sample. The tissue is **chorionic plate**, not villi |
| GSE148241 (RNA editome) | not used | it is named in our hypothesis card as a next cohort, but its metadata show 9 EOSPE vs 34 "Normal" with no preterm control group, so it cannot test this design |
| GSE25906, GSE10588 | not used | controls are mostly term deliveries, the same problem that made L-C2 hard to read |

**Overlap with cohorts already used.** By GEO contributors and PubMed IDs, none of the three shares a group
with the cohorts already in this repository:

- already used: GSE203507 (Ackerman), GSE114691 (Han), GSE75010 (Cox/Leavey), GSE35574 (Tsai/Piedrahita),
  GSE24129 (Nishizawa), GSE147776 (Medina-Bastidas), GSE216799;
- new: GSE44711 (Robinson), GSE74341 (Liang/Niu), GSE66273 (Than/Romero).

The new cohorts come from three different groups.

## Coverage gate (GPL annotation only, passed)

Rule: FLT1 must be present, and at least 80% of the genes of each set must map.

| Platform | FLT1 | hyp_only | hyp_all | HIF2dep |
|---|---|---|---|---|
| GPL10558 (GSE44711) | yes | 125/135 | 188/200 | 87/100 |
| GPL16699 (GSE74341) | yes | 127/135 | 191/200 | 88/100 |
| GPL4133 (GSE66273) | yes | 125/135 | 188/200 | 85/100 |

## Definitions (fixed)

**Scores.** Identical to `artifacts/addendum-l/l3_confirm.py`:

- probes are averaged per gene symbol;
- each gene is z-scored across **all** arrays of the cohort;
- the set score is the mean z, standardised across arrays;
- `FLT1` is the gene z.

**Sets.**

- `hyp_only` = Hallmark hypoxia minus Hallmark glycolysis (v2024.1);
- `hyp_all` = Hallmark hypoxia;
- `hif2dep` = `artifacts/addendum-l/L1_HIF2dep_genes.txt`.

**Values.**

- GSE74341 values are the authors' LOWESS log2(Cy3/Cy5) ratio, which is **sample over pooled reference**
  (Cy3 = sample, ch2; Cy5 = pooled reference, ch1). The script asserts that every array has `label_ch2 = Cy3`,
  so no dye swaps are present. Group labels come from `characteristics_ch2`.
- For the single-channel cohorts, log2(x + 1) is applied only if the matrix maximum is above 100 (the L3 rule).

**Groups.**

| Cohort | PE | Control | Excluded from M1/M2 | Covariates in M1 |
|---|---|---|---|---|
| GSE44711 | EOPET (8) | Control (8) | none | GA |
| GSE74341 | EOPE (7) | preterm controls <34 wk (5) | LOPE (8) and term controls (5); they are still used in z-scoring | none (no per-sample GA; both groups <34 wk) |
| GSE66273 | preterm PE ± HELLP (12) | preterm controls (5) | none | GA + fetal sex |

## Tests

**M1 (confirmatory).** `hyp_only` is higher in PE than in GA-matched or preterm controls.

- Per cohort: OLS `hyp_only ~ PE + covariates`.
- Pooled: fixed-effect inverse-variance estimate across the 3 cohorts; τ² is reported.
- **Decision p:** one-sided stratified permutation p (labels permuted within each cohort, 10,000 permutations,
  seed 42; the statistic is the pooled fixed-effect z).
- The normal-theory p is reported but does not decide.

**M2 (confirmatory).** Within PE, `hyp_all` correlates positively with FLT1.

- Per cohort: Spearman correlation.
- Pooled: Fisher z with n − 3 weights.
- **Decision p:** one-sided stratified permutation p (FLT1 permuted within PE within each cohort, 10,000
  permutations, seed 42).

**M3 (descriptive only, no verdict).** Within PE, the partial Spearman correlation of HIF2dep with FLT1 given
`hyp_all`. It is pooled as in M2 and reported with its CI.

- The subset is **PE only**, because these cohorts have no FGR-only placentas; in L-C4 the subset was PE + FGR.
- With about 27 PE placentas this test is underpowered, so it is not in the Holm family.

**Multiplicity.** Holm across M1 and M2. A test is **supported** if its Holm-adjusted permutation p is below
0.05 and the pooled estimate is positive.

**Pre-specified sensitivity analyses (they do not change the verdict):**

- (a) leave one cohort out, which includes **M1 without GSE66273** (villous or placental tissue only);
- (b) GSE66273 without the HELLP cases;
- (c) GSE44711 with controls < 34 weeks only;
- (d) M1 without covariates.

## Pre-written readings

- **M1 supported.** "Early-onset PE placentas show a higher hypoxia-only score than GA-matched preterm controls
  in this pre-registered test. Counting earlier cohorts, this holds in N of N cohorts with preterm controls."
  - This is **not** a test of why L-C2 failed. Those cohorts differ in more than control GA, so the term-control
    explanation stays untested.
  - If the result depends on GSE66273 (sensitivity a), we say so.
- **M1 positive but not supported.** "Same direction, not significant in small cohorts (about 27 PE vs 18
  controls); not a confirmation."
- **M1 at or below zero.** "Not replicated in GA-matched cohorts. This weakens the 'PE vs preterm control'
  claim, which then rests on three earlier cohorts."
- **M2.** Supported or not by the same rule. A failure would weaken our most-replicated result and must be
  stated in the article.
- **M3.** An estimate and CI only. It cannot confirm or refute L-C4.
