# Analysis Plan

> Written and committed **before** any differential-expression or pathway test was run (see git history). Anything added after results were seen is marked *post hoc*.

## Research Question

Do hypoxia and glycolysis gene-activity signals in early-onset preeclampsia (EOPE) placentas persist when cases are compared with gestational-age-appropriate (preterm) controls, and which placental cell types plausibly carry those signals?

Motivation: Ackerman et al. (2023) reported a glycolytic expression pattern in a subset of PE placentas in this same cohort (GSE203507). Our contribution is to test how robust that signal is to the choice of controls and to delivery timing, to look at it donor by donor, and to place the genes in cell types.

## Datasets

| Dataset ID | Role | Why This Dataset |
|---|---|---|
| GSE203507 villous bulk RNA-seq (`track01_counts_vt.csv`, `track01_metadata_vt.csv`) | Primary analysis | Required core data; 30 donors, one villous sample each |
| MSigDB Hallmark v2024.1 (HYPOXIA, GLYCOLYSIS) | Predefined gene sets | Fixed, published gene lists chosen before seeing results |
| `track01_marker_matrix.csv` | Cell-mix vs cell-behaviour check | Offline marker genes per broad cell class |
| Zeisel/Admati gene browser | Cell-type context for ≥3 genes | Required; independent single-cell study |
| Admati/Skarbianskis pseudobulk (optional) | Within-cell-type check of selected genes | Only if time allows; exploratory (3 early controls) |

## Methods

1. **Data loading and checks.** Use the notebook's checks. Do not edit data. Each donor contributes one villous sample, so donor = sample = unit of replication.
2. **Primary comparison.** `pe_status`: PE (EOPE + EOPE+FGR, n = 15) vs. preterm controls (`control_type == preterm_birth`, iPTB, n = 5). FGR-only and term controls are excluded from the primary model. Reason: term controls (≈39 wk) are delivered ~6–13 weeks later than PE cases (26–34 wk), so PE and gestational age are confounded; iPTB controls overlap the PE gestational-age range (29.9–33.0 wk).
3. **Statistical model and thresholds.** PyDESeq2, design `~ fetal_sex + pe_status`, Wald test, Benjamini–Hochberg. A gene is a hit if `padj < 0.05` and `|log2FC| ≥ 1`, with `baseMean ≥ 10`. These thresholds are fixed now. We will report how the hit count changes at padj < 0.1 and |log2FC| ≥ 0.58 for transparency, but conclusions use the fixed thresholds.
4. **Pathway test (predefined).** For HALLMARK_HYPOXIA and HALLMARK_GLYCOLYSIS only:
   - Per-donor score = mean z-score of variance-stabilised expression across set genes present in the data.
   - Compare PE vs iPTB scores (Mann–Whitney, two-sided), and fit `score ~ PE + gestational age` (OLS) to see whether the PE effect survives timing.
   - Plot every donor (no averages only) to show heterogeneity within PE.
   - Supporting competitive test: rank of set genes in the DESeq2 Wald statistic (Wilcoxon rank-sum vs other genes; noted as anti-conservative because genes are correlated).
   - Known PE anchor genes (FLT1, ENG, LEP, PGF) reported as a positive-control sanity check, not as discoveries.
5. **Gene browser lookups.** At least three genes that (a) are hits or are leading members of a set, (b) have baseMean ≥ 100, (c) survive the robustness checks below. Pick some up and some down where possible.
6. **Robustness / sensitivity analyses** (all predefined):
   - S1: primary model + `ga_delivery_weeks` covariate.
   - S2: EOPE without FGR (n = 8) vs iPTB, to separate PE from growth restriction.
   - S3: PE vs term controls (TB, n = 5), expected to inflate hits because of timing; used to show the size of the confound.
   - S4: leave-one-control-out (5 refits of the primary model); report hit overlap.
   - Report how many primary hits keep the same sign and padj < 0.05 in each.
7. **Interpretation.** Use the marker table: if markers of one cell class move together, a cell-mix explanation is favoured (Campbell et al. 2023); if markers are flat but class-specific genes change, a within-cell behaviour change is favoured. Keep observations separate from interpretations.

## Confounders and Controls

| Confounder | How we handle it |
|---|---|
| Gestational age at delivery | iPTB controls; S1 GA covariate; S3 shows the term-control contrast |
| Fetal growth restriction | 7/15 PE cases have FGR; S2 excludes them |
| Fetal sex | Covariate in all models |
| Cell-type composition | Marker-table check; browser; optional pseudobulk |
| Small n (5 controls) | Treat as exploratory; S4 leave-one-out; state limits |
| QC flags (S12 PE without BP ≥140/90) | Keep the sample as labelled by the study; report flag |

## Expected Figures

1. PCA coloured by group and by gestational age (shows the timing confound)
2. Volcano plot, PE vs iPTB
3. Per-donor hypoxia and glycolysis scores vs gestational age, coloured by group (headline candidate)
4. Robustness: hit overlap across S1–S4
5. Marker-class summary (mean log2FC of markers per cell class)

---

## Addendum A (2026-11-14): written after seeing the primary results, before running these tests

These analyses were proposed **after** the results of notebook 01 were seen (in particular the
post hoc observation that FGR-only placentas also have high glycolysis scores). They are follow-up
checks, not part of the original pre-registration, and are reported as such.

### A1. Independent single-cell check (Admati/Skarbianskis pseudobulk, different donors)
- Cell classes chosen in advance: `villous_cytotrophoblast`, `syncytiotrophoblast` (the trophoblast
  compartment where FLT1/PGF dysregulation is reported). No other class will be examined for A1.
- Comparison: `PE_early` (n = 10) vs `control_early` (n = 3), one column per donor. Exploratory.
- (a) Direction concordance: for the 207 robust bulk hits expressed in the class (≥10 reads in ≥3 donors),
  the fraction whose PE-minus-control difference in mean log2 CPM has the same sign as in bulk.
  Null: exact permutation of the PE/control labels (all 286 splits of 13 donors into 10 + 3).
- (b) Hypoxia and glycolysis scores per donor within each class (mean z of log2 CPM, z across the 13
  donors). Statistic: difference in mean score PE − control; exact permutation p (286 splits).
  Leave-one-donor-out to check dependence on a single donor.
- Secondary (exploratory): all 26 donors, `score ~ PE + gestational age`.

### A2. Is glycolysis a shared placental stress signal while hypoxia is more PE-linked?
- Bulk data, 25 preterm samples (TB excluded). OLS `score ~ C(group, ref = iPTB) + gestational age + fetal sex`.
- Report effect sizes with 95% CIs for EOPE, EOPE+FGR and FGR-only versus iPTB.
- **Direct** contrasts (not "significant vs not significant"): EOPE − FGR-only, and
  mean(EOPE, EOPE+FGR) − FGR-only, each with 95% CI. HC3 robust standard errors as a check.

### A3. Are hypoxia and glycolysis separate signals?
- Split the two Hallmark sets into: shared genes (65), hypoxia-only, glycolysis-only.
- Repeat A2 for each of the three scores, and report the per-sample correlation between
  hypoxia-only and glycolysis-only scores.

Multiple testing: A1–A3 involve about 20 tests; all are reported, none is used alone as a
confirmatory claim.

---

## Addendum B (2026-11-14): written after Addendum A results, before running these tests

Main claim being stress-tested: *a hypoxia-response program is higher in early-onset PE villous
tissue than in preterm controls and FGR-only placentas.* Score definitions are as in notebooks 01/02
(mean z of variance-stabilised expression; Hallmark v2024.1).

### B1. Null calibration with random gene sets
Is the hypoxia effect larger than for an arbitrary gene set of the same size and expression level?
- 1,000 random gene sets per score, matched to the real set by expression decile (bins of mean VST).
- Same OLS model as A2 on the 25 preterm samples; statistics: PE (EOPE) vs iPTB coefficient, and the
  direct PE (both) − FGR-only contrast.
- Empirical one-sided p = (1 + #random ≥ observed) / 1001. Scores: hypoxia (all), shared, hypoxia-only.
- Seed 42.

### B2. Link to a clinical biomarker gene (FLT1, encoding sFlt-1)
FLT1 is **not** in HALLMARK_HYPOXIA, so the test is not circular.
- Spearman correlation of FLT1 (VST) with the hypoxia score: within PE only (n = 15, primary) and
  within the 20 PE + iPTB samples.
- Exploratory replication in pseudobulk syncytiotrophoblast, 13 early donors (log2 CPM).

### B3. Heterogeneity within PE (descriptive)
Within the 15 PE placentas, Spearman correlation of the hypoxia score with: max systolic BP, max
diastolic BP, birthweight (lower bound of the 100 g bin), gestational age, and FGR status
(Mann–Whitney). All five reported; no multiplicity correction claimed as confirmatory.

### B4. Where does the signal sit: placental villi vs maternal decidua?
- Decidua basalis bulk (same study, 27 samples). Score computed within decidua (z across decidua samples).
- Decidua PE (EOPE + EOPE+FGR) vs iPTB, same OLS covariates; compare with the villous estimate.
- Direct test: mixed model `score ~ PE * tissue + gestational age + fetal sex + (1 | donor)` on
  PE + iPTB samples from both tissues; the PE×tissue interaction is the statistic of interest.

### B5. Exploratory: adjust for estimated cell composition
- Signature: for each of the 14 pseudobulk cell classes, the 50 genes most specific to that class
  (log2 CPM class mean minus max of other classes, across all 26 donors).
- Estimate bulk proportions by non-negative least squares on CPM; normalise to sum 1.
- Refit the A2 hypoxia model adding estimated syncytiotrophoblast and villous cytotrophoblast
  proportions. Report whether the PE coefficient shrinks. Labelled exploratory (deconvolution across
  platforms is approximate).

---

## Addendum C (2026-11-14): external replication and final robustness, written before running

### C1. Replication in a different patient cohort: GSE114691 (Awamleh et al. 2019)
79 villous placentas (decidua removed): control = spontaneous preterm labour < 34 wk (n = 21),
EO-PE (n = 20), EO-IUGR (n = 18), EO-PE + IUGR (n = 20); the paper reports no difference in
gestational age at delivery between groups. GEO provides transcript-level counts and the group
label only (no per-sample GA, sex or delivery mode), so models are unadjusted. This cohort was also
used by Ackerman et al. 2023 for their module discovery, so it is "different patients", not a fully
untouched validation set. Our specific claims (hypoxia-only vs glycolysis-only split, PE vs FGR
contrast, FLT1 link) were derived only from GSE203507.

Processing (fixed now): sum transcript counts to gene symbols using Ensembl GRCh37 BioMart
(`ensembl_transcript_id` → `external_gene_name`); same low-count filter (≥10 reads in ≥5 samples),
PyDESeq2 VST, z across all 79 samples, the same five gene-set scores as notebook 02.
Model: `score ~ C(group, ref = control)`; contrasts as in A2 (EOPE vs control, EOPE+IUGR vs
control, IUGR vs control, EOPE vs IUGR, mean(PE groups) vs IUGR).

**Pre-specified predictions** (the claim fails if P1 or P2 fails):
- P1: hypoxia-only score: EOPE − control > 0 (one-sided α = 0.05).
- P2: hypoxia-only score: mean(PE groups) − IUGR > 0 (one-sided α = 0.05).
- P3: shared-gene score behaves like P1/P2.
- P4: glycolysis-only score: EOPE − control smaller than the hypoxia-only effect (report both CIs).
- P5: within the 40 PE placentas, Spearman(hypoxia score, FLT1) > 0.
- P6: random-gene-set null (1,000 expression-matched sets, seed 42) for P1 and P2.
Labour in controls can raise hypoxia genes, which would bias against P1.

### C2. FLT1 link within PE only, in the independent single-cell data
Spearman(hypoxia score, FLT1) in syncytiotrophoblast pseudobulk among the 10 early-PE donors only
(the A/B analysis pooled PE and controls). Secondary: all 17 PE donors (early + late).

### C3. Multiple-testing summary
Benjamini–Hochberg FDR across the family of headline tests (listed in notebook 05), reported next
to raw p-values. Final full re-run from a fresh clone.

---

## Addendum D (2026-11-14): written after Addendum C results, before running

Post hoc follow-ups (all labelled post hoc / exploratory in the write-up).

### D1. In which cell types is the hypoxia program coupled to FLT1? (Admati/Skarbianskis pseudobulk)
- All 14 cell classes. A donor-class column is used only if `n_cells ≥ 20`; a class is analysed only
  if ≥ 6 donors remain in the donor set and FLT1 has ≥ 10 reads in ≥ 3 of them.
- Hypoxia score: mean z (across included donors of that class) of log2 CPM of Hallmark HYPOXIA genes.
- **Primary:** Spearman(hypoxia score, FLT1) within **early-PE donors** (up to 10). Secondary: all PE donors (up to 17).
- **Specificity control:** in the same class and donors, compute Spearman(hypoxia score, gene) for every
  expressed gene not in the hypoxia set; report FLT1's percentile. Coupling is called "FLT1-specific"
  only if FLT1 is in the top 5% *and* rho > 0 with p < 0.05.
- Caveat stated in advance: ambient RNA from syncytiotrophoblast can contaminate other cell classes.
- D1b: per class, hypoxia score PE vs control on all donors, `score ~ PE + gestational age`; BH across classes.

### D2. Two-cohort meta-analysis
Fixed-effect inverse-variance pooling of the hypoxia-only contrasts (GSE203507 adjusted model from A2;
GSE114691 unadjusted model from C1) for: EOPE vs preterm control, and PE (both) vs FGR/IUGR-only.
Report pooled estimate, 95% CI, Cochran's Q and I². Random-effects (DerSimonian–Laird) as a check.

### D3. Literature check (not a computation)
PubMed search to calibrate novelty wording for (a) hypoxia vs glycolysis decomposition in placenta and
(b) cell-type-resolved hypoxia–FLT1 coupling.

---

## Addendum E (2026-11-14): corrections + one new primary question, written before running

### E0. Corrections to earlier code (methods fixes, not new hypotheses)
- **P4** was judged on point estimates only. Now tested directly: per-sample difference
  `glycolysis-only score − hypoxia-only score`, OLS on group, contrast EOPE vs control, two-sided p.
- **Meta-analysis (D2)** derived standard errors from CIs with z = 1.96 although the CIs were
  t-based (discovery df = 19). Now uses the model standard errors saved directly from each fit.
  Note: discovery model is adjusted (GA, sex); replication is unadjusted (no covariates available).
- **FDR summary (C3)** mixed one-sided, two-sided and empirical p-values. Now: two-sided p-values for
  every parametric test in one BH family; empirical random-gene-set p-values reported separately.

### E1. Primary new test: is the hypoxia–FLT1 coupling stronger in PE than in growth restriction?
- **Primary cohort: GSE114691** (power: 40 PE [EO-PE + EO-PE+IUGR] vs 18 IUGR-only). Samples: PE and
  IUGR-only only. Model: `FLT1_vst ~ hyp_z * PE` where `hyp_z` = Hallmark HYPOXIA score (FLT1 not in set)
  standardised across these samples. Statistic: the `hyp_z:PE` interaction, two-sided, 95% CI;
  HC3 robust SE as a check. Interpretation fixed in advance: interaction > 0 with p < 0.05 → coupling
  stronger in PE; CI including 0 and excluding effects larger than the PE slope → similar coupling;
  otherwise → inconclusive.
- **Secondary: GSE203507**, PE (15) + FGR-only (5), same model + gestational age + fetal sex.
  Underpowered (5 FGR); reported for direction only.
- Also reported: Spearman(hypoxia, FLT1) separately within PE and within growth restriction in each cohort.

### E2. Robustness (not separate hypotheses)
- **Same genes in both cohorts:** recompute hypoxia-only, shared, glycolysis-only and hypoxia (all)
  scores using only genes measured in both cohorts; repeat the A2/C1 headline contrasts and E1.
- **Leave-one-out:** E1 interaction (each of 58 samples dropped); within-PE Spearman(hypoxia, FLT1) in
  discovery (15), replication (40) and single-cell STB early-PE donors (10). Report range and any sign flip.
- **Cell-composition alternative:** add a syncytiotrophoblast marker score (mean z of the marker-table
  syncytiotrophoblast genes) as a covariate in E1.
