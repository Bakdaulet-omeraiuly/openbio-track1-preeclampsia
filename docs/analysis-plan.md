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

---

## Addendum F (2026-11-14): discovery → validation of a hypoxia-independent FLT1 program, written before running

Motivation (post hoc, from E1): in GSE114691, PE placentas carried ~2× more FLT1 than growth-restricted
placentas at the same hypoxia score. Question: which gene program accompanies FLT1 *beyond* hypoxia,
and does it validate in an independent cohort?

### F1. Discovery (GSE114691, all 79 samples; nothing from GSE203507 is used)
- Hypoxia score = Hallmark HYPOXIA (all), as before. For FLT1 and for every expressed gene g, take the
  residual from `expression ~ hypoxia score`. Spearman correlation of each gene's residual with the FLT1
  residual (partial correlation given hypoxia).
- **Module** = the 50 genes with the highest partial correlation, excluding FLT1 itself and any
  Hallmark HYPOXIA gene. Frozen and saved to `results/F1_module_genes.csv` before validation.
- Descriptive: overlap of the module with Hallmark sets (hypergeometric, BH) and with the known PE genes
  LEP, ENG, INHA, PAPPA2, HTRA4, PGF.

### F2. Validation in GSE203507 (independent patients; 25 preterm samples, TB excluded)
Module score = mean z (VST, z across the 30 samples) of module genes measured in GSE203507.
- **V1 (primary):** partial Spearman(module score, FLT1 | hypoxia score) across the 25 preterm samples > 0
  (one-sided α = 0.05).
- **V2:** `module ~ C(group) + hypoxia score + GA + fetal sex`; contrast mean(PE groups) − FGR-only > 0
  (one-sided α = 0.05).
- **V3 (null check):** V1 exceeds the 95th percentile of 1,000 random 50-gene sets matched by
  expression decile (seed 42).

### F3. Cell-type location (single-cell pseudobulk; descriptive)
Mean CPM of module genes per cell class (all donors); fraction of module genes whose highest-expressing
class is syncytiotrophoblast, EVT or another class. Within early-PE syncytiotrophoblast (10 donors):
partial Spearman(module, FLT1 | hypoxia).

The claim "a hypoxia-independent FLT1 program" is supported only if V1 and V3 pass.
Known risk: the module may simply track syncytiotrophoblast abundance; reported via F3 and by repeating
V1 with the STB marker score as an extra covariate.

---

## Addendum G (2026-11-14): villous stress → decidual response, written before running

Question: at comparable villous hypoxia-response levels, do PE pregnancies show a stronger decidual
inflammatory/vascular response than FGR-only pregnancies, and does a donor's villous stress carry
information about *the same donor's* decidual response?

Data: GSE203507 paired villous + decidua basalis samples (27 donors: EOPE 8, EOPE+FGR 4, FGR 5, iPTB 5,
TB 5). Single cohort, no replication available → exploratory. Decidua basalis contains invading
extravillous trophoblast, so trophoblast contamination is modelled explicitly.

Scores (fixed now; no gene selection from results):
- **Stress (villi):** Hallmark HYPOXIA score in villous tissue (as in notebooks 01–03).
- **Response (decidua), primary:** Hallmark TNFA_SIGNALING_VIA_NFKB ∪ INFLAMMATORY_RESPONSE, **minus every
  Hallmark HYPOXIA gene**. Secondary: Hallmark ANGIOGENESIS minus HYPOXIA genes ("vascular").
  Decidual scores = mean z of VST across the 27 decidua samples.
- **Trophoblast content (decidua):** mean z of marker-table genes for villous cytotrophoblast,
  syncytiotrophoblast and extravillous trophoblast, measured in decidua.

### G1 (primary): is the real donor pairing informative?
Preterm pairs only (n = 22; TB excluded). Statistic: Spearman correlation between villous stress and
decidual response after removing group means from both (group-centred) and regressing the decidual
response on decidual trophoblast content. Null: permute decidua samples **within clinical group**
10,000 times (seed 42). One-sided p (coupling > 0). "Cross-tissue coupling" is claimed only if p < 0.05.

### G2 (secondary, underpowered: 12 PE vs 5 FGR-only)
`decidual response ~ PE + villous stress + decidual trophoblast content + GA`, PE + FGR-only pairs.
Report the PE coefficient with 95% CI, and report whether villous-stress ranges of the two groups overlap
(common support). Not used for any confirmatory claim.

### G3 (descriptive)
Same as G1/G2 for the secondary (angiogenesis) response score; leave-one-donor-out range for G1.

---

## Addendum H (2026-11-14): selection-process null for the Addendum F module, written before running

Concern (from review): gene–gene co-expression replicates across cohorts in general, so a module built
around *any* gene might show a high partial correlation in another cohort. The random-gene-set null in V3
does not repeat the selection step. Clarification: GSE203507 was not used to build the module, but it had
been analysed extensively before; it is a previously studied cohort of different patients, not an
unopened holdout.

### H1. Anchor-gene null
- Draw 300 "anchor" genes at random (seed 42) from genes that are (i) not in Hallmark HYPOXIA, (ii) not
  FLT1, (iii) measured in both cohorts, and (iv) in the same expression decile (GSE114691 mean VST) as FLT1.
- For each anchor, repeat Addendum F exactly: in GSE114691, residualise all genes and the anchor on the
  hypoxia score (rank-based), take the 50 non-hypoxia genes with the highest partial correlation with the
  anchor, then compute V1 in GSE203507 (partial Spearman(module, anchor | hypoxia), 25 preterm samples).
- Report FLT1's percentile among the anchors. Interpretation fixed now:
  - FLT1 ≥ 95th percentile → its module replicates unusually well;
  - otherwise → the F1 result mainly reflects general co-expression reproducibility, and the claim is
    limited to "FLT1 sits in a reproducible co-expression module (the canonical PE genes)".
- Also report, for the same anchors, the V2-style PE − FGR-only module contrast (adjusted for hypoxia,
  GA, sex) and FLT1's percentile there (descriptive).

### H2. Wording corrections (no analysis)
- The +1.07 offset in E1 is on the VST scale; it is not reported as a fold change.
- "Hypoxia-independent" is replaced by "not explained by the Hallmark hypoxia score".

---

## Addendum I (2026-11-14): does the FLT1 module add anything beyond known programs and differentiation? Written before running

All analyses are exploratory (both bulk cohorts have been analysed many times). The frozen 50-gene module
from Addendum F (`results/F1_module_genes.csv`) is used unchanged.

### I1. Published Ackerman et al. 2023 modules (contingent on obtaining Table S4)
- Gene lists for VT Blue and VT Black from Table S4; symbols matched to our data as-is; background = genes
  measured in each cohort.
- Overlap of our module with Blue and Black (counts, Jaccard, hypergeometric).
- Scores: Blue score and Black score = mean z of member genes, **excluding FLT1 and excluding our module's
  genes** (so predictors do not share genes with each other or with the outcome).
- Model in each cohort (GSE203507 preterm n = 25; GSE114691 n = 79):
  `FLT1 ~ hypoxia + Blue + Black + group (+ STB marker score)` vs the same + `module`; report ΔR² and the
  nested F-test for adding the module. If |Spearman(module, Blue)| > 0.9 → "not identifiable".
- Stop rule: if the module adds nothing beyond Blue/Black (F-test p ≥ 0.05 in both cohorts), the
  "distinct module" claim is dropped; the module is described as a compact proxy of the published program.

### I1/I3 correction (2026-11-14): documented after external review and after the first I1 run
- The original I1 score definition above is retained and labelled **overlap-free sensitivity**. It removes
  information from the published modules and therefore is not the fairest primary baseline.
- Corrected **primary** scores use every measured VT Blue/VT Black member except FLT1, including genes that
  overlap our module. The same nested models, STB sensitivity and `|rho| > 0.9` non-identifiability rule apply.
  Results from this correction are post-review, not pre-specified.
- Terminology correction after review: the implementation calls `|rho| > 0.9` **high collinearity**, not
  non-identifiability. The original wording above is preserved as history; this threshold is a caution rule,
  not a formal mathematical identifiability test.
- Gene symbols for I3 are audited against each target dataset. Exact matches are preferred; a `.` ↔ `-`
  substitution is accepted only when it gives one unique target symbol. The audit separately records genes
  absent from the raw dataset and genes present in raw data but removed by the expression filter. No broader
  alias is inferred without an explicit mapping source.

### I2. Residual-permutation selection null (selection step repeated)
- Discovery nuisance model (GSE114691): FLT1 and all genes residualised (rank-based) on hypoxia score +
  clinical group. Observed: top-50 non-hypoxia genes by correlation with FLT1's residual (this spec may
  differ slightly from F1, which used hypoxia only; both reported).
- Null: permute FLT1's residual **within clinical group** (500 permutations, seed 42), repeat selection,
  freeze, evaluate V1 in GSE203507 (partial Spearman(module, FLT1 | hypoxia), 25 preterm). p = (1 + #null ≥ obs)/501.
- Exchangeability assumption: within a clinical group, FLT1 residuals are exchangeable across donors under H0.
- Separate inference on the frozen F1 module in GSE203507: donor bootstrap 95% CI (2,000, seed 42) and
  leave-one-donor-out range for V1.

### I3. Normal CTB → EVT differentiation reference (GSE294617; 3 donors × day 0/3/6; CPM, not counts)
- Score rule: log2(CPM + 1); module score = mean over module genes present of the gene-centred values
  (centred across the 9 samples). Same for Hallmark hypoxia score and FLT1.
- Report each donor's trajectory and the day-6 minus day-0 change; direction consistency across 3 donors.
  Context: fraction of 1,000 expression-matched random 50-gene sets (seed 42) with a larger mean d6 − d0.
- **Differentiation score:** genes up at day 6 vs day 0 in all three donors, top 100 by mean log2 change,
  **excluding FLT1, our module genes and Hallmark hypoxia genes**. In each bulk cohort, test whether the
  module still tracks FLT1 after adding this score: partial Spearman(module, FLT1 | hypoxia, differentiation).
- Interpretation fixed: if the module rises with normal CTB→EVT differentiation and the partial
  correlation drops substantially (by > 50%) after adjusting for the differentiation score, differentiation
  state is the leading explanation. A weak overlap is **not** evidence of PE specificity (early, in-vitro reference).

### I4. Mechanism (text only)
Candidate regulators are discussed, not tested; prior work (Iriyama 2015 HIF without hypoxia; Yong 2025
FOXM1/CEBPB–sFLT1) is cited and not claimed.

Stop rule for the whole line: whatever I1–I3 show, no further hypotheses are added in this project.

---

## Addendum J (2026-11-14): FLT1 transcript-form feasibility, written before any PE–FGR isoform contrast

This is a separate molecular question from the FLT1 co-expression-module line. No group difference is
tested until the measurement gate below passes.

### J0. Annotation and processed-count audit
- GRCh37 FLT1 transcript structures are compared across every annotated FLT1 transcript; exact protein
  matches to curated UniProt isoforms provide names only, not RNA quantification validity.
- Transcript-specific terminal exon bodies and junctions are listed. A coordinate-distinct exon is not
  considered sequence-specific when its body lies inside another FLT1 exon.
- GSE114691 processed transcript rows are inspected for presence, count structure and correlation, but are
  not treated as isoform fractions because shared-read allocation, multi-mapping and length correction are
  not documented.
- Primary candidates are canonical Flt1 (`ENST00000282397`) and sFlt1-14 (`ENST00000541932`). sFlt1-14 is
  measured by its unique terminal junction; canonical Flt1 uses the mean of two unique internal junctions
  because its terminal splice boundary is shared with i15 and isoform 8. Other forms remain exploratory.
  No `total soluble FLT1` quantity is defined because the annotation does not contain
  exact matches for every curated soluble isoform.

### J1. Raw-read feasibility gate, evaluated without group effects
- Raw source: GSE114691 / SRP148556, 79 single-end runs, approximately 50-nt reads.
- Primary markers must be uniquely mappable against the full GRCh37 genome/transcriptome, not merely unique
  among the five FLT1 transcripts.
- At least 70% of the 58 disease donors must have at least 5 reads for the sFlt1-14 junction and each of the
  two canonical-specific internal junctions, plus at least 20 combined primary junction reads. Neither clinical group may have less than 60%
  usable samples.
- Length-normalized unique-body and junction ratios must correlate at rho ≥ 0.70, or one measure must be
  fixed as primary with a documented technical reason before group effects are examined.
- If the gate fails, the line stops with a measurement-limit conclusion; processed rows are not substituted.

### J2. Conditional primary contrast, run only if J1 passes
- Outcome: `log2((sFlt1-14 terminal-junction reads + 0.5) / (mean canonical-specific internal-junction reads + 0.5))`.
- Population: PE (`PEONLY` + `PEIUGR`, n = 40) versus IUGR-only (n = 18).
- Model: `isoform_log_ratio ~ PE + log2(total_FLT1_CPM + 1) + hypoxia_score`, HC3 robust standard errors;
  two-sided PE coefficient, alpha 0.05.
- Controls provide descriptive context only. Body-read ratios and probabilistic transcript abundance are
  sensitivities. The other three annotated forms are exploratory and use BH correction as one family.
- Interpretation is limited to RNA isoform composition at the resolution of 50-nt single-end data; it is
  not equated with secreted protein abundance. A null result does not prove equal isoform biology.

Implementation details and compute constraints are frozen in
`artifacts/flt1-isoform-feasibility/raw-reanalysis-plan.md`.

### Addendum J1 pointer (2026-11-15)

The local execution protocol (Addendum J1) was committed before any cohort run was processed. It is
recorded in `artifacts/flt1-isoform-feasibility/raw-reanalysis-plan.md`, commit `d90d8ce`. The gate passed,
and the result of the primary contrast was null. See `artifacts/flt1-isoform-feasibility/feasibility-audit.md`.

### Addendum J2 pointer (2026-11-15)

Addendum J2, the GSE203507 disease-versus-preterm-control isoform replication, was committed (`cec7e47`)
before any GSE203507 read was fetched. Its blinded coverage gate failed on criterion 4 (rho 0.677 < 0.70),
so the contrast was not run. See `artifacts/flt1-isoform-feasibility/feasibility-audit.md`.

## Addendum K (2026-11-15, 07:50 EST): mechanism screen in the two used cohorts, then one confirmation in an untouched cohort (GSE75010)

**Status of data at writing.** GSE75010 (Leavey et al. 2016; 157 placentas, Affymetrix HuGene 1.0 ST,
GPL6244) has never been used in this repository. Only its `!Sample_` metadata lines were read, to count
groups:

| Group | n | Detail |
|---|---:|---|
| PE | 80 | 26 with IUGR; 49 delivered <34 wk |
| non-PE | 77 | 14 with IUGR; 27 non-IUGR delivered <34 wk; 15 chorioamnionitis |

No expression value was viewed. The series matrix file sits in a scratch folder and is not opened
beyond the metadata lines until the screen winner is committed.

**Processing of GSE75010 (fixed now).**

- Use the deposited series-matrix values. If the maximum value is above 100, apply log2(x + 1);
  otherwise use the values as they are.
- Map probesets to genes with the GPL6244 `gene_assignment` column, taking the first gene symbol.
- Drop probesets without a symbol, then average the probesets of each gene.
- Score = the mean of per-gene z-scores across all 157 samples, the same rule as notebooks 01/05.
- GA = `ga (week)` + `ga (day)`/7; sex = `infant gender`; labour proxy = `attempted vaginal delivery`.

**Family A: replication of article claims, unchanged.** Gene sets are as in notebook 01: `hyp_only` =
HALLMARK_HYPOXIA − HALLMARK_GLYCOLYSIS, `gly_only` = the reverse, `hyp_all` = HALLMARK_HYPOXIA, all from
MSigDB v2024.1. The covariates are GA + sex + labour proxy (OLS).

| Test | Comparison | Prediction / criterion |
|---|---|---|
| A1 | `hyp_only`, PE delivered <34 wk vs non-PE non-IUGR delivered <34 wk | PE higher; one-sided p < 0.05 |
| A2 | `gly_only`, same contrast as A1 | Not raised: one-sided p(>0) ≥ 0.05, reported descriptively |
| A3 | Spearman(`hyp_all`, FLT1) within all 80 PE | rho > 0; one-sided p < 0.05 |
| A4 | `hyp_only`, all PE (80) vs non-PE IUGR (14) | PE higher; one-sided p < 0.05. This is the test that did not replicate in GSE114691 |

Sensitivity analyses for Family A:

- excluding chorioamnionitis cases;
- dropping the labour covariate.

**Family B: screen, then one confirmation.**

- **Closed candidate list.** Seven sets from MSigDB Hallmark v2024.1, each with HALLMARK_HYPOXIA genes
  removed so that they are not redundant with Family A:
  - UNFOLDED_PROTEIN_RESPONSE
  - REACTIVE_OXYGEN_SPECIES_PATHWAY
  - TNFA_SIGNALING_VIA_NFKB
  - OXIDATIVE_PHOSPHORYLATION
  - MTORC1_SIGNALING
  - P53_PATHWAY
  - COMPLEMENT
- **Screen target.** PE (with or without FGR) versus FGR-only, i.e. what separates PE from growth
  restriction:
  - GSE203507: `score ~ PE + ga_delivery_weeks + C(fetal_sex)` (15 vs 5);
  - GSE114691: `score ~ PE` (40 vs 18), on the same VST/z-score pipeline as notebook 05.
- **Metric.** The t statistic of the PE coefficient in each cohort.
- **Selection rule.** Among candidates whose t has the same sign in both cohorts, choose the one with the
  largest min(|t|). If no candidate has the same sign in both, Family B stops.
- **Reporting.** The screen is reported only as a ranking. These two cohorts have been heavily viewed,
  so screen results are not findings.
- **Confirmation (B1).** In GSE75010, the winner is tested as all PE (80) vs non-PE IUGR (14), with
  covariates GA + sex + labour proxy (OLS). The test is one-sided in the screen's direction.

**Multiplicity.** Holm correction across A1, A3, A4 and B1 at a one-sided alpha of 0.05. A2 and the
sensitivity analyses are descriptive.

**Order.**

1. Commit this plan.
2. Run the screen.
3. Commit the winner.
4. Open the GSE75010 expression values once and run Families A and B in a single script.

**Cut-off.** If the run and the write-up are not finished by 14:00 EST, the plan is pushed as registered
but not run, and the article is not changed.

**Interpretation.**

- A positive A1 or A3 replicates claims already in the article in a third cohort.
- A positive A4 or B1 would be the first replicated PE-versus-FGR difference in this project. It would
  still be called a replication of association, not mechanism. Novelty is not claimed until the cohort's
  own paper (Leavey et al. 2016) has been checked for the same result.
- 14 IUGR-only placentas make a null result for A4 and B1 underpowered, not a refutation.

_Correction (added after the run): the heading time "07:50 EST" is a typo. The git timestamp of commit `167cddc` is authoritative for when this plan was committed; the screen and confirmation commits follow it._

## Addendum L (2026-11-15): HIF-2 program and PE-vs-FGR in untouched cohorts

Background and data search: `artifacts/addendum-l/research-memo.md`.

**Data status at writing.** For GSE216799, GSE35574, GSE24129 and GSE147776, only `!Sample_` metadata
lines and the GPL annotation headers have been read; no expression value has been viewed. All earlier
results stay in the repository unchanged; this addendum adds to them and replaces nothing.

**Step 1: perturbation-defined HIF-2 program (GSE216799; raw counts; 3 placentas × DMSO/PT2385 at 2.5% O2).**

- Donor = the `!Sample_description` prefix (P1/P2/P3). Keep genes with ≥10 counts in ≥3 samples.
- Run PyDESeq2 with design `~ donor + treatment`, contrast PT2385 vs DMSO.
- Define **HIF2dep** as the 100 genes with the most negative Wald statistic. FLT1 and all
  HALLMARK_HYPOXIA genes are excluded before ranking, so HIF2dep cannot overlap the generic hypoxia score.
- Descriptive only (n = 3):
  - FLT1 log2 fold change, with the expectation that it is negative;
  - the paired difference in the F1 module score (the 50 genes in
    `notebooks/results/F1_module_genes.csv`) between PT2385 and DMSO, with z-scores taken over the
    6 samples (tests H2);
  - the fraction of F1 genes with negative fold change, compared with all tested genes.

**Step 2: screen of H1 in the used cohorts (GSE114691, GSE203507 villous, GSE75010).**

- Population: disease placentas only, i.e. PE (with or without FGR) plus FGR-only (58, 20 and 94
  placentas).
- Scores use the same per-gene z-score/mean rule as before (VST for RNA-seq; the Addendum K processing
  for GSE75010).
- Statistic: the partial Spearman correlation between the HIF2dep score and FLT1, after residualising
  both on `hyp_all` with OLS.
- **Screen passes** only if the correlation is positive in all three cohorts *and* the fixed-effect
  Fisher-z pooled one-sided p is below 0.05. Otherwise H1 stops and is reported as not supported.
- These cohorts are heavily viewed, so the screen is not a finding.

**Step 3: confirmation in the untouched cohorts, opened once, in one script.**

*Processing.*

- Use the series-matrix values; apply log2(x+1) if the maximum is above 100.
- Map probes to genes and average per gene. Symbol columns: GPL6102 `Symbol`; GPL6244 the first symbol
  of `gene_assignment`; GPL20844 `GENE_SYMBOL`.
- GSE35574 technical repeats: base id = `Sample id` with trailing letters removed. Arrays sharing a base
  id are averaged.
- Within each cohort: z-score genes, take set means, then standardise each score to SD 1 across all
  samples of that cohort. Coefficients are therefore in within-cohort SD units.

*Groups.* PE = all PE, including PE+IUGR. FGR-only = normotensive IUGR/FGR. Control = control.

*Models within each cohort* (OLS, PE coefficient):

| Cohort | Model |
|---|---|
| GSE35574 | `score ~ PE + GA + C(sex) + C(labour) + C(batch)` |
| GSE24129 | `score ~ PE` |
| GSE147776 | `score ~ PE` |

*Pooling.* Inverse-variance fixed effect across cohorts; DerSimonian–Laird random effects as a sensitivity
analysis.

*Tests.*

| Test | Question | Contrast / statistic | Sidedness |
|---|---|---|---|
| **L-C1** (H3, primary) | PE specificity | `hyp_only`, PE vs FGR-only | Two-sided; earlier cohorts disagreed on direction |
| **L-C2** | Replication | `hyp_only`, PE vs control | One-sided, PE higher |
| **L-C3** | Replication | Spearman(`hyp_all`, FLT1) within PE, pooled by Fisher z | One-sided, positive |
| **L-C4** (H1) | Mechanism | HIF2dep–FLT1 partial Spearman within PE + FGR-only, pooled | One-sided, positive |

L-C4 is run only if the step 2 screen passed.

*Multiplicity.* Holm correction across the tests that are run (alpha 0.05).

*Sensitivity.* GSE35574 without the labour and batch covariates; leave-one-cohort-out pooled estimates.

**Interpretation, fixed now.**

- A null L-C1 with a 95% interval is reported as **no clear difference**, with the interval as the bound.
  It does not establish equality.
- L-C4, if positive, is evidence that a trophoblast HIF-2-dependent program tracks placental FLT1 beyond
  generic hypoxia. It is an association with a perturbation-defined program, not proof of HIF-2
  dependence in vivo.
- RNA is not equated with secreted sFlt-1 protein.

**Cut-off.** If the run and the write-up are not finished by 13:00 EST, the plan is pushed as registered
but not run.

**Addendum L amendment, before step 3 was run.** The GSE35574 technical-repeat rule above, which strips
trailing letters from `Sample id`, does not work for this file: most `Sample id` values are array
barcodes. The repeats are instead marked by a `.rep` suffix on `!Sample_title`. There are 8 such arrays,
which matches the series description. The rule is therefore: base id = `!Sample_title` without `.rep`,
giving 86 placentas. This was found from metadata lines only, before any expression value was read.

### Addendum L outcome pointer

Results, post hoc calibration and disclosures: `artifacts/addendum-l/results.md`.

## Addendum M (2026-11-15): hypoxia program against GA-matched preterm controls in three untouched cohorts

The plan and the confirmation script were committed before any expression value was downloaded (`0b12ad5`):
`artifacts/addendum-m/plan.md` and `artifacts/addendum-m/m_confirm.py`. The cohorts are GSE44711,
GSE74341 and GSE66273 (chorionic plate, with a pre-specified leave-it-out sensitivity).

- **Confirmatory, Holm across the two:**
  - M1: hyp_only, PE vs preterm controls.
  - M2: hyp_all–FLT1 rho within PE.
- **Descriptive:** M3, the HIF2dep partial rho within PE.

**Outcome** (`c3df572`): M1 and M2 were supported (M1 leans on GSE44711). FLT1:PGF was added post hoc as descriptive only. See `artifacts/addendum-m/results.md`.

## Addendum N (2026-11-15): clinical anchors in GSE75010, exploratory

Plan and script: `artifacts/addendum-n/` (`494bde1`, committed before running). Results:
`artifacts/addendum-n/results.md`. The Doppler and severity associations were null.

## Addendum O (2026-11-15): separating the HIF2dep–FLT1 link from cell state, and testing regulator dependence

- **O1** (`artifacts/addendum-o/`; plan and script committed before running, `a823f57`). Primary: GSE80995
  ARNT shRNA in first-trimester CTB lowered HIF2dep at 2% O2 (5/5 donors, p = 0.031, beyond all random sets).
  This is supported. A sample-label amendment changes only the normoxia numbers.
- **O2** (`artifacts/addendum-o2/`, run by the second session; plan `8abfcc1`). Primary: within-STB
  pseudobulk partial rho 0.26 [−0.14, 0.59]. This is not supported.

## Addendum P (2026-11-15): GSE186257, fourth untouched cohort

The plan and script were committed before running (`d2be780`); the cohort is GSE186257.

- **P1 (replicates L-C4):** 0.40 [0.06, 0.65], Holm p = 0.011.
- **P2:** 0.69.
- **P1 adjusted for STB/EVT markers:** 0.63.
- **Specificity not shown:** the unadjusted result sits at the 83rd percentile of random sets.

Details: `artifacts/addendum-p/results.md`.
