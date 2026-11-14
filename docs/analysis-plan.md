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
