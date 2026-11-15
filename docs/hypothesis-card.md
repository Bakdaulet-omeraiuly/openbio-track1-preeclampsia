# Hypothesis Card

## Title

A trophoblast low-oxygen (HIF-type) program, not a separate glycolytic switch, marks a subset of early-onset preeclampsia placentas in two cohorts and tracks FLT1, the gene for the sFlt-1 biomarker; it is stronger in, but probably not exclusive to, preeclampsia compared with growth restriction alone.

## Clinical/Biological Problem

Early-onset preeclampsia (EOPE) is thought to start with a poorly perfused placenta, but placentas from EOPE are heterogeneous, and the same pregnancies are often complicated by fetal growth restriction (FGR) and early delivery. Ackerman et al. (2023) reported a "glycolytic" placental subphenotype in a subset of EOPE ± FGR cases in this cohort. It is unclear whether that signal (i) survives gestational-age-appropriate controls, (ii) is specific to preeclampsia rather than to growth restriction, (iii) is a metabolic program in its own right or part of a hypoxia response, and (iv) which cells carry it.

## Dataset(s) and Relevant Cells

- GSE203507 villous bulk RNA-seq, 30 donors (8 EOPE, 7 EOPE+FGR, 5 FGR-only, 5 preterm controls [iPTB], 5 term controls), one sample per donor.
- GSE203507 decidua basalis bulk RNA-seq (27 donors, same pregnancies).
- Admati/Skarbianskis single-cell pseudobulk (independent cohort, 26 donors): villous cytotrophoblast and syncytiotrophoblast, chosen in advance.
- Zeisel/Admati gene browser for FLT1, LEP, HK2, NDRG1, P4HA1.
- MSigDB Hallmark HYPOXIA and GLYCOLYSIS (v2024.1).
- Replication: GSE114691 villous bulk RNA-seq, 79 donors (20 EO-PE, 20 EO-PE+IUGR, 18 EO-IUGR, 21 preterm-labour controls; Awamleh et al. 2019, [doi:10.1186/s12920-019-0548-x](https://doi.org/10.1186/s12920-019-0548-x)).

## Observed Result

1. **PE vs preterm controls** (design `~ fetal sex + PE`): 424 genes pass padj < 0.05, |log2FC| ≥ 1; 207 of them survive a gestational-age covariate and every leave-one-control-out refit. Known PE genes behave as expected (FLT1 log2FC 3.0, LEP 6.0, ENG 2.1).
2. **Hypoxia score** is higher in EOPE than in iPTB after adjusting for gestational age and fetal sex: +0.54 (95% CI 0.24 to 0.84). In a **direct** contrast it is also higher in PE than in FGR-only placentas (+0.35, 0.08 to 0.63; HC3 p = 0.018). FGR-only is not different from iPTB (+0.09, −0.24 to 0.42).
3. **The "glycolysis" signal is the hypoxia signal.** Splitting the two Hallmark sets: genes shared by both sets (+0.76) and hypoxia-only genes (+0.42) are higher in EOPE, but glycolysis-only genes are not (+0.07, −0.19 to 0.32). Hypoxia-only and glycolysis-only scores are uncorrelated across placentas (rho = −0.14).
4. **Not an artifact of many changed genes**: the observed effects exceed all 1,000 expression-matched random gene sets (empirical p ≤ 0.001 for each score and contrast).
5. **Heterogeneous**: about 6 of 15 PE placentas carry most of the signal; it does not clearly follow blood pressure, birthweight or delivery week within PE (all p > 0.07, n = 15).
6. **FLT1 link**: the hypoxia score correlates with FLT1 within PE (rho = 0.78, p = 0.0006). FLT1 is not in the hypoxia set. The same relation appears in independent single-cell syncytiotrophoblast data (rho = 0.72, p = 0.006, 13 donors).
7. **Cell context**: in the independent pseudobulk data, the hypoxia score is higher in early PE than early controls in both villous cytotrophoblast (+0.27) and syncytiotrophoblast (+0.36); exact permutation p = 0.15 and 0.08 (only 3 control donors), positive in every leave-one-donor-out refit, and p = 0.02 / 0.013 when all 26 donors are modelled with gestational age. Glycolysis scores do not differ. In the gene browser, FLT1 and LEP are concentrated in syncytiotrophoblast (LEP almost only in the SCT-PSG subtype) and EVT and appear higher in early PE; NDRG1 and P4HA1 are broad across trophoblast. Marker genes of no single cell class move together in bulk (fig4).
8. **Not villus-only**: the decidua (maternal side) shows a similar PE effect (+0.48 vs +0.45 in villi; interaction p = 0.45).
9. **Replication in different patients (GSE114691, 79 placentas, preterm-labour controls; predictions registered before running, Addendum C):**
   - P1 replicated: hypoxia-only score higher in EOPE than controls (+0.21, 95% CI 0.06 to 0.37, one-sided p = 0.004), larger than 984/1,000 matched random gene sets (p = 0.017).
   - **P2 not replicated**: PE vs IUGR-only +0.11 (−0.03 to 0.25, one-sided p = 0.059). IUGR-only placentas were themselves modestly elevated (+0.12). By our own pre-registered rule, the "PE-specific versus growth restriction" part of the claim fails in this cohort; the direction is the same.
   - Shared genes replicated for both contrasts; glycolysis-only genes were *lower* in EOPE (−0.20), reinforcing that the "glycolytic" signal is the hypoxia program.
   - FLT1 link replicated: within 40 PE placentas rho = 0.59 (p = 6×10⁻⁵). In the single-cell syncytiotrophoblast data, the link holds within early-PE donors alone (rho = 0.84, n = 10, p = 0.002), so it is not just a PE-vs-control group difference.
   - Across 17 headline tests, Benjamini–Hochberg q < 0.05 for all discovery and replication tests except the PE-vs-IUGR replication (q = 0.08), its random-set check (q = 0.11), and the small single-cell group comparisons (q = 0.11–0.17) (`notebooks/results/C3_fdr_summary.csv`).

10. **Post hoc, exploratory (Addendum D):**
   - *Two-cohort meta-analysis* of the hypoxia-only PE − FGR/IUGR contrast: fixed effect +0.16 (95% CI 0.04 to 0.28), random effects +0.18 (0.00 to 0.35), I² = 40%. Pooled evidence leans toward a modest PE excess, but the discovery cohort is included (winner's curse), so this does not overturn the failed replication test.
   - *Cell-type map* (independent single-cell data, within early-PE donors): the hypoxia score tracks FLT1 with rho > 0.6, p < 0.05 **and** FLT1 ranks in the top 5% of all ~10,000 genes in syncytiotrophoblast (rho 0.84), extravillous trophoblast (0.93), endothelial cells (0.77), smooth muscle (0.68) and Hofbauer cells (0.64). Syncytiotrophoblast and endothelial cells also pass this in all 17 PE donors. Villous cytotrophoblast shows the alarm (PE +0.24) but weaker FLT1 coupling (rho 0.41); however, a paired bootstrap does not show the STB−VCT difference to be reliable (difference 0.44, 95% CI −0.13 to 1.24). Low-FLT1 immune classes may reflect ambient RNA from syncytium.
   - The PE-vs-control rise in hypoxia score is largest in the three trophoblast classes (EVT +0.31, STB +0.30, VCT +0.24; nominal p < 0.05, BH q = 0.07).

11. **Is the hypoxia → FLT1 coupling stronger in PE than in growth restriction? (Addendum E, primary test pre-specified)**
   - GSE114691 (40 PE vs 18 IUGR-only), `FLT1 ~ hypoxia × PE`: interaction −0.31 (95% CI −0.92 to 0.29, p = 0.30). Within-group Spearman: PE 0.59, IUGR-only 0.58. Pre-specified verdict: **coupling similar, not stronger in PE**. Adding a syncytiotrophoblast marker score made the interaction negative (−0.61, p = 0.02), i.e. if anything weaker in PE. No sign flip in leave-one-out (range −0.81 to −0.18).
   - Discovery cohort (15 PE vs 5 FGR): same direction (−1.79, p = 0.09), underpowered.
   - *Descriptive, not pre-specified:* in GSE114691, at the same hypoxia score PE placentas have higher FLT1 (offset +1.07 on the VST scale, 95% CI 0.50 to 1.63; +0.62 after the STB adjustment; not converted to a fold change). This offset was **not** seen in the small discovery cohort (−0.46, −2.32 to 1.40).
   - Robustness: using only genes measured in both cohorts left all headline contrasts unchanged; within-PE hypoxia–FLT1 correlations never changed sign when any single donor was dropped (bulk 15 and 40 donors, single-cell 10 donors).
   - *Methods corrections (Addendum E0):* P4 now has a direct test (glycolysis-only minus hypoxia-only, EOPE vs control: −0.42, p < 10⁻⁵); meta-analysis uses model standard errors (PE vs growth restriction pooled +0.16, 0.04 to 0.28); the BH summary now uses two-sided p-values throughout.

12. **FLT1 co-expression not captured by the Hallmark hypoxia score (Addendum F; module built in GSE114691 → tested in GSE203507, predictions pre-specified; Addendum H selection-process null):**
   - In GSE114691 we removed the hypoxia score from FLT1 and from every gene and froze the 50 genes whose residuals best tracked FLT1's residual. The module is the "canonical" PE placental gene set: LEP, ENG, PAPPA2, HTRA4, INHBA, FSTL3, PROCR, EBI3, HTRA1, ADAM12, SH3PXD2A, CYP11A1, SPAG4, among others.
   - **Tested in different patients (GSE203507; not used to build the module, but previously analysed, so not an untouched holdout):** V1 partial Spearman(module, FLT1 | hypoxia) = 0.89 across 25 preterm placentas (p < 10⁻⁷); V3 above all 1,000 expression-matched random 50-gene sets (95th percentile 0.49, p = 0.001); V2 module higher in PE than FGR-only after adjusting for hypoxia, gestational age and sex (+0.59, one-sided p = 0.024).
   - Post hoc checks: unchanged after removing pseudogenes/antisense (0.87), after adjusting for EVT and STB marker scores (0.82), and within PE only (0.93, n = 15). The module does not track EVT markers (rho −0.22).
   - In single-cell data, 56% of module genes peak in extravillous trophoblast and 18% in syncytiotrophoblast (FLT1 also peaks in EVT). Within early-PE syncytiotrophoblast the partial correlation was weak (0.33, p = 0.20, n = 10), so the cellular source is unresolved.
   - **Selection-process null (Addendum H, pre-specified):** repeating the whole procedure for 300 random anchor genes matched to FLT1's expression level gave a median replication partial rho of 0.63 (95th percentile 0.875). FLT1's 0.887 ranked at the **95.7th percentile** (13/300 anchors higher), i.e. just above the pre-set 95% bar. The PE − FGR module difference was **not** FLT1-specific (76th percentile; anchor median +0.30). So much of the V1/V2 signal reflects general co-expression reproducibility and broad PE shifts; FLT1's module replicates only modestly better than typical.
   - **Caveat:** "not captured by the hypoxia score" refers only to the generic Hallmark HIF signature. Several module genes (e.g. LEP, ENG) have been reported as hypoxia-responsive in trophoblast, so a trophoblast-specific oxygen response cannot be excluded.

13. **Does placental stress translate into a stronger maternal-side response in PE? (Addendum G, paired villous + decidua from 27 pregnancies; not supported)**
   - G1 (primary, pre-specified): villous hypoxia score vs the same pregnancy's decidual inflammatory score (TNFα/NF-κB + inflammatory-response genes, no overlap with hypoxia genes), group-centred and adjusted for decidual trophoblast content: rho −0.05; within-group permutation p = 0.61 (22 preterm pairs). Real donor pairing carried no information, so we do not claim cross-tissue coupling.
   - G2 (underpowered): PE vs FGR-only decidual response at similar villous stress −0.11 (−0.47 to 0.25). Villous stress ranges of PE and FGR-only barely overlap (6/17 pairs in common support), so this data cannot answer the question. Descriptively, the two highest decidual inflammation scores were FGR-only pregnancies.
   - Angiogenesis response (secondary): rho −0.41, opposite to the pre-specified direction (one-sided p = 0.97).

14. **Is the FLT1 module more than generic co-expression or normal differentiation? (Addendum I, exploratory)**
   - *I2, selection repeated on shuffled data:* permuting FLT1's hypoxia- and group-adjusted residual within clinical group and re-selecting 50 genes each time gave a median test-cohort partial rho of −0.14 (95th percentile 0.44); the real residual gave 0.84 (p = 0.002). Donor bootstrap 95% CI for the frozen module 0.60–0.96; leave-one-out 0.86–0.91. Together with Addendum H: the replication is not a selection artefact, but modules built around other real genes also replicate well, and FLT1's is near (not far beyond) the top.
   - *I3, normal CTB → EVT differentiation (GSE294617, 3 first-trimester donors, in vitro):* the module rose in all three donors (+1.8 to +2.0 log2 units, larger than all 1,000 matched random sets), as did FLT1 (+3.2 to +6.9 log2). **The module is part of the normal trophoblast differentiation program.** A symbol audit recovered `HLA.G` as `HLA-G`, increasing both bulk differentiation scores from 89 to 90 genes; the remaining exclusions are recorded as absent from the dataset vocabulary or removed by expression filtering. The corrected score still did not explain the module–FLT1 association (partial rho 0.887 → 0.887 and 0.939 → 0.924). By the pre-specified rule, differentiation state is not the leading explanation of the bulk association, with the caveat that an early, in-vitro reference may not capture third-trimester differentiation.
   - *I1 (Ackerman VT Blue/Black comparison; post-review baseline correction):* 31 of 44 module genes measured in GSE203507 were assigned to VT Blue (70%; hypergeometric p = 3.3×10⁻²⁸), versus 3 to VT Black (7%; p = 0.082). The corrected primary baseline uses the full published scores excluding only FLT1. In GSE203507 the frozen module and full Blue score correlated at rho = 0.936, crossing the pre-set |rho| > 0.9 high-collinearity caution threshold; the nested-model ΔR² = 0.053 (p = 1.4×10⁻⁴) therefore does not establish a distinct biological contribution. This threshold is not a formal identifiability test. In GSE114691 rho = 0.802 and ΔR² = 0.123 (p = 8.1×10⁻²⁰; 0.111 with an STB marker score), but this is descriptive because that cohort selected the module. The original overlap-free scores are retained only as a sensitivity. **Reading:** the list is predominantly a compact VT Blue-type signal; reproducible information distinct from the full Blue program remains unresolved.
   - *I4, regulators (not tested):* transcription factors FOXM1 and CEBPB regulate sFLT1 during CTB → EVT differentiation (Yong et al. 2025, [doi:10.1093/molehr/gaaf031](https://doi.org/10.1093/molehr/gaaf031)); HIF-1α can be induced without hypoxia by AT1 autoantibodies or LIGHT (Iriyama et al. 2015). These are candidate explanations for the co-expression, not findings of ours.

## Proposed Mechanism

In a subset of early-onset PE pregnancies, the trophoblast experiences low or fluctuating oxygen (poor spiral-artery remodelling and malperfusion), stabilising HIF and switching on a canonical hypoxia program. Part of that program is glycolytic enzyme genes (robust hits in our data: HK2, LDHA, TPI1, P4HA1, ERO1A), which is why the subtype looks "glycolytic". The same HIF activation drives FLT1 transcription in syncytiotrophoblast (FLT1 has a HIF-responsive enhancer; Gerber et al. 1997), contributing to the excess sFlt-1 that causes maternal endothelial dysfunction (Maynard et al. 2003). Growth restriction without preeclampsia engages the same program, and converts it to FLT1 at a similar rate (Addendum E), consistent with shared placental malperfusion. What may distinguish preeclampsia is therefore **how much** hypoxia signalling there is, plus FLT1 variation that the generic hypoxia score does not explain and that FLT1 shares with the canonical PE gene module (LEP, ENG, PAPPA2, HTRA4, INHBA, FSTL3 …) in two cohorts (Addendum F). Its replication is only modestly better than for modules built around random genes, and its PE > FGR shift is not FLT1-specific (Addendum H), so it is a lead about co-regulation, not evidence of a separate mechanism.

**Working model (speculative):** placental FLT1 variation = (i) a component that follows the generic hypoxia/HIF signature and is shared by PE and growth restriction + (ii) a component shared with the canonical PE genes that this signature does not capture. Prediction: in trophoblast cultured at normal oxygen, perturbing a driver of component (ii) should change FLT1 together with LEP/ENG/PAPPA2 without changing Hallmark hypoxia genes; HIF knockdown should lower component (i) but leave (ii) largely intact.

## Why the Result Is Biologically Plausible

- Placental ischemia/hypoxia is a long-standing model of PE onset, and excess placental sFlt-1 causes PE-like disease in rats (Maynard et al. 2003, [doi:10.1172/JCI17189](https://doi.org/10.1172/JCI17189)).
- FLT1, but not KDR, is transcriptionally induced by hypoxia through a HIF site (Gerber et al. 1997, [doi:10.1074/jbc.272.38.23659](https://doi.org/10.1074/jbc.272.38.23659)).
- Single-cell work shows cell-autonomous FLT1/PGF dysregulation in the syncytium in early, not late, PE (Admati et al. 2023, [doi:10.1016/j.medj.2023.07.005](https://doi.org/10.1016/j.medj.2023.07.005)).
- The serum sFlt-1:PlGF ratio is a validated clinical test for suspected PE (Zeisler et al. 2016, [doi:10.1056/NEJMoa1414838](https://doi.org/10.1056/NEJMoa1414838)).
- The same bulk cohort showed a hypoxia/metabolic module in a subset of EOPE ± FGR (Ackerman et al. 2023, [doi:10.1161/HYPERTENSIONAHA.122.20807](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807)); our contribution is the robustness, specificity and decomposition tests above.
- A very recent exploratory study (one EOPE-FGR vs one control placenta by scRNA-seq) frames EOPE-FGR changes as "glycolysis-related" and prioritises FLT1 in EVT (Zhang et al. 2026, [doi:10.2174/0115665232535282261029093348](https://doi.org/10.2174/0115665232535282261029093348)). Our decomposition across two cohorts suggests such glycolysis signals are better read as part of the hypoxia program.
- In a 321-placenta multi-omics study, a hypoxia-induced microRNA (miR-210) was the hub of the FGR + hypertension network together with FLT1, LEP and HTRA4 (Piekos et al. 2025, [doi:10.1038/s42003-025-08631-6](https://doi.org/10.1038/s42003-025-08631-6)), consistent with our hypoxia–FLT1 coupling.

- Support for a non-hypoxic FLT1/ENG component: in trophoblast cell lines, soluble endoglin release did not increase with hypoxia while sFlt-1 did (Zhao et al. 2020, [doi:10.1096/fba.2020-00057](https://doi.org/10.1096/fba.2020-00057)), consistent with ENG sitting in our module that the generic hypoxia score does not explain. Angiotensin-receptor autoantibodies and the cytokine LIGHT raise placental FLT1 without hypoxia, but through HIF-1α (Iriyama et al. 2015, [doi:10.1161/HYPERTENSIONAHA.115.05314](https://doi.org/10.1161/HYPERTENSIONAHA.115.05314)); our module is not explained by the generic HIF target signature, which would be compatible with a different route, but expression data cannot show this.

**What appears new (based on a targeted, not systematic, PubMed search on 2026-11-14):** (1) separating Hallmark hypoxia and glycolysis genes into shared / hypoxia-only / glycolysis-only and showing, in two cohorts, that glycolysis-only genes do not rise; (2) a pre-registered cross-cohort test of PE-vs-FGR specificity; (3) donor-level, cell-type-resolved hypoxia–FLT1 coupling with a genome-wide specificity ranking; **Not new:** that FLT1 travels with canonical PE genes separately from a HIF/glucose-metabolism program. Ackerman et al. 2023 already described this split in the same cohort: a large VT Blue module containing FLT1, ENG and LEP, and a separate VT Black module enriched for HIF signalling and glucose metabolism, noting the Black module's "independence from the larger VT Blue cluster" and its conservation in GSE114691. Our Addendum F/H module re-derives the Blue-type signal by a different route; we report it as a confirmation, not a discovery. We claim these as contributions of this analysis, not as established discoveries.

## Alternative Explanations/Confounders

| Alternative | Evidence for/against in our data |
|---|---|
| Delivery timing (gestational age) | Preterm controls used; GA covariate; term controls do not drive PCA; effect survives. Labour/mode of delivery not available in bulk metadata and could add acute hypoxia signal. |
| Cell composition (more syncytium in PE) | Campbell et al. 2023 ([doi:10.1038/s42003-023-04623-6](https://doi.org/10.1038/s42003-023-04623-6)) found composition explains ~35–38% of FLT1/LEP/ENG overexpression. Our markers do not move together, the effect persists within single trophoblast classes (pseudobulk), and adjusting for estimated STB/CTB proportions leaves it unchanged (exploratory; deconvolution estimates are noisy). Not fully excluded. |
| Growth restriction rather than PE | Discovery: direct PE − FGR-only contrast positive, FGR-only ≈ preterm controls. Replication: same direction but not significant, and IUGR-only modestly elevated. **Partly supported alternative**: the program is shared with FGR to some degree. |
| Labour in controls (replication) | GSE114691 controls are spontaneous preterm labour, which can itself induce hypoxia genes; this would bias toward the null. Per-sample delivery mode is not public, so we could not adjust. |
| Gene-set artifact | Beats 1,000 expression-matched random sets. |
| Signal is maternal/systemic, not trophoblast-specific | Decidua shows a similar effect; decidua basalis contains invading extravillous trophoblast, so this does not separate the two. |
| Single influential donors | Leave-one-out in bulk and pseudobulk; direction stable. |

## Proposed Validation Experiment or Follow-up Analysis

1. **Further cohorts (computational):** we already did this for GSE114691 (above). Next, a third cohort with per-sample gestational age and delivery mode (e.g. GSE148241, also used by Ackerman et al.), with labour included as a covariate, and a pooled PE − FGR estimate across cohorts to decide whether the PE-specific part is real but small or absent.
2. **Tissue validation:** in the same placentas, immunostain HIF-1α/HIF-2α and FLT1 in syncytiotrophoblast, and measure maternal serum sFlt-1:PlGF; prediction: "hypoxia-high" placentas have higher nuclear HIF in syncytium and higher maternal sFlt-1.
3. **Causal test in vitro:** primary trophoblast or trophoblast stem-cell-derived syncytium cultured at 1–2% vs 8% O2, ± HIF knockdown; prediction: hypoxia raises the hypoxia-only and shared genes and FLT1, and HIF knockdown abolishes both, while glycolysis-only genes change little.
4. **Falsification:** the hypothesis is wrong if, in an independent cohort with preterm controls, the hypoxia-only score is not higher in EOPE, or is equally high in FGR-only, or does not track FLT1.

## Potential Intervention Relevance

If confirmed, a placental hypoxia score could define a "hypoxic" EOPE subgroup expected to have the highest sFlt-1:PlGF ratios, i.e. patients most likely to benefit from sFlt-1-lowering approaches (e.g. apheresis trials) or from future HIF-pathway targeting. This is a hypothesis for stratification, not a treatment recommendation.

## Confidence Level and Limitations

- **Medium–high** for "hypoxia program higher in a subset of EOPE than in preterm controls" and "the glycolytic signal is the hypoxia program" (replicated in a second cohort).
- **Medium–high** for the hypoxia–FLT1 coupling within PE (three datasets).
- **Low–medium** for specificity versus growth restriction (discovery yes, replication same direction but not significant).
- **Low–medium** for the FLT1 co-expression module being special (passes the selection-process null only at the 95.7th percentile; PE > FGR shift not FLT1-specific; 'beyond hypoxia' refers only to the generic HIF signature).
- **Low–medium** for trophoblast localisation and the HIF mechanism.

- Small groups (5 preterm controls, 5 FGR-only; 3 early controls in the single-cell data).
- Bulk tissue mixes cell types; deconvolution across platforms is approximate.
- Gestational age and labour/delivery mode can affect hypoxia genes; only partly controlled.
- Expression shows association, not cause; no protein, oxygen or HIF measurements.
- Several analyses (Addenda A–C) were planned after seeing earlier results; they are labelled as such in `docs/analysis-plan.md`.
- The replication cohort lacks per-sample covariates (GA, sex, delivery mode) and was also used by Ackerman et al. for module discovery, so it is "different patients", not a fully untouched validation set.
- Gene-browser plots are a different study and are used as context only.
