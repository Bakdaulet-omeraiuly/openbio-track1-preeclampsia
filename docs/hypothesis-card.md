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

## Proposed Mechanism

In a subset of early-onset PE pregnancies, the trophoblast experiences low or fluctuating oxygen (poor spiral-artery remodelling and malperfusion), stabilising HIF and switching on a canonical hypoxia program. Part of that program is glycolytic enzyme genes (robust hits in our data: HK2, LDHA, TPI1, P4HA1, ERO1A), which is why the subtype looks "glycolytic". The same HIF activation drives FLT1 transcription in syncytiotrophoblast (FLT1 has a HIF-responsive enhancer; Gerber et al. 1997), contributing to the excess sFlt-1 that causes maternal endothelial dysfunction (Maynard et al. 2003). Growth restriction without preeclampsia may engage a weaker version of the same program (as in the replication cohort), consistent with shared placental malperfusion; what would distinguish preeclampsia is the strength of the response and its coupling to FLT1/sFlt-1 output.

## Why the Result Is Biologically Plausible

- Placental ischemia/hypoxia is a long-standing model of PE onset, and excess placental sFlt-1 causes PE-like disease in rats (Maynard et al. 2003, [doi:10.1172/JCI17189](https://doi.org/10.1172/JCI17189)).
- FLT1, but not KDR, is transcriptionally induced by hypoxia through a HIF site (Gerber et al. 1997, [doi:10.1074/jbc.272.38.23659](https://doi.org/10.1074/jbc.272.38.23659)).
- Single-cell work shows cell-autonomous FLT1/PGF dysregulation in the syncytium in early, not late, PE (Admati et al. 2023, [doi:10.1016/j.medj.2023.07.005](https://doi.org/10.1016/j.medj.2023.07.005)).
- The serum sFlt-1:PlGF ratio is a validated clinical test for suspected PE (Zeisler et al. 2016, [doi:10.1056/NEJMoa1414838](https://doi.org/10.1056/NEJMoa1414838)).
- The same bulk cohort showed a hypoxia/metabolic module in a subset of EOPE ± FGR (Ackerman et al. 2023, [doi:10.1161/HYPERTENSIONAHA.122.20807](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807)); our contribution is the robustness, specificity and decomposition tests above.

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
- **Low–medium** for trophoblast localisation and the HIF mechanism.

- Small groups (5 preterm controls, 5 FGR-only; 3 early controls in the single-cell data).
- Bulk tissue mixes cell types; deconvolution across platforms is approximate.
- Gestational age and labour/delivery mode can affect hypoxia genes; only partly controlled.
- Expression shows association, not cause; no protein, oxygen or HIF measurements.
- Several analyses (Addenda A–C) were planned after seeing earlier results; they are labelled as such in `docs/analysis-plan.md`.
- The replication cohort lacks per-sample covariates (GA, sex, delivery mode) and was also used by Ackerman et al. for module discovery, so it is "different patients", not a fully untouched validation set.
- Gene-browser plots are a different study and are used as context only.
