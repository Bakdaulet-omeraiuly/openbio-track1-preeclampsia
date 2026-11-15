# Hypothesis Card

## Title

A placental low-oxygen program, not a separate glycolytic switch, is higher in a subset of early-onset preeclampsia placentas than in gestational-age-matched preterm births. It tracks FLT1 (the sFlt-1 gene) in 10 bulk cohorts and within single trophoblast cells and intact tissue. A HIF-2-dependent sub-program tracks FLT1 beyond it in 4 cohorts, and FLT1 needs HIF-2 in a trophoblast model. It is shared with growth restriction rather than specific to preeclampsia.

## Clinical/Biological Problem

Early-onset preeclampsia (EOPE) is thought to start with a poorly perfused placenta, but placentas from EOPE are heterogeneous, and the same pregnancies are often complicated by fetal growth restriction (FGR) and early delivery. Ackerman et al. (2023) reported a "glycolytic" placental subphenotype in a subset of EOPE ± FGR cases in this cohort. It is unclear whether that signal (i) survives gestational-age-appropriate controls, (ii) is specific to preeclampsia rather than to growth restriction, (iii) is a metabolic program in its own right or part of a hypoxia response, and (iv) which cells carry it.

## Dataset(s) and Relevant Cells

- GSE203507 villous bulk RNA-seq, 30 donors (8 EOPE, 7 EOPE+FGR, 5 FGR-only, 5 preterm controls [iPTB], 5 term controls), one sample per donor.
- GSE203507 decidua basalis bulk RNA-seq (27 donors, same pregnancies).
- Admati/Skarbianskis single-cell pseudobulk (independent cohort, 26 donors): villous cytotrophoblast and syncytiotrophoblast, chosen in advance.
- Zeisel/Admati gene browser for FLT1, LEP, HK2, NDRG1, P4HA1.
- MSigDB Hallmark HYPOXIA and GLYCOLYSIS (v2024.1).
- Replication: GSE114691 villous bulk RNA-seq, 79 donors (20 EO-PE, 20 EO-PE+IUGR, 18 EO-IUGR, 21 preterm-labour controls; Awamleh et al. 2019, [doi:10.1186/s12920-019-0548-x](https://doi.org/10.1186/s12920-019-0548-x)).
- Later pre-registered cohorts and experiments, each with its plan committed before the data were opened:
  - **Bulk patient cohorts:** GSE75010, GSE35574, GSE24129, GSE147776, GSE44711, GSE74341, GSE66273,
    GSE186257 (hypoxia–FLT1, 10 cohorts in total).
  - **Single-cell and spatial:** GSE319236 scRNA-seq and Visium (within-donor, cell-level and intact-tissue
    link).
  - **Cell experiments:**
    - GSE216799: PT2385 HIF-2 inhibitor; defines the HIF2dep program;
    - GSE80995: ARNT shRNA;
    - GSE209936: ARNT ChIP-seq;
    - GSE288288 + GSE302997: EPAS1- and HIF1A-null iPSC syncytium;
    - GSE324038: HIF-KO T cells.
  - **Pathology:** GSE303463 (MVM) and GSE154829 (ELGAN, infarcts).
  - **Maternal severity:** GSE75010 clinical variables.

## Observed Result

The full, numbered list of all 29 results, with every number, is in
[`supplementary-results.md`](supplementary-results.md#appendix-full-observed-results-list-from-the-hypothesis-card).

1. **The "glycolytic" signal is a hypoxia signal.** In GSE203507, hypoxia-only genes are higher in EOPE than in
   preterm controls (+0.54 SD [0.24, 0.84]). Glycolysis-only genes are not (+0.07). The effect beats 1,000
   matched random sets. This split was post hoc.
2. **EOPE vs GA-matched preterm controls:** higher in 5 of 6 cohorts. The pooled small-cohort test (M) is
   supported but fragile.
3. **PE vs growth restriction:** higher in the discovery cohort (+0.30) only. It is not replicated in any later
   cohort (−0.19 SD [−0.56, 0.19], 42 FGR placentas).
4. **Hypoxia ↔ FLT1 within PE:** positive in **10/10 bulk cohorts**, 9 pre-registered, pooled 0.47–0.69. A
   HIF-2-inhibitor-defined program tracks FLT1 beyond generic hypoxia in 4 cohorts (L-C4 0.39, P 0.40).
5. **Within one placenta (GSE319236):**
   - single trophoblast cells: 19/20 donors pre-registered; 17/19 after a post hoc purity filter;
   - intact Visium tissue, composition-adjusted: 18/19.
   - This makes cell mix and dissociation stress unlikely as the whole explanation. These results use the
     generic hypoxia score; the HIF2dep-program versions (R2, S2, U) are descriptive.
6. **Mechanism in culture:**
   - the program falls with PT2385 and with ARNT shRNA;
   - ARNT binds 24% of program genes vs 7% expected, but not FLT1;
   - EPAS1-null iPSC syncytium loses FLT1 (−5.5 log2) and the program; HIF1A-null does not. This confirms
     Sasagawa et al. 2018 for FLT1, in one iPSC line.
7. **Malperfusion:** no rise in generic hypoxia with MVM (n = 168, 49 MVM) or infarcts (ELGAN), and no link with
   Doppler. FLT1 *is* higher with these lesions (pre-registered in ELGAN: +0.42 SD).
8. **Maternal disease:** placental FLT1 does not track BP, proteinuria or HELLP within PE, and does not separate
   PE from hypertension alone. It is already up in hypertension alone.

## Proposed Mechanism

**Current hypothesis (2026-11-15, after Addenda R–Y).** In a subset of early-onset preeclampsia, and in growth
restriction, trophoblast runs a **HIF-2-dependent program** that keeps FLT1 high.

- **Strongly supported:**
  - links 3–5 (see the table below);
  - the link holds within single cells and intact tissue.
- **Not supported in our data:**
  - what switches the program on in patients: malperfusion did not raise the generic hypoxia score, although
    FLT1 was raised with malperfusion lesions;
  - whether placental FLT1 explains maternal disease severity.

The paragraph below is the original mechanism statement, kept for the record.

In a subset of early-onset PE pregnancies, the trophoblast experiences low or fluctuating oxygen (poor spiral-artery remodelling and malperfusion), stabilising HIF and switching on a canonical hypoxia program. Part of that program is glycolytic enzyme genes (robust hits in our data: HK2, LDHA, TPI1, P4HA1, ERO1A), which is why the subtype looks "glycolytic". The same HIF activation drives FLT1 transcription in syncytiotrophoblast (FLT1 has a HIF-responsive enhancer; Gerber et al. 1997), contributing to the excess sFlt-1 that causes maternal endothelial dysfunction (Maynard et al. 2003). Growth restriction without preeclampsia engages the same program, and converts it to FLT1 at a similar rate (Addendum E), consistent with shared placental malperfusion. What may distinguish preeclampsia is therefore **how much** hypoxia signalling there is, plus FLT1 variation that the generic hypoxia score does not explain and that FLT1 shares with the canonical PE gene module (LEP, ENG, PAPPA2, HTRA4, INHBA, FSTL3 …) in two cohorts (Addendum F). Its replication is only modestly better than for modules built around random genes, and its PE > FGR shift is not FLT1-specific (Addendum H), so it is a lead about co-regulation, not evidence of a separate mechanism.

### Evidence for each link in the chain (2026-11-15)

The model above is a chain, and it is only as strong as its weakest link. Each link rests on different
evidence. "Causal" below means a perturbation experiment in culture, run by other groups.

| # | Link | What we have | Type |
|---|---|---|---|
| 1 | Malperfusion → placental hypoxia | **Tested and not supported in our data.** Histological MVM (Addendum W, 168 placentas) does not raise the generic hypoxia score beyond disease and GA (+0.12 SD, p = 0.25), and uterine-artery Doppler does not track it (Addendum N). FLT1 and the HIF2dep program are higher with MVM, beyond random genes (post hoc) | generic hypoxia **not raised** by malperfusion in W or X, nor by Doppler in N; **FLT1 is raised** (W post hoc, X pre-registered); the HIF2dep program is raised in W only | **not supported** for hypoxia; malperfusion → FLT1 supported in two cohorts |
| 2 | Hypoxia program up in EOPE vs GA-matched preterm controls | Discovery, plus pre-registered replication in GSE114691 and GSE75010, plus Addendum M. Clearly positive in 5 of 6 preterm-control cohorts. M1 is **fragile**: its CI reaches or crosses 0 in 4 of 6 pre-specified sensitivities | association, replicated |
| 3 | HIF activity → the HIF2dep program | **Causal in culture, twice and in agreement:** PT2385 lowers it (L1, 3 donors) and ARNT shRNA lowers it (O1, 5/5 donors, p = 0.031, beyond all random sets). ARNT also partners HIF-1α and AhR, and both treatments shift differentiation. **ARNT binds near 24% of program genes vs 7% expected (Addendum T, ChIP-seq), so part of the program is direct** | causal in culture + direct binding; isoform not resolved |
| 4 | **HIF-2 (rather than HIF-1) → FLT1** | **Now supported in a trophoblast model (Addendum V):** EPAS1-null syncytium loses FLT1 (−5.46 log2) and the program, while HIF1A-null does not (genetic double dissociation, 3% O2). Earlier evidence was **conflicting**: PT2385 lowers FLT1 (L1, log2FC −1.93), but ARNT shRNA did **not** lower FLT1 consistently (O1: 2 up, 2 down, 1 flat), and hypoxia itself lowers FLT1 in those first-trimester cells. Addendum Q found **no HIF isoform preference** for the program in hypoxic human T cells (57th percentile of random sets), though FLT1 itself was HIF-2-preferring there (descriptive). Prior literature supports HIF-2α (Sasagawa et al. 2018). **Addendum T: in trophoblast at 1% O2, ARNT does not bind near FLT1 (nearest peak 123 kb), although it binds a quarter of the program, so FLT1 may be an indirect target there** | **HIF-2-dependent in iPSC syncytium (genetic); possibly indirect; not yet shown in term STB** |
| 5 | Hypoxia ↔ FLT1 in placental tissue | **Our strongest result.** Positive in 10/10 bulk cohorts, pre-registered in 9; pooled estimates 0.47–0.69. The HIF2dep–FLT1 partial link replicated in four untouched cohorts (L-C4 0.39; P1 0.40) | association, strongly replicated |
| 5a | …within one cell type, rather than cell mix? | **Supported within donor at cell level, in one independent scRNA-seq cohort (Addendum R, GSE319236).** In marker-pure trophoblast barcodes (KRT7+, negative for 5 lineage markers): ρ 0.099, 17/19 donors, sign-flip p = 6×10⁻⁴, above all 1,000 random sets. Contamination inflated the first estimate (0.137 in marker-selected barcodes) but did not create it. **Dissociation stress made unlikely as the whole explanation (not excluded):** in intact, never-dissociated tissue of the same donors (Addendum S, Visium), the link holds after adjusting for spot composition (ρ 0.118, 18/19, 99.9th percentile of random sets). One cohort, two designs. At pseudobulk level (O2) the link was 0.26 [−0.14, 0.59], not significant. Marker and deconvolution adjustments left the bulk estimate unchanged (L5, L6; P1-cs 0.63) | association, supported in one cohort |
| 6 | sFlt-1 → maternal endothelial disease | Causal step: literature only (Maynard et al. 2003). **Placental end tested (Addendum Y): not supported.** FLT1 does not clearly track BP, proteinuria or HELLP within PE, nor separate PE from hypertension alone (p = 0.093); it is already up with hypertension alone | literature; placental-end association **not supported** |

**Where the chain breaks.** Links 1 and 6 are borrowed from the literature. Link 4 is contradicted by our own
data in one system. So we have a well-replicated **association** (link 5) that also holds
inside single cells of one purified cell type and in intact tissue (link 5a), and a **culture-level causal step** (link 3) — but
these are not joined into a causal chain in human placental tissue, because the step that would join them
(link 4) is the one that failed. Addenda T and U narrow the gap. ARNT (the HIF-β subunit shared by HIF-1 and HIF-2) binds near a quarter of the program genes in
trophoblast, but not near FLT1. In intact placenta those direct targets travel with FLT1 within each donor. That fits
FLT1 being an **indirect** target of the HIF program in trophoblast. It is an association, and it does not tell
direct-target and differentiation routes apart. **Addendum V closes most of link 4 in a trophoblast model:** genetic
loss of HIF-2α, but not HIF-1α, removes FLT1 and the program in hypoxic iPSC syncytium. The chain is now supported
link by link (3, 4, 5, 5a), but in different systems: culture, iPSC model and tissue association. It has not been
tested in a single system in term placenta.

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

**What appears new (based on a targeted, not systematic, PubMed search on 2026-11-14):** (1) separating Hallmark hypoxia and glycolysis genes into shared / hypoxia-only / glycolysis-only and showing, in two cohorts, that glycolysis-only genes do not rise; (2) a pre-registered cross-cohort test of PE-vs-FGR specificity; (3) donor-level, cell-type-resolved hypoxia–FLT1 coupling with a genome-wide specificity ranking; **Not new:** that FLT1 travels with canonical PE genes separately from a HIF/glucose-metabolism program. Ackerman et al. 2023 already described this split in the same cohort: a large VT Blue module containing FLT1, ENG and LEP, and a separate VT Black module enriched for HIF signalling and glucose metabolism, noting the Black module's "independence from the larger VT Blue cluster" and its conservation in GSE114691. Our Addendum F/H module re-derives the Blue-type signal by a different route; we report it as a confirmation, not a discovery. We claim these as contributions of this analysis, not as established discoveries. **Update 2026-11-15 (see `docs/novelty-audit.md`):** the most plausible new contribution is Addendum L. It links a perturbation-defined trophoblast program (genes lowered by the HIF-2 inhibitor PT2385; GSE216799; Colson et al. 2023, doi:10.1161/HYPERTENSIONAHA.122.20739) to FLT1 in clinical placenta cohorts beyond the generic hypoxia score. A targeted PubMed search found no prior study doing this. The HIF-2 → FLT1 biology itself is known (Sasagawa et al. 2018). The multi-cohort hypoxia–FLT1 results are a replication contribution. We do not claim a new PE mechanism or a PE-specific pathway.

## Alternative Explanations/Confounders

| Alternative | Evidence for/against in our data |
|---|---|
| Delivery timing (gestational age) | Preterm controls used; GA covariate; term controls do not drive PCA; effect survives. Labour/mode of delivery not available in bulk metadata and could add acute hypoxia signal. |
| Cell composition (more syncytium in PE) | **Largely addressed for the FLT1 link:** the link holds within donor in pure single trophoblast cells (R) and composition-adjusted intact tissue (S). Earlier evidence: Campbell et al. 2023 ([doi:10.1038/s42003-023-04623-6](https://doi.org/10.1038/s42003-023-04623-6)) found composition explains ~35–38% of FLT1/LEP/ENG overexpression. Our markers do not move together, the effect persists within single trophoblast classes (pseudobulk), and adjusting for estimated STB/CTB proportions leaves it unchanged (exploratory; deconvolution estimates are noisy). Not fully excluded. |
| Growth restriction rather than PE | Discovery: direct PE − FGR-only contrast positive, FGR-only ≈ preterm controls. Replication: same direction but not significant, and IUGR-only modestly elevated. **Partly supported alternative**: the program is shared with FGR to some degree. |
| Labour in controls (replication) | GSE114691 controls are spontaneous preterm labour, which can itself induce hypoxia genes; this would bias toward the null. Per-sample delivery mode is not public, so we could not adjust. |
| Gene-set artifact | Beats 1,000 expression-matched random sets. |
| Dissociation stress (single-cell data) | Addendum S: the within-donor link also holds in never-dissociated tissue (Visium), after adjusting for spot composition. Made unlikely as the whole explanation, though not excluded. |
| Signal is maternal/systemic, not trophoblast-specific | Decidua shows a similar effect; decidua basalis contains invading extravillous trophoblast, so this does not separate the two. |
| Single influential donors | Leave-one-out in bulk and pseudobulk; direction stable. |

## Proposed Validation Experiment or Follow-up Analysis

1. **Further cohorts (computational):** we already did this for GSE114691 (above). Done since: GSE75010 (Addendum K), three untouched cohorts (L), three GA-matched cohorts (M) and GSE186257 (P). GSE148241 was screened and not usable (no preterm controls). Still open: a cohort, with labour included as a covariate, and a pooled PE − FGR estimate across cohorts to decide whether the PE-specific part is real but small or absent.
2. **Tissue validation:** in the same placentas, immunostain HIF-1α/HIF-2α and FLT1 in syncytiotrophoblast, and measure maternal serum sFlt-1:PlGF; prediction: "hypoxia-high" placentas have higher nuclear HIF in syncytium and higher maternal sFlt-1.
3. **Causal test in vitro** (partly done: Addendum V used EPAS1- and HIF1A-null iPSC syncytium at 3% O2; still needed in primary term trophoblast, with sFlt-1 protein measured): primary trophoblast or trophoblast stem-cell-derived syncytium cultured at 1–2% vs 8% O2, ± HIF knockdown; prediction: hypoxia raises the hypoxia-only and shared genes and FLT1, and HIF knockdown abolishes both, while glycolysis-only genes change little.
4. **Falsification.** PE vs FGR-only showed no difference in any replication cohort. No difference is not the same as equality, so the PE-specific part is **not supported**, rather than refuted. The remaining conditions still apply: the hypothesis is wrong if, in an independent cohort with preterm controls, the hypoxia-only score is not higher in EOPE, or is equally high in FGR-only, or does not track FLT1.

## Potential Intervention Relevance

Our data **weaken** the original stratification idea:

- the hypoxia score did not track Doppler, maternal blood pressure or other clinical severity (Addenda N, Y);
- placental FLT1 did not separate PE from hypertension alone.

What remains is mechanistic. HIF-2 inhibitors exist (e.g. PT2385/belzutifan), so HIF-2 is a candidate lever on
placental FLT1, to be tested first in term trophoblast (validation item 3). This is not a treatment
recommendation.

## Confidence Level and Limitations

- **Medium–high** for "hypoxia program higher in a subset of EOPE than in preterm controls" (pre-registered replications in GSE114691 and GSE75010; Addendum M supported but fragile) and "the glycolytic signal is the hypoxia program" (replicated in a second cohort).
- **High** for the hypoxia–FLT1 coupling within PE (pooled estimates 0.47–0.62 in the pre-registered tests; sign count positive in 10/10 bulk cohorts, some weak; registered before running in 9 of them: GSE114691, GSE75010, the three Addendum L and three Addendum M cohorts, and GSE186257 (Addendum P, rho 0.69)).
- **Medium–high** for a trophoblast program lowered by HIF-2 inhibition, and also lowered by ARNT shRNA in an independent experiment (Addendum O1; not proven HIF-2-specific, direct or free of cell state), tracking FLT1 beyond the generic hypoxia score (Addendum L; composition not excluded).
- **Low** for the hypoxia score as a marker of clinical severity or measured perfusion (Addendum N null).
- **Low** for specificity versus growth restriction (discovery yes; not clear in GSE114691, GSE75010 or the 42 FGR-only placentas of Addendum L).
- **Low–medium** for the FLT1 co-expression module being special (passes the selection-process null only at the 95.7th percentile; PE > FGR shift not FLT1-specific; 'beyond hypoxia' refers only to the generic HIF signature).
- **Medium** for trophoblast localisation: the coupling holds within marker-pure trophoblast cells of the same donor in an independent scRNA-seq cohort (Addendum R, 17/19 donors, above all random sets), and in intact tissue of the same donors (Addendum S, Visium, 18/19), so neither cell mix nor dissociation explains it; one cohort only.
- **Medium** for the HIF mechanism and for HIF-2 specificity in trophoblast. There is a genetic double dissociation in iPSC syncytium (Addendum V), and HIF-2 inhibition (L1) agrees. ARNT knockdown (O1) and T cells (Q) gave mixed or no isoform signal, and term STB is untested.

- Small groups (5 preterm controls, 5 FGR-only; 3 early controls in the single-cell data).
- Bulk tissue mixes cell types; deconvolution across platforms is approximate.
- Gestational age and labour/delivery mode can affect hypoxia genes; only partly controlled.
- Expression shows association, not cause; no protein, oxygen or HIF measurements.
- Several analyses (Addenda A–C) were planned after seeing earlier results; they are labelled as such in `docs/analysis-plan.md`.
- The replication cohort lacks per-sample covariates (GA, sex, delivery mode) and was also used by Ackerman et al. for module discovery, so it is "different patients", not a fully untouched validation set.
- Gene-browser plots are a different study and are used as context only.
