# 404: Team Not Found — A Low-Oxygen Alarm in the Early-Preeclampsia Placenta

OpenBio 2026 Datathon, Track 1: Maternal-Fetal Cell Communication in Preeclampsia

![Headline](figures/headline.png)

## Team

Team name: **404: Team Not Found** (one-person team)

| Name | Role |
|---|---|
| Bakdaulet Omertayev | Everything: analysis, pre-registration, article, video |

## Research Question

The question started from the "glycolytic" signal in early-onset preeclampsia (EOPE) placentas. Is it real
against preterm controls, is it specific to preeclampsia, and which cells carry it?

It became a test of each link in the textbook chain:

malperfusion → placental hypoxia program → HIF → FLT1 (sFlt-1) → maternal disease.

## Results at a Glance

Every row below was a pre-registered test, with its plan committed before the data were opened, unless marked
*post hoc*.

| Link | Result | Where |
|---|---|---|
| "Glycolytic" signal = hypoxia program | Hypoxia-only genes up in EOPE vs preterm controls; glycolysis-only genes not; beats 1,000 random sets | notebooks 01–03 |
| EOPE > GA-matched preterm controls | 5 of 6 cohorts; **fragile** in the pooled small cohorts (Addendum M) | 05, K, M |
| Specific to PE vs growth restriction? | **No** (−0.19 SD [−0.56, 0.19], 42 FGR placentas) | L |
| Hypoxia ↔ FLT1 within PE | **10/10 cohorts**, 9 pre-registered; pooled 0.47–0.69 (**headline figure**) | C, K, L, M, P |
| HIF-2 program ↔ FLT1 beyond generic hypoxia | 0.39 [0.19, 0.57] in 3 untouched cohorts; replicated (0.40) | L, P |
| Within one placenta: single cells | 17/19 donors, pure trophoblast barcodes, above all random sets | R |
| Within one placenta: intact tissue | 18/19 donors, composition-adjusted, 99.9th percentile | S |
| HIF binds the program (ChIP) | 24% vs 7% expected; **FLT1 itself not bound** | T |
| HIF-2 vs HIF-1 (iPSC syncytium, 3% O2) | EPAS1-null: FLT1 −5.5 log2, program down; HIF1A-null: neither (confirms Sasagawa 2018 for FLT1) | V |
| Malperfusion → hypoxia program | **Not supported** (MVM n = 168; infarct, ELGAN; Doppler). FLT1 *is* up with malperfusion lesions (pre-registered in ELGAN) | N, W, X |
| Placental FLT1 → maternal severity | **Not supported** (BP, proteinuria, HELLP, PE vs hypertension alone) | Y |

Article: [`docs/article.md`](docs/article.md) (short). Full log: [`docs/supplementary-results.md`](docs/supplementary-results.md).
Hypothesis: [`docs/hypothesis-card.md`](docs/hypothesis-card.md).

## Data Used

- [x] GSE203507 bulk RNA-seq, villous tissue, 30 donors (required)
- [x] Zeisel/Admati gene browser: FLT1, LEP, HK2, NDRG1, P4HA1 (required) → [`notebooks/results/browser_lookups.csv`](notebooks/results/browser_lookups.csv)
- [x] Admati/Skarbianskis single-cell pseudobulk files (optional extension)
- [x] GSE203507 decidua basalis bulk RNA-seq (optional extension)
- [x] Other public data: GSE114691 villous RNA-seq, 79 placentas ([GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)), used for replication
- [x] Other public data: GSE294617 CTB→EVT differentiation reference (CPM; Yong et al. 2025), used descriptively
- [x] Other public data: Ackerman et al. 2023 Table S4 VT Blue/VT Black co-expression module assignments
- [x] Other public data (Addenda K–N): GSE75010 (157 placentas), GSE35574, GSE24129, GSE147776, GSE44711, GSE74341, GSE66273 (array cohorts), GSE216799 (PT2385 vs DMSO in hypoxic cytotrophoblasts; Colson et al. 2023). Downloaded by `artifacts/fetch_external_cohorts.sh`
- [x] Other public data (Addenda O–Y): GSE80995 (ARNT shRNA), GSE186257, GSE324038 (HIF KO T cells), GSE319236 (placental scRNA-seq + Visium), GSE209936 (ARNT ChIP-seq), GSE288288 + GSE302997 (EPAS1/HIF1A-null iPSC syncytium), GSE303463 (MVM pathology), GSE154829 (ELGAN, infarcts). Each addendum's script downloads its own files (`fetch` command) and records SHA-256 in `inputs.sha256`.
- [x] Other public data: MSigDB Hallmark gene sets v2024.1 (HYPOXIA, GLYCOLYSIS), downloaded by the notebook from the [Broad data server](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/)

## How to Run

```bash
git clone https://github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia.git
cd openbio-track1-preeclampsia
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt jupyter
cd notebooks
jupyter nbconvert --to notebook --execute --inplace 01_metabolic_stress_analysis.ipynb
jupyter nbconvert --to notebook --execute --inplace 02_followup_pseudobulk_and_specificity.ipynb
jupyter nbconvert --to notebook --execute --inplace 03_stress_tests.ipynb
jupyter nbconvert --to notebook --execute --inplace 05_external_replication.ipynb
jupyter nbconvert --to notebook --execute --inplace 06_headline_figure.ipynb
jupyter nbconvert --to notebook --execute --inplace 07_celltype_coupling_and_meta.ipynb
jupyter nbconvert --to notebook --execute --inplace 08_coupling_PE_vs_FGR.ipynb
jupyter nbconvert --to notebook --execute --inplace 09_hypoxia_independent_flt1_program.ipynb
jupyter nbconvert --to notebook --execute --inplace 10_villous_to_decidua.ipynb
jupyter nbconvert --to notebook --execute --inplace 11_anchor_gene_null.ipynb
jupyter nbconvert --to notebook --execute --inplace 12_selection_null_and_differentiation.ipynb
```

The first notebook clones the official datathon data package at a **pinned commit** (`fb751e9`) next to this repository and downloads the Hallmark gene sets, so no data files need to be copied by hand. Run the notebooks in order: later notebooks read result files written by earlier ones. Notebook 05 also downloads GSE114691 counts from GEO and a transcript-to-gene table from Ensembl GRCh37 BioMart.

**External cohorts (Addenda K and L).** These scripts read public GEO files that are not in the datathon package:

```bash
bash artifacts/fetch_external_cohorts.sh external_cohorts   # ~270 MB on disk; SHA-256 in artifacts/external_cohorts.sha256
python artifacts/addendum-k/k_confirm_gse75010.py external_cohorts/GSE75010_series_matrix.txt external_cohorts/GPL6244.txt
python artifacts/addendum-l/l1_hif2_program.py external_cohorts/GSE216799_raw_counts.tsv.gz external_cohorts/GSE216799_series_matrix.txt.gz
python artifacts/addendum-l/l2_screen.py external_cohorts/GSE75010_series_matrix.txt external_cohorts/GPL6244.txt
python artifacts/addendum-l/l3_confirm.py external_cohorts                 # pre-specified confirmation
python artifacts/addendum-l/l4_calibration_posthoc.py external_cohorts     # post hoc checks (l4, l5, l6)
python artifacts/addendum-l/l5_stb_adjustment_posthoc.py external_cohorts
python artifacts/addendum-l/l6_deconvolution_posthoc.py external_cohorts
```

**Addenda R–Y (single-cell, spatial, ChIP, knockouts, pathology, maternal disease).** Each script has a
`fetch DIR` and a `run DIR` step; inputs are checksummed in its folder:

```bash
python artifacts/addendum-r/r_within_donor.py fetch r_data && python artifacts/addendum-r/r_within_donor.py run r_data
python artifacts/addendum-r/r3_purity_posthoc.py r_data                      # post hoc purity check
python artifacts/addendum-s/s_spatial.py fetch s_data && python artifacts/addendum-s/s_spatial.py run s_data
python artifacts/addendum-t/t_chip.py fetch t_data && python artifacts/addendum-t/t_chip.py run t_data
python artifacts/addendum-u/u0_make_lists.py t_data && python artifacts/addendum-u/u_direct_targets.py s_data
python artifacts/addendum-v/v_isoform.py fetch v_data && python artifacts/addendum-v/v_isoform.py run v_data
python artifacts/addendum-w/w_mvm.py fetch w_data && python artifacts/addendum-w/w_mvm.py run w_data
python artifacts/addendum-x/x_infarct.py fetch x_data && python artifacts/addendum-x/x_infarct.py run x_data
python artifacts/addendum-y/y_maternal.py external_cohorts w_data
python artifacts/summary-figures/make_figures.py                             # fig17-fig20 and headline.png
```

On 2026-11-15 we re-downloaded every file. Each one matched the checksums of the files we had analysed. Re-running K and L1–L6 on the fresh files reproduced every committed result file byte for byte. The FLT1 isoform analysis (J1/J2) re-processes raw SRA reads and needs HISAT2, samtools and sra-tools. Its steps are in `artifacts/flt1-isoform-feasibility/`.

Reproducibility details:

- Notebook files: `notebooks/01_…` to `12_…`, run in order (about 10–15 minutes in total on a laptop)
- Where we ran it: our own computer (macOS, Python 3.13, packages pinned in `requirements.txt`)
- Random seeds: `numpy` seed 42 (random gene sets, notebook 03); PyDESeq2 is deterministic
- Analysis plan committed **before** tests were run; follow-up analyses registered as Addenda A–Y before they were run (see git history of `docs/analysis-plan.md`)

## Files

| Path | Contents |
|---|---|
| `notebooks/01_metabolic_stress_analysis.ipynb` | Pre-registered analysis: QC, PCA, PyDESeq2 PE vs preterm controls, sensitivity analyses S1–S4, hypoxia/glycolysis scores, marker table |
| `notebooks/02_followup_pseudobulk_and_specificity.ipynb` | Addendum A: independent single-cell (pseudobulk) check, PE vs FGR-only direct contrasts, shared vs set-specific genes |
| `notebooks/03_stress_tests.ipynb` | Addendum B: random-gene-set null, FLT1 link, heterogeneity within PE, decidua, composition adjustment |
| `notebooks/05_external_replication.ipynb` | Addendum C: pre-specified replication in GSE114691 (79 placentas), FLT1 within PE in single-cell data, FDR summary |
| `notebooks/06_headline_figure.ipynb` | Per-placenta hypoxia scores in two cohorts (`figures/fig0_group_scores.png`) |
| `notebooks/07_celltype_coupling_and_meta.ipynb` | Addendum D (post hoc): hypoxia–FLT1 coupling in 14 cell classes, two-cohort meta-analysis |
| `notebooks/08_coupling_PE_vs_FGR.ipynb` | Addendum E: pre-specified test of whether hypoxia–FLT1 coupling differs between PE and growth restriction; same-gene and leave-one-out robustness |
| `notebooks/09_hypoxia_independent_flt1_program.ipynb` | Addendum F: FLT1 co-expression module beyond the Hallmark hypoxia score, built in GSE114691 and tested in GSE203507 |
| `notebooks/10_villous_to_decidua.ipynb` | Addendum G: paired villous stress vs decidual (maternal-side) inflammatory response; negative result |
| `notebooks/11_anchor_gene_null.ipynb` | Addendum H: selection-process null (same procedure for 300 random anchor genes) |
| `notebooks/12_selection_null_and_differentiation.ipynb` | Addendum I: comparison with Ackerman VT Blue/Black modules; residual-permutation selection null; normal CTB→EVT differentiation reference (GSE294617) |
| `notebooks/results/` | All result tables (DE results, robust hits, scores, contrasts, browser lookups) |
| `figures/` | Exported graphs. `headline.png` (= fig17, hypoxia–FLT1 in 10 cohorts) and fig18–fig20 are made by `artifacts/summary-figures/make_figures.py` |
| `docs/article.md` | Science communication article (short, ~2,000 words) |
| `docs/supplementary-results.md` | Full analysis log (the long article), every planned and unplanned test |
| `docs/analysis-plan.md` | Comparison, model, thresholds, and addenda |
| `docs/hypothesis-card.md` | Our hypothesis and how to test it |
| `docs/novelty-audit.md` | What is known, what is replication, what may be new (targeted PubMed search) |
| `docs/presentation-outline.md` | Video script |
| `artifacts/addendum-k/` | Addendum K: seven-program mechanism screen and one-time GSE75010 confirmation |
| `artifacts/addendum-l/` | Addendum L: HIF-2 program from GSE216799, screen, one-time confirmation in GSE35574/GSE24129/GSE147776, post hoc calibration and composition checks (`results.md`) |
| `artifacts/addendum-m/` … `artifacts/addendum-y/` | One folder per later addendum: `plan.md` (committed before data), script, `results.md`, result JSON, input checksums |
| `artifacts/summary-figures/` | Script for the summary figures (fig17–fig20, headline) |
| `artifacts/fetch_external_cohorts.sh` | Downloads the GEO files used by Addenda K and L |
| `artifacts/flt1-isoform-feasibility/` | FLT1 transcript/exon audit, processed-count diagnostics, SRA metadata, raw-reanalysis plan (Addendum J1), cohort marker counts, coverage gate and the pre-specified PE–FGR isoform contrast (null) |

## Key Findings

**How to read the novelty** (`docs/novelty-audit.md`): the HIF-2 → FLT1 biology is known. Our contributions are (1) a pre-registered, multi-cohort replication of the hypoxia–FLT1 association, with its weak points reported, and (2) linking a perturbation-defined trophoblast HIF-2 program to placental FLT1 in clinical cohorts (Addendum L), which a targeted PubMed search did not find done before. This is not a new PE mechanism, and the program is not PE-specific.


- Comparing 15 EOPE placentas with 5 gestational-age-appropriate preterm controls (`~ fetal sex + PE`), 424 genes pass padj < 0.05 and |log2FC| ≥ 1; 207 survive a gestational-age covariate and every leave-one-control-out refit.
- Discovery: a hypoxia score is higher in EOPE than in preterm controls (+0.54, 95% CI 0.24–0.84) and, in a direct contrast, than in FGR-only placentas (+0.35, 0.08–0.63). The "glycolysis" signal comes from genes shared with the hypoxia set; glycolysis-only genes do not rise in PE. The effect exceeds all 1,000 expression-matched random gene sets.
- **Replication in 79 different patients (GSE114691), with predictions registered before running:** the hypoxia-only score is again higher in EOPE than preterm controls (+0.21, 0.06–0.37; beats 98% of random gene sets), glycolysis-only genes are lower, and the score tracks FLT1 within 40 PE placentas (rho = 0.59). **Not replicated:** PE vs growth-restriction-only (+0.11, −0.03 to 0.25); growth-restricted placentas show a weaker version of the signal.
- The hypoxia score tracks FLT1 (sFlt-1 gene; not in the set) within PE (rho = 0.78) and in independent single-cell syncytiotrophoblast (rho = 0.72). In the independent pseudobulk data the hypoxia score is higher in early PE in both cytotrophoblast and syncytiotrophoblast (exploratory, 3 control donors). The gene browser places FLT1 and LEP mainly in syncytiotrophoblast/EVT.

- Post hoc (Addendum D): pooling both cohorts gives a modest PE excess over growth restriction (+0.16, 0.04–0.28; I² = 40%). In the single-cell data, the hypoxia score tracks FLT1 specifically (FLT1 in the top 5% of all genes) in syncytiotrophoblast, extravillous trophoblast and endothelial cells.
- Pre-specified (Addendum E): the hypoxia → FLT1 coupling is **similar** in PE and growth restriction (GSE114691 interaction −0.31, 95% CI −0.92 to 0.29; within-group rho 0.59 vs 0.58). Descriptively, in GSE114691 PE placentas carry higher FLT1 at the same hypoxia score (+1.07 on the VST scale, not a fold change); this offset was not seen in the small discovery cohort. Robust to using identical genes in both cohorts and to leave-one-out.
- **FLT1 co-expression beyond the Hallmark hypoxia score (Addenda F, H and I):** a 50-gene module built in GSE114691 and tested in GSE203507 (previously analysed, not used for building) tracks FLT1 after removing the hypoxia score (partial rho = 0.89). It is the canonical PE gene set (LEP, ENG, PAPPA2, HTRA4, INHBA, FSTL3…). Direct comparison with Ackerman Table S4 places 31/44 measured module genes in VT Blue (hypergeometric p = 3.3×10⁻²⁸) and only 3/44 in VT Black. A post-review correction used the full published scores, excluding only FLT1: in GSE203507 the module and VT Blue score crossed the pre-set high-collinearity caution threshold (rho = 0.94), so the significant nested model (ΔR² = 0.053) cannot establish a distinct biological contribution. In the module-building cohort the module added ΔR² = 0.123, but this is descriptive; the original overlap-free analysis is retained only as a sensitivity. Against a selection-process null (same procedure for 300 random genes of similar expression) FLT1 ranks at the 95.7th percentile, just above the pre-set 95% bar; its PE > FGR shift is not FLT1-specific (76th percentile). Interpreted as reproducible co-expression and a compact Blue-type signal, not a separate mechanism or a new discovery.
- **Negative (Addendum G):** in 22 paired preterm pregnancies, villous hypoxia did not predict the same pregnancy's decidual inflammatory response beyond shuffled pairs (p = 0.61); PE vs FGR-only could not be answered (little overlap in stress levels, 5 FGR).
- **HIF-2 program and three untouched cohorts (Addendum L):** we defined the 100 genes most lowered by the HIF-2 inhibitor PT2385 in hypoxic trophoblasts (GSE216799), screened them in our three used cohorts, then tested them once in GSE35574, GSE24129 and GSE147776 (138 placentas).
  - **Supported:** the program tracks FLT1 beyond the Hallmark hypoxia score (pooled partial rho 0.39 [0.19, 0.57], Holm p = 5e-4; above 1,000/1,000 matched random sets, post hoc; unchanged after post hoc syncytiotrophoblast-marker (0.39 [0.18, 0.57]) and single-cell deconvolution (0.43 [0.21, 0.60]) adjustments; in GSE24129, n = 16, full composition adjustment turned that cohort's estimate negative). Hypoxia tracks FLT1 within PE (rho 0.62).
  - **Not supported:** PE vs control (Holm p = 0.09; heterogeneous). PE vs FGR-only showed no clear difference (−0.19 SD [−0.56, 0.19]; 42 FGR-only).
  - Figure: `figures/fig16_hif2_program_flt1.png` (script `artifacts/addendum-l/l7_figure.py`). Details: `artifacts/addendum-l/results.md`.
- **Third cohort (Addendum K, GSE75010, 157 placentas, microarray; opened once after the plan was committed):** the hypoxia-only score is again higher in early PE than in <34-week preterm non-PE, non-FGR deliveries: +0.29 unadjusted for labour; +0.23 [0.06, 0.41] adjusted for GA, sex and attempted labour, Holm p = 0.013. Labour is nearly separated by group (8/49 vs 24/27), so the adjustment is not decisive. The hypoxia–FLT1 link within PE holds (rho = 0.58). PE vs FGR-only shows no clear difference (0.04; 14 FGR-only, underpowered). A seven-program screen winner (ROS pathway) did not replicate. Details: `artifacts/addendum-k/`.
- **Matched-preterm-control replication (Addendum M, three untouched cohorts, run once; plan `0b12ad5`):** **M1 supported:** hyp_only PE vs preterm controls +0.67 SD [0.06, 1.28], Holm p = 0.033. **Fragile:** it leans on GSE44711 (+1.05); without it +0.29 [−0.58, 1.16]; without the GSE66273 HELLP cases +0.29 [−0.16, 0.74]; the CI reaches or crosses 0 in 4 of 6 pre-specified sensitivities. GSE66273 (chorionic plate) +0.07; 5 of 6 preterm-control cohorts clearly positive. **M2 supported:** hypoxia–FLT1 within PE 0.47 [0.05, 0.75], Holm p = 0.031 (sign positive in 9/9 bulk cohorts, including the discovery cohort and a weak 0.26 in GSE66273; pooled estimates 0.47–0.62). M3 HIF2dep partial rho 0.45 is descriptive and PE-only, not a replication of L-C4. Post hoc descriptive: placental FLT1:PGF log2 ratio PE − control +1.67 [0.51, 2.83] (villous); cohort ordering is the inverse of M1 (GSE44711 +0.59, GSE66273 +3.44), so it does not support the hypoxia result. Details: `artifacts/addendum-m/results.md`.
- **Cell-level link (Addendum R, pre-registered, GSE319236 scRNA-seq, 20 untouched donors):** **supported.** Within donor and within STB, hypoxia tracks FLT1: marker-selected barcodes ρ 0.137, 19/20 donors, sign-flip p ≤ 1×10⁻⁴, above all 1,000 decile-matched random sets; 0.135 after a post hoc STB-state adjustment. A post hoc purity check (rule committed first, `d4d2c18`) found a median 33% of those barcodes also carried a lineage marker and only 65% were KRT7+; **restricted to marker-pure trophoblast barcodes the link stands at ρ 0.099, 17/19 donors, p = 6×10⁻⁴, above all random sets**, so contamination inflated but did not create it. Descriptive: HIF2dep | hypoxia 0.217 (18/19); FLT1 higher in PE STB (p = 6×10⁻⁴) but mean hypoxia not (p = 0.19). Caveats: dissociation stress is a real alternative (FLT1 is hypoxia-inducible); one cohort; the PE-vs-control split is not interpretable (FLT1 detected in 37% vs 16% of pure barcodes). Details: `artifacts/addendum-r/results.md`.
- **Intact-tissue link (Addendum S, pre-registered, GSE319236 Visium, 19 donors, 46,453 never-dissociated spots):**
  - **Supported.** Within donor, hypoxia tracks FLT1 across spots after adjusting for spot composition (STB,
    endothelium, stroma, immune, erythrocyte markers) and technical covariates.
  - ρ 0.118, 18/19 donors, sign-flip p ≤ 1×10⁻⁴, 99.9th percentile of 1,000 random sets.
  - All 3 sensitivities pass, including the STB-rich half of spots (0.078).
  - **This rules out dissociation stress as the explanation for Addendum R.**
  - Limits: spots mix several cells, so composition is adjusted for, not removed; same donors as R.
  - Details: `artifacts/addendum-s/results.md`.
- **Direct HIF binding (Addendum T, pre-registered, GSE209936 ARNT ChIP-seq, BeWo at 1% O2):**
  - The positive control passed: hypoxia genes are bound 34.5% vs 7.5% expected.
  - **The HIF2dep program is enriched for ARNT binding:** 24.4% vs 7.3%, 3.3×, p = 0.001, at all windows. The
    program therefore includes direct HIF targets.
  - **FLT1 is not bound** (nearest peak 123 kb), so FLT1's link to hypoxia in trophoblast may be indirect.
  - HIF isoform (1 vs 2) is not resolvable with ARNT.
  - Details: `artifacts/addendum-t/results.md`.
- **Direct HIF targets in intact tissue (Addendum U, pre-registered, GSE319236 Visium):**
  - **Supported.** The 20 ARNT-bound program genes track FLT1 within donor, beyond hypoxia, composition and
    technical covariates: ρ 0.155, 18/19 donors, above all random sets.
  - The unbound genes track it too (0.125, 19/19), so the route (direct target vs differentiation) is not
    separated.
  - Details: `artifacts/addendum-u/results.md`.
- **HIF-2 vs HIF-1 (Addendum V, pre-registered, JAX iPSC-derived syncytium, GSE288288 + GSE302997, 3% O2):**
  - **HIF-2-preferring, a genetic double dissociation.**
  - EPAS1-null: FLT1 −5.46 log2 (p = 4×10⁻⁵), program −1.28 SD (p = 3×10⁻⁵), Hallmark hypoxia unchanged.
  - HIF1A-null: Hallmark hypoxia −0.39 (p = 0.009), FLT1 and program not lowered.
  - Not a differentiation failure: EPAS1-null cells are more STB-like, and FLT1 is still −6.2 after adjustment.
  - For FLT1 this independently confirms Sasagawa et al. 2018 (siRNA); the program-level double dissociation is new.
  - Limits: one parental iPSC line (KOLF2.2J), so replicates are clones, not donors; an early embryonic-like model; two batches; HIF1A n = 3.
  - Details: `artifacts/addendum-v/results.md`.
- **Malperfusion → hypoxia (Addendum W, pre-registered, GSE303463, 168 placentas, 49 with MVM):**
  - **Not supported.** MVM does not raise the generic hypoxia score beyond disease, GA, sex and batch: +0.12 SD,
    p = 0.25, 7.9th percentile of random sets. Within hypertensive disorders the difference is −0.01.
  - Descriptive: FLT1 (+0.61) and the HIF2dep program (+0.46) are higher with MVM.
  - Post hoc calibration (committed first): HIF2dep is above all 1,000 random sets, and FLT1 is at the 98.6th percentile of expression-matched genes (rank 121/17,128). MVM tracks the HIF-2/FLT1 axis, not generic hypoxia. A lead, not confirmed.
  - Details: `artifacts/addendum-w/results.md`.
- **Replication of the MVM lead (Addendum X, pre-registered, ELGAN GSE154829, singletons, 40 infarcted vs 163):**
  - **FLT1 replicates:** +0.42 SD, Holm p = 0.011, above all random genes.
  - **The HIF2dep program does not** (+0.03, p = 0.43). Generic hypoxia is flat.
  - Malperfusion lesions link to FLT1 in two cohorts, but to no hypoxia program.
  - Details: `artifacts/addendum-x/results.md`.
- **Placental FLT1 vs maternal disease (Addendum Y, pre-registered, link 6): not supported.**
  - Within PE (GSE75010), FLT1 vs BP, proteinuria and HELLP is positive but not significant after Holm.
  - PE vs hypertension alone (GSE303463): +0.32 SD, p = 0.093.
  - FLT1 is already raised in hypertension alone (+0.68, p = 0.0007).
  - Details: `artifacts/addendum-y/results.md`.
- **Regulator and cell-state tests (Addendum O, pre-registered):** **O1 supported.** In GSE80995, ARNT (HIF-1β) shRNA in first-trimester CTB lowered the PT2385-defined program at 2% O2 (Δ −1.15 z, 5/5 donors, exact sign-flip p = 0.031, the minimum for 5 pairs; below all 1,000 random sets). The two independent perturbations agree. Every drug off-target effect is not excluded, since ARNT partners (HIF-1α, AhR) and differentiation are not separated. HIF-2 specificity is not shown. The program also fell at 20% O2 (−0.52 as labelled; −0.85 under the alternative, title-based labels, a sensitivity analysis). There was no consistent FLT1 fall after knockdown. **O2 not supported** (second session): within-STB pseudobulk partial rho 0.26 [−0.14, 0.59], so cell mix is not excluded. Details: `artifacts/addendum-o/results.md`, `artifacts/addendum-o2/results.md`.
- **Fourth untouched cohort (Addendum P, GSE186257, RNA-seq, plan `d2be780`):** HIF2dep–FLT1 partial rho | hyp_all 0.40 [0.06, 0.65], Holm p = 0.011 (replicates L-C4). With STB + EVT marker adjustment it is 0.63 [0.35, 0.80]. Hypoxia–FLT1 within PE is 0.69. **Specificity not shown:** the unadjusted result is at the 83rd percentile of matched random sets (random median +0.23). Post hoc, after the marker adjustment it is above all random sets (median 0.00). Within PE only: 0.26 (not clear). Details: `artifacts/addendum-p/results.md`.
- **HIF-2 vs HIF-1 specificity (Addendum Q, GSE324038, hypoxic human T cells):** **not shown in T cells.** In trophoblast it was later tested directly in Addendum V (HIF-2-preferring). The knockouts worked (HIF1A −1.19, EPAS1 −0.95) and Hallmark hypoxia was HIF-1-dominated. But the 59 expressed program genes fell with neither knockout (D = −0.09, 57th percentile of random sets). FLT1 fell with EPAS1-KO (−0.66) and rose with HIF1A-KO (+0.36), descriptively. Details: `artifacts/addendum-q/results.md`.
- **Clinical anchors (Addendum N, exploratory, GSE75010):** the hypoxia score did not track uterine-artery Doppler PI (rho −0.01, n = 52) or, within PE, systolic BP, newborn or placental weight z, or HELLP (all Holm p ≥ 0.38). Late-onset PE was also raised vs ≥34-week controls (+0.16 [0.02, 0.31]); early > late was not clear (+0.17 [−0.03, 0.37]). Details: `artifacts/addendum-n/results.md`.
- **FLT1 transcript forms (Addenda J/J1): pre-specified null.** We re-processed the raw reads of all 79
  GSE114691 placentas: candidate extraction, then full-genome HISAT2, counting unique reads only. We then
  counted the junction unique to sFlt1-14 and the two junctions unique to canonical Flt1. A coverage gate
  ran before any group comparison, without group labels, and passed: median 1,042 sFlt1-14 junction reads
  and 38–52 canonical junction reads per donor. The model adjusted for total FLT1 and the hypoxia score.
  The sFlt1-14 : Flt1 RNA ratio did not differ between PE (n=40) and growth-restriction-only (n=18)
  placentas: −0.06 log2, 95% CI −0.52 to 0.40, p=0.79. Both disease groups have higher ratios than
  controls descriptively. A pre-specified GSE203507 replication (Addendum J2, disease vs GA-matched preterm controls) was not run because its blinded coverage gate failed (body-vs-junction rho 0.68 < 0.70), so this remains untested. Details: `artifacts/flt1-isoform-feasibility/feasibility-audit.md`.

## Limitations

5 preterm controls and 5 FGR-only placentas; 3 early control donors in the single-cell data. Bulk tissue mixes cell types and our deconvolution is approximate. Labour and delivery mode (possible acute hypoxia) are not recorded for the bulk samples. A similar signal in decidua means we cannot call it villus-specific. The replication cohort has no public per-sample covariates (GA, sex, delivery mode), its controls laboured, and it was also used by Ackerman et al. for module discovery. Expression data show association, not cause. Addenda A–I were planned after seeing earlier results and are labelled as such. Gene-browser plots are from a different study and are used as context only; we read the condition legend (early/late, control/PE) from the browser's legend bars.

## Presentation

Video link: [add your unlisted YouTube link]

## Citations

- Ackerman WE et al. (2023) Hypertension 80:1363–1374. [doi:10.1161/HYPERTENSIONAHA.122.20807](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807) (GSE203507)
- Admati I, Skarbianskis N et al. (2023) Med 4:687–709. [doi:10.1016/j.medj.2023.07.005](https://doi.org/10.1016/j.medj.2023.07.005) (single-cell data and gene browser)
- Awamleh Z, Gloor GB, Han VKM (2019) BMC Med Genomics 12:91. [doi:10.1186/s12920-019-0548-x](https://doi.org/10.1186/s12920-019-0548-x) (GSE114691)
- Zhang R et al. (2026) Curr Gene Ther. [doi:10.2174/0115665232535282261029093348](https://doi.org/10.2174/0115665232535282261029093348)
- Piekos SN et al. (2025) Commun Biol 8:1239. [doi:10.1038/s42003-025-08631-6](https://doi.org/10.1038/s42003-025-08631-6)
- Leavey K et al. (2016) Hypertension 68:137–47. [doi:10.1161/HYPERTENSIONAHA.116.07293](https://doi.org/10.1161/HYPERTENSIONAHA.116.07293)
- Zhao J et al. (2020) FASEB BioAdv 3:23–35. [doi:10.1096/fba.2020-00057](https://doi.org/10.1096/fba.2020-00057)
- Iriyama T et al. (2015) Hypertension 65:1307–15. [doi:10.1161/HYPERTENSIONAHA.115.05314](https://doi.org/10.1161/HYPERTENSIONAHA.115.05314)
- Yong Q et al. (2025) Mol Hum Reprod 31(3). [doi:10.1093/molehr/gaaf031](https://doi.org/10.1093/molehr/gaaf031) (GSE294617)
- Campbell KA et al. (2023) Commun Biol 6:264. [doi:10.1038/s42003-023-04623-6](https://doi.org/10.1038/s42003-023-04623-6)
- Gerber HP et al. (1997) J Biol Chem 272:23659–67. [doi:10.1074/jbc.272.38.23659](https://doi.org/10.1074/jbc.272.38.23659)
- Maynard SE et al. (2003) J Clin Invest 111:649–58. [doi:10.1172/JCI17189](https://doi.org/10.1172/JCI17189)
- Zeisler H et al. (2016) N Engl J Med 374:13–22. [doi:10.1056/NEJMoa1414838](https://doi.org/10.1056/NEJMoa1414838)
- Liberzon A et al. (2015) The Molecular Signatures Database Hallmark gene set collection. Cell Syst 1:417–425. [doi:10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004)
- Muzellec B et al. (2023) PyDESeq2. Bioinformatics 39:btad547. [doi:10.1093/bioinformatics/btad547](https://doi.org/10.1093/bioinformatics/btad547)
