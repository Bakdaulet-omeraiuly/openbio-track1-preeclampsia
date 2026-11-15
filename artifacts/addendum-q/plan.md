# Addendum Q: is the PT2385-defined program more HIF-2- than HIF-1-dependent? (written before running)

Written on 2026-11-15 at about 14:05 EST, before any expression value of GSE324038 was read. Only the
header row, the first two gene IDs and the GEO sample titles were viewed. The plan is committed together with
`q_hif1_vs_hif2.py`.

## Why and the main limit

We searched GEO and found no public trophoblast experiment that removes HIF-1α and HIF-2α separately. The
best available test is in another cell type.

[GSE324038](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE324038): primary human T cells
activated and kept at 1% O2, with CRISPR knockout of HIF1A or EPAS1 (HIF-2α) and a non-targeting control.
There are 4 donors (batch 1 = the supplementary RSEM count table: LRS-36, -38, -42, -43).

**This is not trophoblast.** Many program genes are trophoblast-specific and will not be expressed. The test
asks only whether, *for the program genes that T cells express*, loss of HIF-2α lowers them more than loss of
HIF-1α. That is evidence about HIF isoform preference of these genes. It is not evidence about the placenta.

## Processing

- **Counts:** RSEM expected counts. Ensembl IDs (version stripped) are mapped to symbols with
  `data_external/ensembl_gene_id_to_name_GRCh37.87.tsv`; counts are summed per symbol.
- **Samples:** the 12 HIF1A-KO, EPAS1-KO and non-targeting (NT) samples. The "No CRISPR" samples are not used.
- **Expression filter:** keep genes with CPM ≥ 1 in ≥ 4 of the 12 samples.
- **Scores:** log2(CPM + 1), gene z across the 12 samples, then the mean z of each set.
  - HIF2dep: frozen list from Addendum L.
  - hyp_all: HALLMARK_HYPOXIA.

## Tests

Within each donor, Δ1 = HIF1A-KO − NT and Δ2 = EPAS1-KO − NT.

| Item | Statistic | Prediction / rule |
|---|---|---|
| Coverage gate | number of HIF2dep genes passing the filter | ≥ 30, otherwise the test is uninterpretable |
| Knockout check (descriptive) | HIF1A in Δ1, EPAS1 in Δ2 | CRISPR may not lower mRNA, so this is reported only |
| Positive control | hyp_all Δ1 | < 0 (the Hallmark set is mainly HIF-1-driven) |
| **Q1 primary** | D = mean over donors of (Δ2 − Δ1) for HIF2dep | **D < 0**: the program falls more with HIF-2α loss than with HIF-1α loss |
| Q1 inference | D against 1,000 random sets matched on mean-expression deciles (excluding HYP, HIF2dep and FLT1, seed 0) | primary: HIF2dep below the random 5th percentile |
| Q1 also reported | sign count over 4 donors; exact sign-flip p | minimum p with 4 donors is 0.0625 |
| Q1 comparator (descriptive) | the same D for hyp_all | expected ≥ 0 |
| FLT1 (descriptive) | FLT1 Δ1 and Δ2, if FLT1 is expressed | — |

## Reading rule

- **Gate passes, D < 0, and HIF2dep is below the random 5th percentile:** the expressed part of the
  program is HIF-2-preferring in human T cells.
- **D ≥ 0 or not beyond random:** not shown, and reported as such.
- **In either case**, trophoblast HIF-2 specificity remains untested.
