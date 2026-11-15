# Addendum P: fourth untouched cohort for the HIF2dep–FLT1 link, with a cell-state adjustment (written before running)

Written on 2026-11-15 at about 13:35 EST, before any expression value of GSE186257 was read. Only the
series-matrix metadata, the header row and the first three gene IDs of the supplementary table were viewed.
The plan is committed together with `p_confirm.py`.

## Data

[GSE186257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186257) (Horii et al. 2023, Mod Pathol,
[10.1016/j.modpat.2022.100035](https://doi.org/10.1016/j.modpat.2022.100035), Parast lab, UCSD):

- **Samples:** bulk placental total RNA-seq (NovaSeq) from 26 placentas with severe PE (25 SGA) and 18
  non-hypertensive controls (8 SGA, 10 not SGA).
- **Values:** the authors' DESeq2-normalised, filtered table (17,536 genes); we use log2(x + 1).
- **Covariates:** fetal sex only. GEO gives no gestational age.
- **Prior use:** the cohort has not been used anywhere in this repository. The paper reports hypoxia
  pathways in PE with maternal vascular malperfusion, and STB-high and EVT-high PE subgroups. So a hypoxia
  difference itself is not new here.

## Groups

| Group | Definition | n |
|---|---|---:|
| PE | severe preeclampsia (any SGA status) | 26 |
| SGA-only | non-hypertensive control with SGA = yes | 8 |
| Control | non-hypertensive, SGA = no | 10 |

## Scores

Scores are as in L3: gene z across all 44 placentas, then the mean z of each set.

- HIF2dep: `artifacts/addendum-l/L1_HIF2dep_genes.txt`
- hyp_all and hyp_only: Hallmark hypoxia, and hypoxia minus glycolysis
- FLT1: the gene z score
- STB marker score: CGA, the CGB family (CGB/CGB3/CGB5/CGB8, as in L5), ERVW-1 and PSG1
- EVT marker score: CEACAM1, HLA-G, ITGA1 and MMP2

## Tests

**Confirmatory:** Holm across P1 and P2; one-sided; permutation p from 10,000 shuffles of FLT1 residuals.

| Test | Set | Statistic | Prediction |
|---|---|---|---|
| **P1** (replicates L-C4) | PE + SGA-only (n = 34) | partial Spearman(HIF2dep, FLT1 \| hyp_all) | > 0 |
| **P2** (replicates C3/M2) | within PE (n = 26) | Spearman(hyp_all, FLT1) | > 0 |

**Secondary (cell state; reported, not gated):**

| Test | Statistic |
|---|---|
| P1-cs | P1 additionally adjusted for the STB and EVT marker scores |
| P1-PE | P1 within PE only |

**Descriptive:**

- hyp_only, PE vs non-SGA controls (adjusted for fetal sex; no GA available, so PE vs control is
  confounded by delivery timing);
- PE vs SGA-only;
- the 1,000-random-set null for P1, as in L4 (decile-matched, seed 0).

## Reading rule

- **P1 Holm p ≤ 0.05:** the HIF2dep–FLT1 link replicates in a fourth untouched cohort.
- **P1-cs:** if it stays > 0 with a CI excluding 0, the link is not explained by these marker scores (a
  coarse proxy). If not, cell state could explain it in this cohort.
- **Otherwise:** not replicated, and the result is reported.

## Limits

- No GA is available, and the controls are probably later deliveries.
- Small SGA-only group.
- The marker scores are coarse proxies for cell state.
- The same lab produced GSE80995 (Addendum O1). The cohorts are different, but the lab is the same.
