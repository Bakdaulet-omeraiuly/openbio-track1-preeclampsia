# Addendum T: does the HIF complex bind FLT1 and the HIF2dep program directly in trophoblast? (GSE209936 ARNT ChIP-seq)

Written on 2026-11-15, before the peak file or the gene annotation was downloaded. The plan and `t_chip.py` are
committed together, then run **once**.

## Why

Link 4 in the hypothesis card (HIF-2 → FLT1) is the weakest link. It has two parts:

- **4a.** Does the HIF complex bind near FLT1 and near the HIF2dep program genes in trophoblast at low oxygen?
  This would mean direct targets, not downstream effects of differentiation.
- **4b.** Is that HIF-2 rather than HIF-1?

Only 4a can be tested here. ARNT (HIF-1β) is the shared partner of HIF-1α and HIF-2α, so ARNT ChIP cannot
separate the isoforms. **4b stays untested.** No public trophoblast data separate the isoforms; see Addendum Q
and the GEO screen.

## Data

- **GSE209936** (Kanai lab; PMID 37028593). BeWo choriocarcinoma trophoblast cells cultured at 1% O2 for 24 h.
- ARNT ChIP-seq in 3 biological replicates, hg38.
- We use the authors' reproducible peak set across replicates, `GSE209936_output_optimal.bed.gz` (ChIP-R).
- Gene annotation: UCSC hg38 `ncbiRefSeq` transcripts. TSS = transcript start, strand-aware, using NM_
  transcripts only.
- **Universe:** genes expressed in primary trophoblast, i.e. baseMean ≥ 10 in
  `artifacts/addendum-l/L1_GSE216799_PT2385_vs_DMSO.csv`, which also have a RefSeq NM_ TSS.

## Definition (fixed)

A gene is **bound** if any reproducible ARNT peak overlaps ±10 kb around any of its NM_ TSSs.

## Tests

**T0, positive control (gate).** Hallmark hypoxia genes (v2024.1, in the universe) are bound more often than
1,000 random sets.

- The random sets have the same size and are drawn within baseMean deciles matched to the hypoxia set.
- Seed 0.
- If the empirical p ≥ 0.05, the assay is not informative here, and T1 is reported as "not testable".

**T1 (confirmatory).** The HIF2dep program genes (`artifacts/addendum-l/L1_HIF2dep_genes.txt`, in the universe)
are bound more often than 1,000 random sets.

- The sets are matched as in T0, and drawn without hypoxia or HIF2dep genes.
- **Supported** if the empirical one-sided p < 0.05 **and** T0 passed.

**T2 (descriptive).** Is FLT1 bound (±10 kb)?

- Also reported: the distance from FLT1's TSS to the nearest peak, and the base rate of bound genes in the
  universe.
- A single gene being bound is weak evidence on its own, so T2 has no verdict.

**Sensitivity analyses:** windows of ±5 kb and ±50 kb, for T0 and T1.

## Pre-written readings

- **T1 supported.** "In trophoblast cells at low oxygen, the HIF complex binds near the genes that HIF-2
  inhibition lowers more often than expected. This suggests they include direct HIF targets, not only
  differentiation effects."
  - It does not show which HIF isoform.
  - It does not show that binding regulates expression.
  - BeWo is a cell line.
- **T1 not supported (T0 passed).** "The HIF2dep program is not enriched for direct HIF binding; its response to
  PT2385 may be mostly indirect (differentiation)."
- **T2.** Reported as is, with the base rate.
