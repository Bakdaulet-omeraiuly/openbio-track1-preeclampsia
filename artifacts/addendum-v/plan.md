# Addendum V: HIF-1 or HIF-2 in trophoblast? EPAS1-null vs HIF1A-null iPSC-derived syncytium at 3% O2 (link 4b)

Written on 2026-11-15, before any count file of these studies was downloaded. The plan, `v_isoform.py` and
`samples.tsv` are committed together, then run **once**.

## Why

Link 4b asks whether the trophoblast HIF effect on the program and on FLT1 runs through HIF-2α (EPAS1) or
HIF-1α (HIF1A). This has not been tested in trophoblast so far:

- Addendum Q used T cells.
- ARNT, used in Addenda O1 and T, is shared by both isoforms.

## Data (viewed so far: GEO sample metadata only)

The JAX iPSC trophoblast studies. Both use the same parental line (KOLF2.2J, male) and the same "primitive
syncytium" differentiation (day 6), at **3% O2 ("hypoxia")** and 20% O2 ("normoxia").

| Study | Null line | Null samples | Wild-type samples in the same study |
|---|---|---|---|
| GSE288288 (Study 1) | **EPAS1**: 3 strategies (full coding deletion, critical exon, termination codon) × 3 | 9 hypoxia + 9 normoxia | 3 hypoxia + 3 normoxia |
| GSE302997 (Study 8) | **HIF1A**: termination codon × 3 | 3 hypoxia + 3 normoxia | 3 hypoxia + 6 normoxia |

Only samples labelled hypoxia or normoxia are used (`samples.tsv`). Wild-type samples without an oxygen label
are excluded.

**Each null is compared only with the wild type of its own study.** The two studies are separate batches.

## Processing (fixed)

- Counts: the authors' `genesCounts` table for each study.
- log2 CPM + 1. Keep genes with CPM ≥ 1 in ≥ 3 samples of that study.
- Scores:
  - **`program`:** mean z of HIF2dep genes (`artifacts/addendum-l/L1_HIF2dep_genes.txt`), z-scored across that
    study's included samples.
  - **`hyp_all`:** the same for Hallmark hypoxia.
  - **FLT1:** log2 CPM.

## Gates (positive controls; the isoform test is "not testable" for a study whose gate 1 fails)

1. Wild type: hypoxia raises `hyp_all` vs normoxia (difference > 0; Welch one-sided p < 0.05).
2. Descriptive check that the null took: target mRNA (EPAS1 or HIF1A) in null vs wild type. Point mutations can
   leave mRNA in place, so a failure here is not disqualifying.

## Primary estimands (at 3% O2: null vs its own wild type)

- **V1 (HIF-2 → FLT1):** EPAS1-null lowers FLT1 (log2 difference < 0; Welch one-sided p < 0.05; 9 vs 3).
- **V2 (HIF-1 → FLT1):** HIF1A-null lowers FLT1 (same test; 3 vs 3).
- **V3 / V4:** the same for the `program` score (EPAS1 / HIF1A).
- **Positive control for isoform separation (descriptive, as in Addendum Q):** HIF1A-null should lower
  `hyp_all` at 3% more than EPAS1-null, because canonical glycolytic hypoxia genes are HIF-1-dominated.

## Isoform reading (fixed)

For each outcome (FLT1, then the program):

- **"HIF-2-preferring":** the EPAS1 effect is significant (p < 0.05), **and** the HIF1A effect is not
  significant with |estimate| < 50% of |EPAS1 estimate|.
- **"HIF-1-preferring":** the mirror case.
- **"Both":** both effects significant.
- **"Neither / unresolved":** otherwise.

## Sensitivity (descriptive)

- (a) The EPAS1 full-coding-deletion strategy alone (3 vs 3).
- (b) The KO × oxygen interaction (does the null blunt the hypoxic change?).

## Limits stated in advance

- iPSC-derived "primitive syncytium" is an early, embryonic-like model, not term villous STB.
- The two isoforms come from two studies (two batches), so they are compared through effect sizes, not directly.
- 3 vs 3 gives low power for HIF1A.
