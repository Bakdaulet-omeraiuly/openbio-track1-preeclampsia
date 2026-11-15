# Addendum T results (2026-11-15)

The plan and script were committed and pushed (`9a1c812`) before the peak file or the annotation was
downloaded. The script was run **once**. Outputs: `T_results.json`; input checksums: `inputs.sha256`.

## Setup

- 1,812 reproducible ARNT peaks (ChIP-R, 3 replicates) in BeWo cells at 1% O2.
- Universe: 12,819 genes expressed in primary trophoblast (GSE216799, baseMean ≥ 10).
- Base rate: 6.8% of these genes are bound within ±10 kb of a TSS.

## Results (±10 kb, primary)

| Test | Genes | Fraction bound | Random sets (expression-matched): median / 95th pct | Empirical p | Verdict |
|---|---|---|---|---|---|
| **T0** Hallmark hypoxia (positive control) | 174 | **34.5%** | 7.5% / 10.9% | 0.001 | **passed** (4.6× the null) |
| **T1** HIF2dep program | 82 | **24.4%** | 7.3% / 13.4% | 0.001 | **supported** (3.3× the null) |
| T2 FLT1 (descriptive) | 1 | **not bound** | — | — | nearest peak **123 kb** from any FLT1 TSS |

**Sensitivity.** T1 holds at both other windows:

- ±5 kb: 14.6% vs 6.1%, p = 0.009;
- ±50 kb: 42.7% vs 19.5%, p = 0.001.

FLT1 is unbound at ±50 kb.

**Bound HIF2dep genes** (20 of 82; PT2385 log2FC in parentheses): SPNS2 (−3.6), CCDC183 (−3.4), PNCK (−2.1),
CRISPLD2 (−2.0), **FSTL3 (−2.0)**, SYT8 (−1.9), TNNI2 (−1.9), TMEM141 (−1.7), TRIM14 (−1.7), ERVH48-1 (−1.6),
TRAM2 (−1.4), UBASH3B (−1.4), TNFSF10 (−1.4), DSP (−1.2), COL17A1 (−1.2), PXN (−1.1), CLIC3 (−1.1), DPP7 (−1.1),
FHL2 (−0.9), MTUS1 (−0.8).

- FSTL3 is a canonical PE placental gene.
- ERVH48-1 is a syncytin-family gene.

## Reading (by the rules written before the run)

- **T1 is supported.** In trophoblast cells at low oxygen, the HIF complex binds near the genes that HIF-2
  inhibition lowers 3.3 times more often than expected for genes of similar expression.
  - The program therefore **includes direct HIF targets**, not only downstream effects of differentiation.
  - Together with PT2385 (L1) and ARNT shRNA (O1), link 3 now has binding, pharmacological and genetic support.
- **FLT1 is not bound** in this dataset (nearest peak 123 kb).
  - FLT1 has a hypoxia-response element described in other cells (Gerber et al. 1997).
  - **In trophoblast, FLT1's hypoxia link may therefore be indirect**, via the HIF-driven program or STB
    differentiation, rather than direct HIF binding at FLT1.
  - This is a new, specific constraint on link 4.
  - Caveats: the peak set is stringent (1,812 peaks), so missing a peak is weaker evidence than finding one.
    BeWo is a choriocarcinoma line, and very distal enhancers (> 50 kb) were not considered.
- **4b (HIF-2 vs HIF-1)** remains untested. ARNT is shared by both isoforms.

## Limits

- One cell line (BeWo), one condition (1% O2, 24 h), and one antibody (ARNT).
- Binding does not show that it regulates expression.
- TSS windows are a simple assignment rule; distal enhancers are assigned to their nearest window or not at all.
