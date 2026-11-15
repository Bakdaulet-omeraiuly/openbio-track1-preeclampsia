# Addendum U: do direct HIF targets (ARNT-bound program genes) track FLT1 in intact human placenta?

Written on 2026-11-15, before this question was run on any tissue data. The plan, `u_direct_targets.py` and the
two gene lists are committed together, then run **once**.

## Why

The chain breaks between:

- link 3: HIF drives the program in culture, now including direct binding (Addendum T); and
- link 5: hypoxia ↔ FLT1 in tissue.

Addendum T also found that ARNT does **not** bind near FLT1 in trophoblast. If FLT1 follows the HIF program
indirectly, the **direct HIF targets** of that program should co-vary with FLT1 in intact tissue beyond the
generic hypoxia score.

Expression correlation cannot prove direct or indirect causation. This test asks only whether trophoblast
direct-HIF-target genes travel with FLT1 in human placenta.

## Gene sets (fixed from ChIP alone, no tissue data used)

- `HIF2dep_ARNT_bound.txt`: the **20** HIF2dep genes with an ARNT peak within ±10 kb of a TSS (Addendum T).
- `HIF2dep_ARNT_unbound.txt`: the **62** HIF2dep genes in the same universe without a peak.

## Data and method

Identical to Addendum S: GSE319236 Visium, intact tissue, 19 donors. We use its committed functions with the
HIF2dep set replaced by each list.

- **Statistic:** within-donor partial Spearman ρ(set score, FLT1 | `hyp_all` + spot composition + technical
  covariates). The S2 statistic of Addendum S.
- **Inference:** one-sided sign-flip p across donors, and 1,000 expression-decile-matched random gene sets of
  the same size (seed 0).
- In Addendum S, the full 82-gene list gave 0.170, 19/19 descriptively. That analysis was seen before this plan,
  so U1 is not independent of it for the union of the two sets.

## Tests

**U1 (confirmatory).** The bound (direct) set tracks FLT1 beyond `hyp_all`.

- **Supported** if the sign-flip p is < 0.05 **and** the observed value is above the random-set 95th percentile.
- The 20 genes are a subset of the 82 already seen in S2. U1 asks whether the **direct** subset alone carries the
  signal, so it is confirmatory for this subset only.

**U2 (descriptive).** The same test for the unbound set. Also reported: bound minus unbound, per donor.

## Pre-written readings

- **U1 supported.** "Genes that HIF binds directly in trophoblast co-vary with FLT1 within the same intact
  placenta, beyond the generic hypoxia score."
  - This is consistent with FLT1 following a direct-HIF-target program rather than direct HIF binding at FLT1.
  - It is still association.
- **U1 not supported while U2 is.** "The FLT1 link is carried by indirect (unbound) program genes. That is
  consistent with a differentiation route."
- **Neither.** Reported as is.
