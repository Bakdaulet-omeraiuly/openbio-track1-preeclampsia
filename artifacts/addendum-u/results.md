# Addendum U results (2026-11-15)

The plan, script and gene lists were committed and pushed (`9d121d3`) before running. The script was run
**once**. Output: `U_results.json`.

Data: GSE319236 Visium, intact tissue, 19 donors, the Addendum S design. The statistic is the within-donor
partial ρ(set, FLT1 | `hyp_all` + spot composition + technical covariates). The number of genes from each list
that passed the spatial ≥ 5% detection filter was not written to the output file.

| Set | Mean ρ | ρ > 0 | Sign-flip p | Random sets: median / 95th pct | Verdict |
|---|---|---|---|---|---|
| **U1** ARNT-bound program genes (direct HIF targets, 20) | **0.155** | **18/19** (negative: FGS) | ≤ 1e-4 | 0.004 / 0.061; observed above all 1,000 | **Supported** |
| U2 unbound program genes (62, descriptive) | 0.125 | 19/19 | ≤ 1e-4 | −0.002 / 0.041; above all | — |
| Bound minus unbound (descriptive) | +0.030 | bound higher in 12/19 donors | — | — | — |

## Reading (by the rules written before the run)

- **U1 is supported.** Genes that HIF binds directly in trophoblast (Addendum T) co-vary with FLT1 within the
  same intact placenta, beyond the generic hypoxia score, spot composition and random gene sets.
  - This **joins the direct-binding part of link 3 to link 5 in human tissue**, as an association.
  - It is consistent with FLT1 following a direct-HIF-target program rather than HIF binding at FLT1 itself.
- **The unbound genes carry the link too** (0.125, 19/19). The bound set is only slightly stronger (+0.03,
  12/19 donors), so this test **does not separate** a direct-target route from a differentiation route.
  Both parts of the program travel with FLT1.

## Limits

- Association within tissue, not cause.
- The 20 bound genes are a subset of the 82-gene program, whose union was seen descriptively in Addendum S (S2).
- Binding comes from one cell line (BeWo), and the tissue data come from one cohort.
