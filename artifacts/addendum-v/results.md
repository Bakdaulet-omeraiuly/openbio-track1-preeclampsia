# Addendum V results (2026-11-15): HIF-2, not HIF-1, is required for FLT1 and the program in iPSC-derived syncytium

The plan (`f9939dc`) and amendments 1–2 (`e2d9103`, `54b72c8`) were committed before any count value was read.
The amendments fixed only the column mapping, the file separator and the ID-to-symbol mapping. The script was
run **once**. Outputs: `V_results.json`; input checksums: `inputs.sha256`.

## Gates (both studies passed)

| | EPAS1 study (GSE288288) | HIF1A study (GSE302997) |
|---|---|---|
| Wild type: 3% vs 20% O2 raises Hallmark hypoxia | +0.58 SD, p = 0.021 | +0.84 SD, p = 0.0002 |
| Target mRNA, null vs wild type at 3% | EPAS1 −3.15 log2, p = 9e-6 | HIF1A −2.11 log2, p = 3e-5 |
| Wild type: hypoxia effect on FLT1 (descriptive) | +3.60 log2 (about 12×) | +3.64 log2 (about 12×) |

## Primary: null vs its own wild type at 3% O2 (Welch, one-sided "lower")

| Outcome | **EPAS1-null (HIF-2α)**, 9 vs 3 | **HIF1A-null (HIF-1α)**, 3 vs 3 | Isoform reading (fixed rule) |
|---|---|---|---|
| **FLT1** | **−5.46 log2, p = 4e-5** | +0.87 (not lower; p = 0.999) | **HIF-2-preferring** |
| **HIF2dep program** | **−1.28 SD, p = 3e-5** | +0.40 (not lower) | **HIF-2-preferring** |
| Hallmark hypoxia (isoform positive control) | +0.34 (not lower) | **−0.39 SD, p = 0.009** | HIF-1-dominated, as expected |

- **Interaction** (null blunts the hypoxic change): FLT1 −5.32 (EPAS1) vs +0.57 (HIF1A); program −0.83 vs +0.30.
- **Sensitivity (a)**, the EPAS1 full-coding deletion alone: FLT1 −5.62 (p = 3e-4); program −1.22 (p = 3e-4);
  Hallmark hypoxia unchanged.

## Post hoc differentiation check (`v2_differentiation_posthoc.py`, committed `3ab0810` before running)

EPAS1-null cells are **more** syncytium-like, not less.

- The STB marker score is +1.65 SD vs wild type. CGA +1.70, LGALS13 +1.86, SDC1 +1.13, CYP19A1 +0.81 and
  ERVFRD-1 +0.75 log2.
- Yet FLT1 falls. Adjusted for the STB score, the EPAS1-null effect on FLT1 is **−6.22 log2** (SE 1.54,
  p = 0.003).
- So the FLT1 loss is **not** a failure to differentiate. This also fits the reported repression of
  syncytialisation by HIF/ARNT (Kanai lab, PMID 37028593).
- HIF1A-null: STB score −1.10, FLT1 +0.68 after adjustment.

## Reading (by the rules written before the run)

- **Link 4b: HIF-2-preferring for both FLT1 and the program.** In trophoblast-lineage syncytium at 3% O2, loss of
  HIF-2α (EPAS1):
  - abolishes the hypoxic induction of FLT1 and lowers it further;
  - lowers the HIF2dep program;
  - does not touch canonical HIF-1 targets.

  Loss of HIF-1α does the reverse: it lowers canonical hypoxia genes but neither FLT1 nor the program.
- **This is a double dissociation in a trophoblast model, with genetic nulls.** It resolves the conflict in link 4
  (O1 ARNT shRNA; Q T cells) in favour of HIF-2, consistent with Sasagawa et al. 2018 and with PT2385 (L1).
- With Addendum T, which found no ARNT peak near FLT1 in BeWo, the most parsimonious reading is that **HIF-2 is
  required for FLT1 in trophoblast, possibly indirectly.** The two datasets come from different models, so this
  is not settled.

## Prior report

- Sasagawa et al. 2018 (Sci Rep; PMID 30478339) showed by siRNA in placental trophoblasts that HIF-2α, not
  HIF-1α, mediates hypoxic FLT1 induction. For FLT1, V is therefore an **independent genetic confirmation**.
- A 2021 study (Mol Hum Reprod; PMID 34665260) found HIF-1β (ARNT) essential for hypoxic FLT1 induction in
  trophoblasts. Addendum T found no ARNT peak near FLT1 in BeWo, so "possibly indirect" stays a hypothesis.
- The JAX GEO series (no linked paper) do not mention FLT1.
- **Our additions:** the HIF2dep-program result and the double dissociation against canonical HIF-1 targets.

## Descriptive, not pre-specified

The data show a compensation pattern:

- HIF1A-null raised FLT1 (+0.87) and the program (+0.40);
- EPAS1-null raised Hallmark hypoxia slightly (+0.34).

## Limits

- **Unit of replication.** All nulls and wild types come from **one parental iPSC line (KOLF2.2J, one male
  donor)**. Replicates are clones and differentiations, not donors.
- **Wild-type FLT1 induction.** Hypoxic induction of FLT1 was p = 0.07 in the EPAS1 study and p = 0.013 in the
  HIF1A study. "About 12×" is the point estimate in both.

- iPSC-derived "primitive syncytium" (day 6) is an early, embryonic-like model, not term villous STB.
- The two isoforms come from two studies or batches. Each is compared with its own wild type, and the isoforms
  are compared through effect sizes.
- HIF1A has n = 3 vs 3. Its FLT1 estimate points the other way (+0.87), so low power does not explain the
  null.
- For EPAS1, three editing strategies agree; HIF1A has only one strategy (termination codon).
- In culture this is a causal statement about gene requirement (HIF-2α is *required for* FLT1 here), not about direct binding.
