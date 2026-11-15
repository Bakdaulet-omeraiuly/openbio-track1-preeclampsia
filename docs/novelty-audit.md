# Novelty audit: what is known, what is replication, what may be new (2026-11-15)

This is a targeted PubMed search, not a systematic review. Searches (all run on 2026-11-15):

- `PT2385 trophoblast`
- `GSE216799`
- `(HIF-2 OR HIF2A OR EPAS1) AND placenta AND (gene signature OR transcriptome OR microarray) AND (FLT1 OR sFlt-1)`
- `(EPAS1 OR HIF-2alpha) AND preeclampsia AND (knockdown OR inhibitor) AND RNA-seq`
- `(HIF-2alpha OR EPAS1 OR HIF2) AND (preeclampsia OR pre-eclampsia) AND placenta`
- `(preeclampsia) AND placenta AND gene expression AND (hypoxia signature OR hypoxia score) AND FLT1 AND (cohorts OR meta-analysis)`
- PubMed "similar articles" for the PT2385 study

## 1. Known biology. We do not claim it.

- **Low oxygen raises FLT1/sFLT1 in trophoblasts mainly through HIF-2α.**
  - Sasagawa et al. 2018, [10.1038/s41598-018-35745-1](https://doi.org/10.1038/s41598-018-35745-1).
  - Sasagawa & Shibuya 2025 showed the same in stem-cell-derived syncytiotrophoblasts, with
    belzutifan, [10.1016/j.yexcr.2025.114773](https://doi.org/10.1016/j.yexcr.2025.114773).
- **PT2385 in hypoxic primary cytotrophoblasts** (the experiment behind GSE216799) normalised
  angiogenic-factor secretion and *enhanced* syncytiotrophoblast differentiation. It also lowered sFLT-1
  in the rat RUPP model.
  - Colson et al. 2023, [10.1161/HYPERTENSIONAHA.122.20739](https://doi.org/10.1161/HYPERTENSIONAHA.122.20739).
  - Our step 1 (FLT1 lower with PT2385) is a re-analysis of their data, not new.
- **HIF-2α protein is higher in early-onset PE placentas**, e.g. Wang et al. 2025 (HIF-2α/HILPDA),
  [10.1016/j.placenta.2025.12.009](https://doi.org/10.1016/j.placenta.2025.12.009).
- **FLT1 travels with canonical PE genes** (VT Blue), separately from a HIF/glucose module.
  - Ackerman et al. 2023, [10.1161/HYPERTENSIONAHA.122.20807](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807).
  - Our Addendum F/H module re-derives this.

## 2. Replication contribution. It adds confidence, not mechanism.

- **The hypoxia–FLT1 association within PE**:
  - nine bulk cohorts, across RNA-seq and three array platforms;
  - pooled estimates 0.47–0.62 in the pre-registered tests;
  - where it weakens is reported (GSE66273 rho 0.26).
- **Hypoxia-only score higher in early PE than in GA-matched preterm controls**:
  - pre-registered replications in GSE114691 and GSE75010;
  - Addendum M is supported but fragile;
  - not clearly higher against mostly term controls (L-C2);
  - not PE-specific versus FGR.

## 3. Possibly new. Worth stating, with its limits.

**What it is.** We linked a perturbation-defined trophoblast program to clinical placental FLT1:

- the genes lowered by PT2385 in hypoxic cytotrophoblasts (GSE216799);
- tested against FLT1 variation in clinical placenta cohorts;
- beyond the generic Hallmark hypoxia score (Addendum L-C4: 0.39 [0.19, 0.57] in three untouched
  cohorts, pre-registered);
- above all of 1,000 matched random gene sets.

**Prior work.** None of the searches above found a study that scored the GSE216799/PT2385-responsive
gene program in clinical placenta cohorts or related it to FLT1. Colson et al. analysed the RNA-seq in
vitro only. Multi-cohort placenta transcriptome studies exist, e.g. Kong et al. 2026,
[10.1016/j.compbiolchem.2026.109323](https://doi.org/10.1016/j.compbiolchem.2026.109323), which uses
GSE44711/GSE74341/GSE75010 among others. They use curated gene panels (EPAS1 as a single gene), not a
perturbation-defined HIF-2 program. Absence from a targeted search is not proof of novelty.

**Limits that cap the claim:**

- The program is "genes lower with PT2385 over 4 days of differentiation", not proven direct HIF-2
  targets. Colson et al. report that PT2385 shifts differentiation, so cell-state effects are built in.
- Composition is not excluded. Two post hoc adjustments left the pooled estimate unchanged, but one
  small cohort turned negative under full adjustment.
- It is an association; it is not PE-specific; and it does not separate PE from FGR.

## How to describe the work

**Say:**

> "A multi-cohort study that links an experimentally defined trophoblast program to placental FLT1
> variation, with pre-registered replication and its weak points reported."

**Do not say:**

- "We discovered a new mechanism of preeclampsia."
- "We found a PE-specific pathway."

## Update after Addendum O (same day)

- **Regulator, partly addressed.**
  - The PT2385-defined program was also lowered when ARNT expression was reduced by shRNA, in a different lab and at a different gestational stage (O1). This is agreement between two independent perturbations. It does not exclude every drug off-target effect, because ARNT partners and differentiation are not separated.
  - It is not shown to be HIF-2α-specific: ARNT also partners HIF-1α and AhR, and EPAS1 mRNA rose.
  - FLT1 was measured after ARNT knockdown, but its fall was not consistent across donors. This limits the claim; hypoxia also lowers FLT1 in those cells.
- **Cell state, not resolved.** The within-STB link was the same direction but not significant (O2: 0.26 [−0.14, 0.59]).
- **The claim is not upgraded.** Suggested wording: "A trophoblast program defined by PT2385, and also lowered in an independent ARNT-knockdown experiment, is associated with FLT1 expression in clinical placentas. Its HIF-2 specificity and its independence from cell state are not established."

- **Addendum P (fourth untouched cohort, RNA-seq).** The link replicates (0.40) and survives marker adjustment (0.63). It is not specific against random sets without adjustment (83rd percentile), and it is not clear within PE alone. This adds replication, not mechanism.

- **Addendum Q (HIF-1α vs HIF-2α knockout, T cells).** No support for HIF-2 specificity of the program: it fell with neither knockout. Either its HIF dependence is trophoblast-specific, or it is not made of direct HIF targets. The claim stays "PT2385- and ARNT-sensitive program", not "HIF-2 program".

## What would raise the novelty

Not another positive correlation, but evidence that separates the link from cell state and ties it to a
regulator:

- **Genetic test:** EPAS1 knockdown or knockout RNA-seq in trophoblasts at fixed differentiation state,
  then re-test the knockdown-defined program against FLT1 in a new, pre-registered cohort.
- **Tissue test:** spatial or single-nucleus data showing the program and FLT1 co-vary *within*
  syncytiotrophoblast nuclei, not only across tissue samples.

## Update after Addenda P–Y and the final review (2026-11-15, afternoon)

**Replication.**

- The hypoxia–FLT1 association within PE now covers **10 bulk cohorts**, 9 of them pre-registered (P added
  GSE186257, ρ 0.69). The pooled estimates are 0.47–0.69.
- It also holds within donor in pure single trophoblast cells (R: 17/19) and in intact, composition-adjusted
  Visium tissue (S: 18/19), in one cohort (GSE319236). This makes cell mix and dissociation stress unlikely as the
  whole explanation; it does not exclude them.
- GSE319236's authors already reported hypoxia signals in PE. The within-placenta test is our addition.

**The possibly-new finding (section 3) after review:**

- L-C4 holds when adjusted for diagnosis (0.35 [0.13, 0.53]) and within PE only (0.36 [0.03, 0.62]) (post hoc
  l8).
- It is attenuated to 0.20 with GA, sex, labour and batch in GSE35574.
- P gives 0.40, but with weaker support: n = 34 including 8 SGA-only, 83rd percentile of random sets, PE-only
  0.26.
- The same program falls with EPAS1 knockout in iPSC syncytium (V).
- **This pairing, a perturbation-defined trophoblast program linked to FLT1 in patients and genetically HIF-2
  dependent in a trophoblast model, remains the most plausible new contribution.** It joins evidence from
  different systems, so it is not causal evidence in patients.

**HIF-2 vs HIF-1 (V): known for FLT1.**

- Sasagawa et al. 2018 (siRNA) and **Sasagawa & Shibuya 2025** (belzutifan in stem-cell-derived
  syncytiotrophoblast) already showed that HIF-2 drives hypoxic FLT1/sFLT1.
- V confirms this with genetic nulls. The program-level result is new.
- The HIF1A comparison comes from a separate study (n = 3). "No such drop" is not proof of no HIF-1 effect.

**Negatives.**

- No adjusted association between malperfusion pathology and the generic hypoxia score (W, X; and Doppler in N).
- No clear association between placental FLT1 RNA and maternal severity (Y).
- FLT1 was higher with malperfusion lesions (pre-registered in ELGAN).
