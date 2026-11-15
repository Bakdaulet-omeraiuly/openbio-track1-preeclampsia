# Addendum Q results: HIF-1α vs HIF-2α knockout and the PT2385-defined program (GSE324038, human T cells)

The plan and script were committed before any expression value was read (`3c98efb`). The script ran once.
Full output is in `Q_results.json`.

| Item | Result |
|---|---|
| Coverage gate | 59 of 100 HIF2dep genes expressed: **passed** |
| Knockout check | HIF1A mRNA −1.19 with HIF1A-KO; EPAS1 mRNA −0.95 with EPAS1-KO (each specific) |
| Positive control | hyp_all Δ: −0.73 with HIF1A-KO vs −0.17 with EPAS1-KO, so the Hallmark set is HIF-1-dominated: **passed** |
| **Q1 primary** | HIF2dep Δ: **+0.21 with HIF1A-KO, +0.12 with EPAS1-KO**; D = −0.09 (3/4 donors negative, sign-flip p = 0.19) |
| Q1 calibration | random median −0.12, 5th percentile −0.32; HIF2dep at the **57th percentile** |
| **Verdict** | **not shown** |
| FLT1 (descriptive) | +0.36 with HIF1A-KO, **−0.66 with EPAS1-KO** |

## Reading

- **No HIF isoform preference in T cells.** In hypoxic human T cells, the program genes that T cells express
  do not fall with either knockout, and they show no HIF-2 preference beyond random genes. So the experiment
  gives no support for HIF-2 specificity of the program.
- **Two explanations remain possible:**
  - the program's HIF dependence is trophoblast-specific (context-dependent targets); or
  - the program is not made of direct HIF targets, and reflects the differentiation shift that PT2385 and
    ARNT loss cause in trophoblast.
- **The data cannot tell these apart.**
- **FLT1 itself** falls with EPAS1 knockout and rises with HIF1A knockout here. This fits the known HIF-2
  preference of FLT1 (Sasagawa et al. 2018), but it is descriptive and comes from T cells.
- **The question stays open.** HIF-2 specificity in trophoblast needs a trophoblast experiment with
  separate HIF1A and EPAS1 loss.
