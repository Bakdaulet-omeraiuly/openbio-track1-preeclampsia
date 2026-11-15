# Addendum O2: does the HIF2dep–FLT1 link (L-C4) hold within one trophoblast cell type?

Written 2026-11-15, before any HIF2dep or FLT1 value was computed in the pseudobulk data. The plan and
`o2_within_celltype.py` are committed together, then run **once**.

## Why

L-C4 found that genes lowered by HIF-2 inhibition (HIF2dep) track placental FLT1 beyond a generic hypoxia score,
in three untouched bulk cohorts. Its main open caveat is cell composition. HIF2dep contains
syncytiotrophoblast (STB) genes, so placentas with more STB could show more of both. A test **within one cell
type** removes between-class mixing.

## Prior use of these data (disclosure)

- The datathon single-cell pseudobulk (`track01_pseudobulk_{counts,metadata}.csv`) was used in Addenda A/D/F
  (notebooks 02, 07, 09: hypoxia score, the F1 module, FLT1 coupling).
- The class-mean profiles were used as the deconvolution reference in L6.
- A grep of `notebooks/*.ipynb` and `artifacts/*/*.py` for "hif2" finds HIF2dep only in Addendum L/M scripts.
  None of those score HIF2dep per pseudobulk, and the only one that touches the pseudobulk (L6) uses class means.
- **The per-donor HIF2dep–FLT1 relation in pseudobulk has not been computed before.**
- Known from earlier work, which this test is therefore not independent of: the hypoxia score tracks FLT1 in STB
  (notebook 02/07; rho 0.72, and 0.84 in early PE).
- Only metadata (donor counts per class) were viewed for this plan.

## Data

- One pseudobulk per donor and class, with `n_cells ≥ 20`:
  - **STB**: 26 donors (10 early PE, 7 late PE, 3 early controls, 6 late controls);
  - villous cytotrophoblast (VCT): 26;
  - extravillous trophoblast (EVT): 15.
- All clinical groups are included.

## Definitions (fixed)

**Normalisation, per class** (as notebook 09 F3):

- keep genes with ≥ 10 counts in ≥ 3 samples of that class;
- log2(CPM + 1);
- z-score each gene across donors.

**Scores.**

- `hif2dep`: mean z of `artifacts/addendum-l/L1_HIF2dep_genes.txt` genes;
- `hyp_all`: mean z of Hallmark hypoxia (v2024.1);
- FLT1 is excluded from both sets (FLT1 is in neither list; the script asserts this);
- the outcome is FLT1 log2(CPM + 1).

**Partial correlation.** As in L-C4:

- linear residuals of hif2dep and FLT1 on [1, hyp_all];
- then the Spearman correlation of the two residuals.

## Tests

**O2 primary (confirmatory, one test).** STB, all 26 donors: partial ρ(hif2dep, FLT1 | hyp_all) > 0.

- **Decision p:** one-sided permutation p. The FLT1 residual is permuted, 10,000 permutations, seed 42.
- **Supported** if p < 0.05 and ρ > 0.

**Pre-specified secondaries (descriptive, no verdict).**

- **S1. STB with clinical covariates.** Residualise on [1, hyp_all, gestational age, PE status, onset group].
  Across all groups, PE status itself could drive both variables; S1 checks that.
- **S2. Other cell types.** The same partial ρ in VCT (n = 26) and EVT (n = 15).
- **S3. Random-set null in STB.** 1,000 random gene sets of the same size as the measured HIF2dep set, drawn
  within mean-expression deciles matched to HIF2dep's decile counts.
  - Excluded from the draws: Hallmark hypoxia genes, HIF2dep genes and FLT1.
  - Seed 0.
  - Reported as the percentile of the observed ρ.
- **S4. STB within PE donors only** (n = 17).

## Pre-written readings

- **Supported.** "Within syncytiotrophoblast alone, HIF2dep genes track FLT1 beyond the hypoxia score, so the L-C4
  link is not only a matter of how much STB a sample contains."
  - It does not exclude STB sub-states, for example more mature or more fusion-active STB.
  - With 26 donors (3 early controls), the CI is wide.
- **Positive but not supported.** "Same direction in STB, not significant (n = 26)." This is not evidence against
  L-C4.
- **At or below zero.** "Within STB we see no HIF2dep–FLT1 link beyond hypoxia. The bulk L-C4 link may then
  reflect composition or STB sub-state." This weakens L-C4 and must be stated in the article.
- **S3.** If the observed ρ is below the 95th percentile of random sets, we say the STB link is not specific
  to HIF2dep, whatever the primary p is.
