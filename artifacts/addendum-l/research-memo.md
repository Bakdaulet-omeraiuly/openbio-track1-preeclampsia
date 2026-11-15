# Addendum L research memo: literature gap, candidate hypotheses and data (2026-11-15)

This covers steps 1–3 of the follow-up brief. Literature was retrieved from PubMed and data from NCBI GEO.
The repository already contains the FLT1 module, the Ackerman Blue/Black comparison, the normal
CTB→EVT reference and the isoform analysis. They are not proposed again.

## 1. What is already known

- **HIF-2α, not HIF-1α, drives hypoxic FLT1/sFlt-1 in trophoblasts.** siRNA or an inhibitor against
  HIF-2α lowers hypoxia-induced FLT1, and HIF-1β is its required partner.
  - Sasagawa et al. 2018, [10.1038/s41598-018-35745-1](https://doi.org/10.1038/s41598-018-35745-1)
  - 2021 follow-up, [10.1093/molehr/gaab065](https://doi.org/10.1093/molehr/gaab065)
- **HIF-2α and FGR.** HIF-2α-positive syncytial nuclei are more frequent in pathological placentas, and
  HIF-2 impairs syncytialisation, which is linked to FGR
  ([10.1096/fj.202001681R](https://doi.org/10.1096/fj.202001681R)). The HIF-2 inhibitor PT2385 reduces
  sFLT-1 and placental dysfunction in vitro and in the rat RUPP model
  ([10.1161/HYPERTENSIONAHA.122.20739](https://doi.org/10.1161/HYPERTENSIONAHA.122.20739)).
- **HIF-1α and HIF-2α protein are raised in both severely preterm PE and idiopathic IUGR placentas**
  compared with spontaneous preterm controls
  ([10.3390/ijms23052772](https://doi.org/10.3390/ijms23052772)). Placental oxygen-sensing changes are
  therefore shared by PE and FGR, consistent with our PE ≈ FGR results.
- **Placental sFLT-1 e15a (= sFlt1-14) protein in maternal plasma is raised in early PE and in term FGR**
  ([10.1080/14767058.2016.1182975](https://doi.org/10.1080/14767058.2016.1182975)).

**Gap.** In our bulk cohorts we do not know how much of the placental FLT1 variation is explained by a
*trophoblast HIF-2-dependent* program, as opposed to the generic Hallmark hypoxia list, which was built
from many cell types. We also do not know whether the FLT1-linked "Blue-type" genes (our F1 module)
depend on HIF-2. And the PE-versus-FGR question has never been tested in this project with more than
18 FGR-only placentas.

## 2. Candidate hypotheses (at most three)

| # | Hypothesis | Type | Supporting observation | Weakening observation |
|---|---|---|---|---|
| H1 | A trophoblast HIF-2-dependent program, defined from a perturbation experiment, tracks placental FLT1 beyond the generic hypoxia score | Mechanism (association with a perturbation-defined program) | Positive partial correlation with FLT1 after removing `hyp_all`, in used cohorts and then untouched cohorts | Partial correlation near zero, or not reproducible |
| H2 | The FLT1-linked Blue-type module (F1, 50 genes) is HIF-2-dependent in trophoblast | Mechanism (perturbation) | The F1 score falls with PT2385 under hypoxia | No fall (n = 3 placentas, so descriptive only) |
| H3 | Placental hypoxia-only activity is higher in PE than in normotensive FGR | PE specificity | PE > FGR-only in pooled untouched cohorts (about 41 FGR-only) | Pooled estimate near zero with a narrow interval |

Each hypothesis has to be separated from other explanations as follows:

- **Generic hypoxia:** H1 uses partial correlation conditional on `hyp_all`.
- **Cell composition:** a limitation. The perturbation is in purified cytotrophoblasts, but the cohorts
  are bulk tissue.
- **Gestational age:** a covariate where it is recorded (GSE35574).
- **Delivery and labour:** GSE24129 is caesarean-only; GSE35574 records spontaneous vs induced labour.

## 3. Data

**Perturbation (mechanism)**

| Accession | Design | Donors | Technology | Limits |
|---|---|---:|---|---|
| GSE216799 | Primary term cytotrophoblasts, 2.5% O2 for 4 days, PT2385 10 µM vs DMSO (HIF-2 inhibitor; the PT2385 rat/in-vitro study above) | 3 placentas, paired | RNA-seq raw counts | Only 3 donors; inhibitor, not genetic knockdown |
| GSE41341 | Cultured CTB/STB at several oxygen levels | many | Illumina array (+ methylation) | Not used in this addendum (time) |
| GSE65271 | HTR8/SVneo at 1/5/20% O2 | cell line | Illumina array | Cell line; not used |

**Untouched clinical cohorts with an FGR-only group** (only metadata lines read so far)

| Accession | Groups | Platform | Covariates | Limits |
|---|---|---|---|---|
| GSE35574 (PMID 23639576) | 19 PE, 35 IUGR arrays (incl. 8 technical repeats), 40 control | Illumina GPL6102 | GA, fetal sex, labour (spontaneous/induced), batch A/B | Mostly term deliveries |
| GSE24129 (PMID 21810232) | 8 PE, 8 normotensive FGR, 8 control | Affymetrix GPL6244 | All caesarean | No GA in metadata |
| GSE147776 (PMID 32443673) | 7 PE, 7 IUGR, 6 PE+IUGR, 8 control | Agilent GPL20844 | Maternal age | No GA, no fetal sex |

Other placenta cohorts were found but not used. GSE100415 and GSE98224 overlap the GSE75010 biobank.
GSE234729 has no FGR group.
