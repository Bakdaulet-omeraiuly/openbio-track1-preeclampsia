## Inspiration

Preeclampsia is a pregnancy disorder in which the mother's blood pressure rises dangerously, often with kidney or
liver damage. The only cure is delivering the placenta. The early-onset form, which appears before 34 weeks, is
the most severe.

The textbook story is a chain of six links:

1. Poor blood supply starves the placenta of oxygen.
2. The placenta switches on a low-oxygen alarm.
3. A sensor protein called HIF runs that alarm.
4. HIF turns on the gene FLT1.
5. More FLT1 means more sFlt-1, the protein measured in the clinical sFlt-1/PlGF blood test
   ([Zeisler et al. 2016](https://doi.org/10.1056/NEJMoa1414838)).
6. sFlt-1 damages the mother's blood vessels.

Gene-activity data cannot test every link as written. They measure FLT1 RNA, not released sFlt-1 protein, and
they cannot show harm to the mother. We therefore turned the chain into **six testable questions** (Figure 6):

1. Is the alarm higher where pathology shows poor blood supply?
2. Is the alarm higher in early preeclampsia than in births at the same stage?
3. Does HIF drive the alarm genes in placental cells?
4. Is it HIF-2 (rather than HIF-1) that FLT1 needs?
5. Do the alarm and FLT1 rise together in placentas?
6. Does placental FLT1 track how sick the mother is?

A 2023 study of our main placentas described a "sugar-burning" (glycolytic) pattern in some early-preeclampsia
placentas ([Ackerman et al. 2023](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807)). We set out to test that
signal, and then each link of the chain, with public data.

## What it does

**Our most novel finding, in one paragraph.**

- Genes that fall when a drug blocks HIF-2 in placental cells track FLT1 in patients' placentas, **beyond the
  general low-oxygen alarm**.
  - We tested this in three patient groups we had never opened, with the prediction written down first.
  - It beat all 1,000 random gene lists.
- The same gene program falls when the HIF-2α gene is removed in stem-cell-derived placental cells. No such drop
  was seen when HIF-1α was removed, in a separate study.
- Connecting a cell-experiment program to FLT1 in patients was not found in our targeted literature search.
- **Limit:** the two kinds of evidence come from different systems, so this does not show that the program
  causes FLT1 in patients.

## How we built it

Our main data are gene-activity measurements (RNA sequencing) from 30 placentas
([GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507)):

- 15 early preeclampsia;
- 5 small-baby (growth restriction) without preeclampsia;
- 5 healthy early births;
- 5 healthy term births.

We then tested every main claim in **other, independent public datasets**:

- ten bulk-tissue patient groups;
- single-cell and intact-tissue ("spatial") data;
- DNA-binding data;
- stem-cell experiments in which HIF genes were switched off;
- two cohorts with placental pathology.

**Main-dataset results.**

- Gene by gene, 424 genes differed at least two-fold between preeclampsia and same-stage healthy placentas
  (adjusted p < 0.05), and 207 survived every robustness check.
- FLT1 was about 8 times higher and leptin about 60 times higher.
- In a public single-cell browser, FLT1 and leptin sat mainly in syncytiotrophoblast: the "trophoblast" cells
  that line the placenta's surface and exchange material with the mother's blood.

**Before opening each new dataset, we committed our prediction, test and pass rule to the public repository.**
Some follow-ups in the main dataset, including the sugar/oxygen split below, were planned only after we saw the
first results, and we label them as such.
We report every result, including failures. As a check against "any gene list would do", each gene-score result
is compared with 1,000 random gene lists of similar activity.

**The main limit is size.** Several groups have fewer than 10 placentas, and most data come from placentas
collected at delivery, so they show association, not cause.

## What we learned

**What held up.** The middle of the chain is strong:

- **In placental cells:** HIF-2 is needed for FLT1 and for a drug-defined HIF-2 gene program, and a HIF protein
  binds many of that program's genes.
- **In patients' placentas:** the general low-oxygen alarm rises and falls with FLT1 in every group we examined,
  and within single cells and intact tissue of the same placentas.
  - The HIF-2 program itself tracked FLT1 beyond the general alarm in three untouched patient groups, with weaker
    support in a fourth.
  - Its within-placenta checks were descriptive only.
  - These are two related but separate findings.

**What did not.** We could not show the two ends of the chain in our data:

- an association between pathology signs of poor blood supply and the alarm;
- that placental FLT1 RNA reflects how sick the mother is.

We also found nothing specific to preeclampsia rather than growth restriction. These are honest negatives, not
failures to look. The two ends rest on earlier animal and clinical studies
([Maynard et al. 2003](https://doi.org/10.1172/JCI17189)).

**What is new.** The low-oxygen → HIF-2 → FLT1 biology was known. Our contributions are:

- a planned, ten-group replication of the alarm–FLT1 link, with its weak points reported;
- showing that the alarm–FLT1 link holds within single cells and intact tissue of the same placentas, which
  makes cell mix and handling stress unlikely as the whole explanation;
- linking a drug-defined HIF-2 program to FLT1 in patients, with genetic evidence that HIF-2 (not HIF-1) is
  needed for both;
- clear negatives for the first and last links.

**Our hypothesis.** In a subset of early preeclampsia, and in growth restriction, placental cells run a
HIF-2-dependent program that keeps FLT1 high. We have not shown what turns it on in patients.

**Limitations.**

- Small groups in several datasets.
- Gene activity measured at delivery.
- One stem-cell line for the HIF-2 versus HIF-1 test.
- Single-cell and intact-tissue results from one study.
- Association, not cause, in all patient data.

**Next steps.**

1. Remove HIF-2α in cells from term placentas at low oxygen and measure released sFlt-1 protein.
2. Stain HIF-2α and FLT1 in placentas with and without poor-blood-supply lesions, with the mother's sFlt-1/PlGF
   ratio.
3. Find a cohort with both placental tissue and maternal blood sFlt-1 to test the last link directly.

**We are wrong if** the alarm stops tracking FLT1 in new placentas, or if removing HIF-2 in term placental cells
leaves FLT1 unchanged.

The full analysis log, with every planned and unplanned test, is in
[`supplementary-results.md`](https://github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia/blob/main/docs/supplementary-results.md).

## Links

- Repository: https://github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia
- Article: docs/article.md · Hypothesis card: docs/hypothesis-card.md · Full log: docs/supplementary-results.md
