# Presentation Script (5:00 max)

Two speakers: **A** and **B**. Every member must appear or speak. Times are cumulative. Speak slowly; this script is ~650 words, about 4:40 at a calm pace.

---

## Introduction (0:00 to 0:45) — Speaker A

**Slide 1: title + headline figure**

> Preeclampsia is a pregnancy disorder where the mother's blood pressure becomes dangerously high. The only cure is delivering the placenta, so the early form — before 34 weeks — often means a very premature baby.
>
> Doctors think it starts in the placenta, possibly when it runs short of oxygen. A 2023 study of the placentas we analysed found that some looked like they had switched how they burn sugar.
>
> We asked: does that signal hold up against the right controls, is it specific to preeclampsia, is it really about sugar or about oxygen, and which cells carry it?
>
> Data: RNA sequencing of 30 placentas (GSE203507), plus an independent single-cell study of 26 pregnancies.

## Methodology (0:45 to 1:45) — Speaker B

**Slide 2: groups table (15 PE, 5 FGR-only, 5 preterm controls, 5 term controls)**

> Preeclampsia forces early delivery, so comparing with full-term placentas confuses disease with timing. Our main comparison was 15 early-preeclampsia placentas against 5 healthy placentas delivered at 30 to 33 weeks, adjusting for the baby's sex — and in a sensitivity analysis, for the exact delivery week.
>
> We used DESeq2 with cut-offs fixed in advance: adjusted p below 0.05 and at least a two-fold change.
>
> **Slide 3: "Pre-registered" screenshot of the git history**
>
> We committed our analysis plan before running any test, and every later analysis was registered before we ran it. Each pregnancy counts once. We checked robustness by dropping each control in turn, by comparing against random gene sets, and in a second, independent dataset.

## Results (1:45 to 3:45)

**Slide 4: PCA (fig1)** — Speaker A

> On a map of overall gene activity, preeclampsia placentas spread out to the right, while early and full-term healthy placentas overlap — so disease, not just timing, drives the big pattern.

**Slide 5: volcano (fig2)** — Speaker A

> 424 genes differ; 207 survive every robustness check, including beating 1,000 random gene lists later on. FLT1 is eight times higher, leptin about sixty times.

**Slide 6: group contrasts (fig6)** — Speaker B

> We split the published "low oxygen" and "sugar burning" gene lists into shared genes, oxygen-only genes and sugar-only genes. Oxygen genes are up in preeclampsia — in this dataset even compared with placentas of small babies *without* preeclampsia. Sugar-only genes are not up at all. So the "glycolytic" placentas are really placentas with a low-oxygen alarm.

**Slide 7: replication (fig9 / headline)** — Speaker B

> Thirty placentas could be a fluke, so we wrote down six predictions and tested them in 79 different patients. The oxygen alarm was again higher in early preeclampsia, again beat random gene lists, and sugar-only genes actually fell. One prediction failed: preeclampsia versus growth restriction alone pointed the same way but wasn't clear — growth-restricted placentas show a weaker version too. We report that as not replicated.

**Slide 8: FLT1 link (fig8)** — Speaker A

> The alarm rises with FLT1, the gene for sFlt-1 — the protein measured in the clinical preeclampsia blood test. FLT1 is not on the oxygen list, and the link repeats in all three datasets — including inside syncytiotrophoblast cells of preeclampsia donors alone.

**Slide 9: single-cell + browser (fig5 + browser table)** — Speaker A

> Inside trophoblast cells from the independent study, the oxygen score is higher in early preeclampsia — suggestive, with only 3 controls. In the gene browser, FLT1 and leptin sit mainly in syncytiotrophoblast; other alarm genes like NDRG1 and P4HA1 are active across all trophoblast.

## Discussion (3:45 to 5:00) — Speaker B

**Slide 10: hypothesis**

> Our hypothesis: in a subset of early preeclampsia, trophoblast cells sense low oxygen and switch on a HIF-driven program that includes FLT1, feeding the excess sFlt-1 that damages the mother's blood vessels. Growth restriction may share a weaker version of this stress.
>
> **Limitations:** only 5 preterm controls in our main data and 3 single-cell controls; no per-patient details in the second study; bulk tissue mixes cell types; labour and delivery type were not recorded; and the maternal-side tissue shows a similar signal, so we cannot call it purely placental. This is association, not cause.
>
> **To test it:** a third cohort with recorded delivery mode, to settle preeclampsia versus growth restriction; stain these placentas for HIF and FLT1 and compare with the mother's sFlt-1/PlGF ratio; and grow trophoblast at low oxygen with and without HIF. If the oxygen score stops tracking FLT1 in new data, we're wrong.

## Final Slide (both on camera)

Team names · Data: GSE203507 (Ackerman et al. 2023), GSE114691 (Awamleh et al. 2019), Admati/Skarbianskis et al. 2023 · Repository: github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia
