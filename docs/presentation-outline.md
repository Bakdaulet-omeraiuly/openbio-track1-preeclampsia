# Presentation Script (5:00 max)

Two speakers: **A** and **B**. Every member must appear or speak. Times are cumulative. Speak slowly; this script is ~784 spoken words, about 5:03 at 155 words per minute — trim one Discussion sentence if the timer runs over — rehearse with a timer at a calm pace.

---

## Introduction (0:00 to 0:45) — Speaker A

**Slide 1: title + headline figure**

> Preeclampsia is a pregnancy disorder where the mother's blood pressure becomes dangerously high. The only cure is delivering the placenta, so the early form — before 34 weeks — often means a very premature baby.
>
> Doctors think it starts in the placenta, possibly when it runs short of oxygen. A 2023 study of the placentas we analysed found that some looked like they had switched how they burn sugar.
>
> We asked: does that signal hold up against the right controls, is it specific to preeclampsia, is it really about sugar or about oxygen, and which cells carry it?
>
> Data: RNA sequencing of 30 placentas (GSE203507), an independent single-cell study of 26 pregnancies, and two more patient groups of 79 and 157 placentas.

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

**Slide 7: replication across cohorts (fig9 / headline)** — Speaker B

> Thirty placentas could be a fluke. So each time, we wrote the predictions down first, then opened a group we had never touched. We did that with 79 patients, then 157, then six more datasets — ten in all. The oxygen alarm was again higher in early preeclampsia than in births at the same stage, in five of our six matched comparisons. But one prediction failed every time: preeclampsia never clearly differed from growth restriction alone, even with 42 growth-restricted placentas. We report that as not replicated.

**Slide 8: FLT1 link (fig8)** — Speaker A

> The alarm rises with FLT1, the gene for sFlt-1 — the protein measured in the clinical preeclampsia blood test. FLT1 is not on the oxygen list. This is our most robust result: it held in all ten datasets, and we had registered the prediction in nine of them.

**Slide 9: the HIF-2 program (fig16)** — Speaker A

> Then we went one step toward mechanism. Textbook oxygen gene lists aren't built from placenta. So we took a published experiment where a drug blocks HIF-2 in placental cells, kept the hundred genes the drug switched off most, and asked whether those genes track FLT1 in patients. They did — in four separate groups of placentas we had never opened, beyond the general oxygen score. And in a different experiment, where the HIF machinery was disabled genetically instead of by drug, the same genes fell in all five donors.

## Discussion (3:45 to 5:00) — Speaker B

**Slide 10: hypothesis**

> Our hypothesis: in a subset of early preeclampsia, trophoblast cells sense low oxygen and switch on a program that includes FLT1, feeding the excess sFlt-1 that damages the mother's blood vessels. Growth restriction uses the same machinery.
>
> One worry about tissue data is that it only reflects which cells are present. So we went down to single cells, in twenty more placentas. Inside individual syncytiotrophoblast cells of the same donor, the alarm still tracked FLT1 — in nineteen of twenty donors, and beyond every random gene list. We call that provisional, because a check on how cleanly those cells were identified is still running.
>
> **We want to be clear about what we did not show.** Preeclampsia never separated from growth restriction. The alarm did not track the mother's blood-flow measurements or how ill she was. And when we tested whether the program prefers HIF-2 over HIF-1, it preferred neither. So we have a strongly replicated association and one causal step in cultured cells — not a causal chain in patients.
>
> **The decisive next test** is a wet-lab one: switch off HIF-2 genetically in placental cells from term pregnancies, at low oxygen, and measure sFlt-1 released into the medium. If the oxygen score stops tracking FLT1 in new data, we're wrong.

## Final Slide (both on camera)

Team names · Data: GSE203507 (Ackerman et al. 2023), GSE114691 (Awamleh et al. 2019), GSE75010 (Leavey et al. 2016), Admati/Skarbianskis et al. 2023, plus GSE35574, GSE24129, GSE147776, GSE44711, GSE74341, GSE66273, GSE186257 and the perturbation sets GSE216799 and GSE80995 · Repository: github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia
