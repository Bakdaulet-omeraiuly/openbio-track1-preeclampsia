# Presentation Script (5:00 max)

**Speaker:** Bakdaulet Omertayev (team 404: Team Not Found, one-person team), on camera at the start and the end.

The spoken text below is about 480 words, roughly 3:15 at 155 words per minute. That leaves time to pause on each
figure. Times are cumulative. Rehearse once with a timer.

---

## Introduction (0:00 to 0:45)

**Slide 1: title + headline figure (`figures/headline.png`, the alarm vs FLT1 in 10 patient groups)**

> Hi, I'm Bakdaulet, and this is team 404: Team Not Found — a team of one.
>
> Preeclampsia is a pregnancy disorder where the mother's blood pressure becomes dangerously high. The only cure
> is delivering the placenta, so the early form often means a very premature baby.
>
> The textbook story is a chain: poor blood supply starves the placenta of oxygen, an oxygen sensor called HIF
> switches on the gene FLT1, and its protein, sFlt-1 — the one in the clinical blood test — damages the
> mother's vessels. I tested that chain, link by link, with public data.

## Methodology (0:45 to 1:45)

**Slide 2: data — 30 main placentas; then 10 patient groups, single-cell, intact tissue, ChIP, knockouts, pathology**

> My main dataset is RNA sequencing of 30 placentas. Preeclampsia forces early delivery, so I compared
> preeclampsia with healthy placentas delivered at the same stage, adjusting for delivery week and the baby's sex.
> Each pregnancy counts once.
>
> **Slide 3: git history — "plan committed before data"**
>
> Then the rule I kept all the way through: before opening any new dataset, I committed my prediction and pass
> rule to the repository, and I report every result, including failures. And every gene score is checked
> against a thousand random gene lists, so "any list would do" can't fool me.

## Results (1:45 to 3:45)

**Slide 4: sugar vs oxygen (fig6)**

> A 2023 study called some early-preeclampsia placentas "glycolytic". When I split the sugar genes from the
> oxygen genes, only the oxygen genes were up. It's a low-oxygen alarm.

**Slide 5: headline — alarm vs FLT1 in 10 groups (fig17)**

> And that alarm rises with FLT1 — in all ten patient groups I tested, nine of them with the prediction written
> down first. That's my strongest result.

**Slide 6: within one placenta (fig18)**

> Tissue mixes cell types, so I looked inside single placentas. Among single trophoblast cells, and in intact
> tissue slices never broken into cells, spots with more alarm carried more FLT1 — in seventeen and eighteen of
> nineteen placentas.

**Slide 7: HIF-2 vs HIF-1 (fig19)**

> Which sensor? In stem-cell-derived placental cells at low oxygen, removing HIF-2 cut FLT1 about forty-fold;
> removing HIF-1 did not. That confirms an earlier study, and HIF binds a quarter of the program genes directly
> — but not FLT1 itself.

**Slide 8: what did not hold (fig20)**

> Now the honest part. Poor blood supply, judged by a pathologist in two separate cohorts, did not switch on
> the general alarm — although FLT1 itself was higher. Placental FLT1 did not track how sick the mother was.
> And preeclampsia never separated from growth restriction.

## Discussion (3:45 to 5:00)

**Slide 9: hypothesis**

> So my hypothesis is narrower than the textbook. In a subset of early preeclampsia, and in growth
> restriction, placental cells run a HIF-2-dependent program that keeps FLT1 high. The middle of the chain is
> strong. The two ends — what switches it on, and whether it explains the mother's illness — are not shown in
> these data.
>
> Limits: small groups, tissue taken at delivery, one stem-cell line for the HIF test, and association, not
> cause, in patients.
>
> **Slide 10: next test**
>
> The decisive next experiment: remove HIF-2 in cells from term placentas at low oxygen and measure the sFlt-1
> they release. If FLT1 stays unchanged, or if the alarm stops tracking FLT1 in new placentas, I'm wrong.
>
> Thanks for watching.

## Final Slide (on camera)

Bakdaulet Omertayev · team 404: Team Not Found · Repository:
github.com/Bakdaulet-omeraiuly/openbio-track1-preeclampsia

Data:

- **Main data:** GSE203507 (Ackerman et al. 2023).
- **Patient groups:** GSE114691, GSE75010, GSE35574, GSE24129, GSE147776, GSE44711, GSE74341, GSE66273,
  GSE186257.
- **Single-cell and spatial:** Admati/Skarbianskis 2023, GSE319236.
- **Cell experiments:** GSE216799, GSE80995, GSE209936, GSE288288, GSE302997.
- **Pathology:** GSE303463, GSE154829.
