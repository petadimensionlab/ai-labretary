# Appendix A — A Record of What We Tried and Why It Failed

- Prepared by Labretary, 2026-10-03 (draft; inclusion in the main report is the PI's decision)
- Sources: `hypotheses/tracker.md`; `logs/learnings-2026-09-25` through `2026-10-03.jsonl`; `.claude/memory/{pi-patterns,guardrails,labretary-level}.md`; `tasks/long-jobs.md`; `experiments/EXP-20261002-00{1,2}/results/criteria_00{6,7}.json`
- **This appendix adds no new conclusions. Every verdict below has been confirmed by the PI** (the last of them, gate G3 on HYP-20261001-01, was confirmed on 2026-10-03).

---

The main body states only the acceptance criteria we finally adopted and the results obtained under them. This appendix records what we tried on the way there, what did not work, and why.

**The unit of counting in this appendix is the verdict per hypothesis.**

> **Of the 10 hypotheses that received a verdict, 5 were inconclusive.**
> **We designed 17 experiments and ran perturbation computations in 15 of them.**
> **These were not consecutive.** Four refutations and one partial support fall in between.

The two experiments that ran no perturbation at all are EXP-20260927-002 (a feasibility check) and EXP-20261002-001 (a re-evaluation).

The breakdown, from `hypotheses/tracker.md`, as of 2026-10-03:

| Verdict | Count | Hypotheses |
|---|---|---|
| Inconclusive | **5** | HYP-20260925-01 / HYP-20260925-02 / HYP-20260926-04 / HYP-20260927-01 / HYP-20260930-01 |
| Refuted | **4** | HYP-20260926-01 / HYP-20260926-02 / HYP-20260928-01 / HYP-20261001-01 |
| Partially supported | 1 | HYP-20260926-03 |
| (For reference) verdict withheld | 1 | HYP-20260929-01 |

**A refutation is an answer; it does not count as "could not decide."** In particular, HYP-20260928-01 (EXP-20260928-002, confirmed as refuted by the PI on 2026-09-29) produced an answer, and the run of inconclusive verdicts is broken there. The last of them, HYP-20261001-01, was likewise confirmed as refuted on 2026-10-03, on the strength of EXP-20261002-002 with its control panel enlarged to 1,000 genes.

**The purpose of this appendix is to let the reader follow what those five inconclusive verdicts actually were.**

A note on units. Unless stated otherwise, the size of a cell-state movement is expressed in **standard deviations**, where one standard deviation is the spread of a large panel of control genes that should have no effect. "0.41 standard deviations" therefore means a shift toward the goal state equal to 0.41 times the scatter of the control genes.

---

## A.1　The sequence at a glance

| Experiment ID | Date | Question | Verdict | Why it stopped | What we learned |
|---|---|---|---|---|---|
| EXP-20260926-001 | 2026-09-26 | In failed-repair proximal tubule cells from chronic kidney disease kidneys rich in tertiary lymphoid tissue, does deleting TP53, CDKN1A or LTBR move cells toward the healthy proximal tubule state? | **Inconclusive** | Deleting any gene moved cells away from all three state centroids by a similar amount (all 90 control genes × all seeds). The positive controls VCAM1 and HAVCR1 could not be distinguished from the housekeeping genes GAPDH and PPIA, so passing meant nothing specific | Deletion carries a uniform movement that does not depend on which gene is deleted. A metric that subtracts it is required |
| EXP-20260926-002 | 2026-09-26 | With the uniform movement subtracted, do the positive controls separate from both the control genes and the housekeeping genes? | **Refuted** | None of 18 candidate metrics satisfied the criteria across all three seeds. Even the best one gave VCAM1 0.67 and HAVCR1 0.78 against a threshold of 0.90 | Subtraction alone is not enough. The strictness of a criterion must itself be quantified (chance pass rate 5.5% per gene) |
| EXP-20260926-003 | 2026-09-26 | Does a model fine-tuned to discriminate failed-repair from healthy cells separate the positive controls? | **Refuted** (secondary predictions supported) | The uniform movement largely vanished (fraction of control genes moving away from the goal: 1.00 → 0.10–0.15), but HAVCR1 failed at 0.45 / 0.75 / 0.75. The one gene that passed, VCAM1, is explained by removal of a marker the classifier itself relies on | Fine-tuning fixes the uniform movement but does not fix the positive-control problem |
| EXP-20260926-004 | 2026-09-26 | Does VCAM1's effect exceed that of 10 expression-matched non-driver genes, and does the verdict reproduce donor by donor? | **Partially supported** | It did exceed them in all three seeds. However, split-half reproducibility of the gene ranking was 0.53, and doubling the donors from 14 to 29 did not reach the 0.7 threshold | The main criterion is not specific to causal drivers: 3 of the 10 expression-matched genes met it as well |
| EXP-20260926-005 | 2026-09-26–27 | If literature-established causal genes are used as positive controls, do they separate? | **Inconclusive** | Deletion: 0 of 5 passed, with the assay demonstrably working — a valid negative result. Overexpression: 3 of 6 passed, but the housekeeping gene GAPDH landed in the top 5% in all three seeds, so a control that should be neutral was not, and the pre-registered rule forced "inconclusive" | Under overexpression, the effect correlates strongly with how low in the input ranking the gene originally sat (0.32–0.71 at gene level). Controls matched on that confounder are required |
| EXP-20260927-001 | 2026-09-27 | With controls matched on original rank and input length, is the confounder removed, and do at least 2 of 6 positive controls separate? | **Inconclusive** | Confounder removal itself succeeded (correlation with original rank 0.32–0.65 → 0.19–0.24). But 0 of 6 passed: two positive controls (ESRRA, PPARGC1A) pointed the wrong way, and the control scatter (SD 1.03–1.36 × 10⁻³) was the same order as the positive-control effects (+0.43 to +1.34 × 10⁻³) | By the rule this is a refutation, but a refutation cannot be written when the positive controls are not working. Statistical power must be estimated in advance from the real data's scatter |
| EXP-20260927-002 | 2026-09-27 | (Feasibility check, no perturbation) Can the goal and starting fibroblast states be told apart in the embedding at all? | Passed (not a hypothesis verdict) | — | Sequencing depth alone (how many genes are detected per cell) separated the groups (AUC 0.676). Depth must be matched before measuring (after matching, AUC 0.704; seeds 0.721 / 0.690 / 0.702) |
| EXP-20260927-003 | 2026-09-27 | Does overexpressing LTBR, RELB or NFKB2 in kidney fibroblasts move cells toward the tertiary-lymphoid-tissue organiser-like state? | **Inconclusive** | The housekeeping gene GAPDH fell inside the neutral range in only 1 of 3 seeds, failing pre-registered rule 2 | **The primary metric was measuring distance moved, not direction.** Across control genes, approach-to-goal and approach-to-an-unrelated-state correlated at 0.946 / 0.913 / 0.915. There was no positive control on the goal axis at all |
| EXP-20260927-004 | 2026-09-27–28 | If the primary metric is changed to a goal-specific component, can we find a gene that moves cells along the goal axis? | "Suggested but not established" (a methods experiment, not a hypothesis verdict; wording weakened by the PI) | One of 10 candidates (SYN3) met the rule in all three seeds — but so did RASA1, a control gene expected to have no effect (2.76 / 2.91 / 2.79 SD). The measured chance pass rate was 1/82 = 0.012, so 10 candidates yield 0.12 passes by chance | A single gene cannot be distinguished from chance. Before reporting a pass, apply the same rule to the control panel and measure the chance passes empirically |
| EXP-20260928-001 | 2026-09-28 | If the unit of judgement changes from one gene to a set of 24, does the goal axis acquire a scale? | Qualified pass (a methods experiment, not a hypothesis verdict; PI decision) | — | Set mean 0.884 / 1.043 / 0.903 SD, the top rank among 200,000 random sets. But the top 3 genes (CEBPD +6.7, CP +5.0, NAMPT +3.2) accounted for 63–70% of the mean. **The single-gene detection floor of 2.73 SD does not fall with more cells or more seeds** |
| EXP-20260928-002 | 2026-09-28–29 | Does overexpressing each of 120 target genes of the LTβR–non-canonical NF-κB pathway move cells toward the goal? | **Refuted** (confirmed by the PI, 2026-09-29. **The run of inconclusive verdicts is broken here**) | — | Set mean 0.157 / 0.139 / 0.227 SD (threshold 0.213; absolute condition 0.245); all three seeds failed. The correct reading is not "no effect" but "no effect of 0.37 SD or larger was detected." **What we could measure was the pathway's entrance and plumbing; its outputs — chemokines and cytokines — dropped out of single-nucleus data** (only 13 of 63 candidates were feasible) |
| EXP-20260929-001 | 2026-09-29–30 | Does overexpressing a 24-gene secretory programme move fibroblasts toward a state high in tertiary-lymphoid-tissue-forming markers? | **Verdict withheld** (PI decision: neither supported, refuted, nor inconclusive) | The positive control on the goal axis — the 7 marker genes used to define the goal itself — reached only 0.718 / 0.659 / 0.674 against a required 0.861–0.896. **With the instrument not working, the primary values (0.024 / −0.036 / 0.053) are not interpreted** | Only 7 genes were used as positive controls. 227 candidates met the eligibility conditions at the time. K = 7 was a choice, not a constraint |
| EXP-20260930-001 | 2026-09-30–10-01 | Do 11 tertiary-lymphoid-tissue organiser genes move fibroblasts toward the "activated stromal state" as defined by the original authors (Lake et al. 2023, PMID 37468583)? | **Inconclusive** (metric-soundness checks failed) | Three control-gene-only regression diagnostics failed: housekeeping gene B2M at 0.003–0.010 (out of range); correlation with the uniform component 0.4006 in seed 81 (threshold 0.40); correlation with displacement magnitude −0.122 / −0.148 / −0.128 (threshold 0.10) | **For the first time in the project, the goal-axis positive control passed** (71-gene set mean 1.002 / 1.093 / 1.067 against a control 99th percentile of 0.291–0.294). The primary result 1.339 / 1.314 / 1.288 missed the absolute condition of 1.36 by 0.02–0.07. With a detection floor of 2.00 SD, "no effect" cannot be written |
| EXP-20261001-001 | 2026-10-01–02 | What happens if the intervention set is enlarged from 11 to 119 genes so that the detection floor falls below the observed effect sizes? | **Inconclusive** (one soundness check failed). **Note**: the experiment itself could not decide, but the same hypothesis (HYP-20261001-01) was later settled by EXP-20261002-002 with its enlarged control panel, and **the hypothesis status became "refuted"** (confirmed by the PI, 2026-10-03) | Four of five diagnostics passed; only the correlation between residual and displacement magnitude failed, at −0.101 / −0.152 / −0.109 against a threshold of 0.10 (seed 81 missed by 0.0013) | **For the first time, the primary comparison was measured with a working instrument and adequate power.** Detection floor 0.604 SD; primary result 0.208 / 0.098 / 0.142. A prediction made on the older data (−0.055 / −0.096 / −0.078) did not transfer to the production run with a redrawn control pool |
| EXP-20261002-001 | 2026-10-02 | (Re-evaluation, no perturbation. A re-reading of HYP-20261001-01, not a verdict on a new hypothesis) Are the five soundness checks valid in the first place? If not, rebuild them and re-read the same data | One soundness check failed (upper confidence bound exceeded; at this point the hypothesis status was left at `inconclusive` and not changed. It was later settled as "refuted" by EXP-20261002-002) | The point estimate of the bias was 0.0364 / 0.0444 / −0.0130, far below the budget of 0.10 SD in all three seeds. What failed was the upper end of the 95% confidence interval: seed 82 came in at 0.1017, **an overshoot of 0.0017** | **The old five checks passed only 7 of 60 (12%) "clean" metrics that cannot contain confounding by construction.** If a single run fails with probability 0.88, six failures occur with probability 0.88⁶ = 46%. The inconclusive verdicts to that point were within the range expected even from a perfectly sound metric |
| EXP-20261002-002 | 2026-10-02–03 | Does doubling the control panel from 500 to 1,000 genes tighten the interval enough to earn the right to state a result? | **Refuted: no effect of 0.604 SD or larger was detected. Gate G3 confirmed by the PI on 2026-10-03** | **It did not stop. For the first time in this project all four soundness checks passed and a verdict was reached** | 95% upper bounds 0.0859 / 0.0970 / 0.0762 against a budget of 0.10. The primary result, 0.261 / 0.146 / 0.199, reaches neither the absolute condition of 0.41 nor the detection floor of 0.604. Enlarging the control panel moved the primary result by +0.048 to +0.057 and the positive control by +0.092 to +0.105, both within the range predicted in advance |

**Note 1.** EXP-20260925-001 was a shake-down run for the environment and procedures, not a test of a hypothesis, and is omitted from the table.

**Note 2.** Five of the rows marked "Inconclusive" above (EXP-20260926-001, -005, EXP-20260927-001, EXP-20260927-003 and EXP-20260930-001) correspond one-to-one with the five inconclusive hypotheses in the opening table. The sixth such row, EXP-20261001-001, ended inconclusive as an experiment, but **its hypothesis (HYP-20261001-01) was settled as "refuted" by EXP-20261002-002**, so in the opening table it sits on the "refuted, 4" side. EXP-20261002-001 is a re-evaluation of the same hypothesis and is not counted separately.

**Note 3.** EXP-20260927-002, EXP-20260927-004 and EXP-20260928-001 were experiments on the method itself and carry no hypothesis verdict. They are not among the 10.

**Note 4.** Only the PI may set a verdict. Every verdict above has been confirmed by the PI.

---

## A.2　The walls we hit repeatedly

Grouped by what the obstacle was, not by date. There were at least six.

### Wall 1　Deletion moves every gene in the same direction

In the first experiment (EXP-20260926-001), whatever gene we deleted, the cell moved away from all three state centroids by a comparable amount — in every one of 90 control genes × seeds. About 56% of this is a geometric second-order term: any displacement of a vector takes it further from every nearby point.

The uniform movement can be subtracted. Doing so brought the fraction of control genes moving away from the goal down from 1.00 to 0.42–0.50. But subtraction alone did not establish the positive controls. What finally removed the uniform movement was fine-tuning: in the fine-tuned model the fraction fell to 0.10–0.15. **That only eliminated the uniform movement, however; the separate problem that the positive controls would not stand up remained** (HAVCR1 failed in EXP-20260926-003).

### Wall 2　The metric measured distance, not direction

Overexpression carries a component that moves a cell closer to *every* state centroid. In EXP-20260927-003, measured across control genes, the correlation between approach-to-goal and approach-to-an-unrelated-state was 0.946 / 0.913 / 0.915 across three seeds. **If approach to a single centroid is the primary metric, what is being measured is not "did it head for the goal" but "how far did it move."**

The fix was to subtract the mean approach to the other three centroids and divide by the displacement magnitude (EXP-20260927-004). This dropped the correlation with the uniform component from 0.979 to 0.106. It also had a side effect: the share of control-gene scatter that is gene-specific and reproduces across seeds rose from 32% to 82%, so **a rank-only decision rule no longer reduced its chance pass rate by adding seeds** (probability of at least one pass among 10 candidates: 0.665). We added an absolute size condition — at least 2.5 control standard deviations — to the rule.

### Wall 3　A single gene does not provide enough power

Pre-analysis for EXP-20260928-001 showed that, to first order, the primary metric equals the inner product of a gene's mean displacement direction with a fixed axis determined by the centroids (direction alone explains 81–85% of the variance). That displacement direction is near-fixed across seeds (cosine 0.985, minimum 0.908).

In other words, **the spread of the control genes is not statistical noise but a scatter of directions on a sphere, so it does not shrink with more cells or more seeds.** Empirically, doubling the cells per gene from 100 to 200 shrank the control standard deviation only 5%, from 3.19 × 10⁻⁵ to 3.04 × 10⁻⁵.

**The single-gene detection floor is 2.73 SD, and there is no way to lower it.** The only remaining move was to raise the unit of judgement. Judging a set of genes brings the floor to 0.54 for 10 genes and 0.32 for 24. The 10 positive-control candidates of EXP-20260927-004 already carried an effect of 1.17 / 1.17 / 1.00 as a set. **They had simply been undetectable one gene at a time.**

### Wall 4　Housekeeping genes refuse to be neutral

Housekeeping genes (ACTB, GAPDH, PPIA) serve as negative controls on the premise that, since they operate uniformly in every cell, perturbing them should not move the cell state. **That premise broke in four experiments** (EXP-20260926-001 and -005, EXP-20260927-003, EXP-20260929-001). Each time, the pre-registered rule forced an inconclusive verdict.

We explained the cause wrongly once. On 2026-09-27 we wrote that it was "because there is more room to push them up," but the three genes sit at nearly the same original relative position (ACTB 0.989, GAPDH 0.981–0.984, PPIA 0.972–0.974), so that explanation does not hold. **The real cause was the magnitude of the perturbation-induced displacement.** The three genes' displacements sit at the 97.8th–99.8th percentile of the control distribution, and the metric retained a correlation of +0.171 with displacement.

As a remedy we added log displacement to the de-confounding regression. On the EXP-20260929-001 data it worked: the three-gene set mean went from 1.071 to 0.335, and the number of seeds inside the neutral range went from 0/3 to 2/3. **But the correction became unnecessary once the design changed — and then made things worse.** In EXP-20260930-001, with no correction at all, ACTB sat at 0.834 / 0.840 / 0.863 and GAPDH at 0.753 / 0.740 / 0.792, both inside the range. The real reason they returned to neutrality was the rebuilt goal definition and subtraction scheme, not the correction. Leaving the correction in degraded two other diagnostics. We removed it in the end.

B2M, by contrast, stayed off-scale (0.0027–0.0054) no matter how the metric was rebuilt. As the light chain of MHC class I it rises sharply with inflammation, so it may simply have been unsuitable as a housekeeping control. An adoption rule fixed before looking at the results relegated it to descriptive reporting.

### Wall 5　The goal-axis positive control would not stand up

This happened four times running. **A positive control that will not stand up means there is no gene known to move cells along the goal direction — that is, the ruler has no scale on it.** Without a scale, a small primary value cannot be told apart from an unmeasurable one.

For these four runs, Labretary initially explained that "the goal is too close to measure." **That explanation was wrong.** The ratio of start-to-goal distance over start-to-myofibroblast distance lies between 0.442 and 0.637 across the four experiments. EXP-20260928-001, which earned a qualified pass, sits at 0.452; EXP-20260929-001, which failed, sits at 0.442 — **essentially the same. This quantity did not separate success from failure.**

The real cause was **the number of genes used as positive controls**: 1 (EXP-20260927-004), 0 (EXP-20260928-002), 7 (EXP-20260929-001) and 24 (EXP-20260928-001). As Wall 3 explains, one gene gives a detection floor of 2.73 SD and seven genes still 1.71 SD. **No positive control in this project ever reached 1.71 SD.** And at the time of EXP-20260929-001 there were 227 eligible candidates. Seven was a choice, not a constraint.

In response we built a pre-commitment feasibility check (`scripts/analysis/axis_feasibility.py`) that decides, **without running a single perturbation**, whether the design can measure the intended axis, and fails any design with fewer than 20 positive-control genes. Applied retrospectively to the four earlier experiments it returns three failures, one hold, and zero passes. **This project had never once had a design in which the goal axis was unconditionally measurable.**

The breakthrough came with EXP-20260930-001. We stopped defining the goal with a home-made score, defined it instead from the original authors' cell annotations (Lake et al. 2023, PMID 37468583), and took 71 positive controls from the same paper's supplementary tables. Their median detection rate in the goal cells in single-nucleus data was 0.321 — four to ten times the 0.039–0.088 of the previous four experiments. The result: **the 71-gene positive control set reached 1.002 / 1.093 / 1.067 SD, and for the first time the goal axis had a scale.** All six sensitivity analyses agreed.

### Wall 6　The acceptance criteria themselves were not valid

This was the last and largest problem to surface.

Once six hypotheses had been declared inconclusive because one of the pre-registered soundness checks had failed — **the count as it stood on 2026-10-02** — the PI directed us to audit the checks themselves. (One of those six, HYP-20261001-01, was later settled as refuted by EXP-20261002-002; every figure in this section is as of the audit on 2026-10-02.) We built 60 "clean" metrics — displacement vectors projected onto random directions in the 1,152-dimensional embedding space, which cannot contain confounding by construction — and ran them through all five checks.

| Check (internal name) | What it measures | Passes out of 60 |
|---|---|---|
| M1 | Calibration of the control rank distribution | 60 / 60 (never fires on any defect; carries no information) |
| M2 (original rank, input length) | Rank correlation of the residual with these | Undetectable in principle: the de-confounding regression is a full quadratic in them, so they are absorbed by definition. Injecting a ±0.30 defect does not change the statistic to four decimal places |
| M2 (relative position) | Rank correlation with relative position | 58 / 60 |
| M3 | Ranks of the four individual neutral controls | 43 / 60 |
| M4 | Rank correlation with the uniform component | 49 / 60 |
| M5 | Rank correlation with displacement magnitude | **8 / 60 (13%)** |
| All five together | — | **7 / 60 (12%)** |

**We had been using a test battery that rejects 88% of confounding-free metrics.** At a per-run failure probability of 0.88, six failures occur with probability 0.88⁶ = 46%. **The six inconclusive verdicts standing on 2026-10-02 were entirely consistent with a perfectly sound metric.**

The thresholds were also out of balance. For the same exposure (a 0.9 quantile within the controls), the distortion imposed on the set mean, expressed against the decision line of 0.41 SD, was:

| Check | Threshold | Distortion at threshold ÷ 0.41 |
|---|---|---|
| M5 (displacement) | 0.10 | 36% |
| M4 (uniform component) | 0.40 | 104% |
| M2 (relative position) | 0.25 | 151% |

**M5 alone was excessively strict, while M4 and M2 were letting through distortions large enough to overturn the conclusion** — and neither M4's nor M5's threshold had any recorded derivation.

They also overlapped. M5 and M2 (original rank) correlate at 0.978: the same thing measured in different coordinates. The reason is that perturbation displacement is nearly the same quantity as "how low the gene originally sat" (rank correlation 0.873–0.878 across 500 controls). M2 is absorbed exactly by the regression and therefore always passes, while M5, which measures the same content, always fails.

Finally, M5's threshold of 0.10 is only **2.1 times** the statistic's own standard error of 0.045–0.047. That is the level at which the pass/fail outcome flips merely by redrawing the control panel — and that is exactly what happened (A.3, item 6).

The remedy was to **re-derive the checks from the decision they are meant to protect**. Instead of asking "how large is this rank correlation," we now estimate directly: **could confounding alone move this set mean by enough to matter against the 0.41 SD decision line?** Under the rebuilt check, the pass rate on clean metrics rose from 12% to **94.6%**, and the budget of 0.10 is 24.4% of the 0.41 decision line, 1.25 times the irreducible noise (SD 0.080 for a 119-gene set mean) and 13.7 times the statistic's standard error of 0.0073 — **large enough to be worth defending and small enough to be measurable.** Sensitivity was confirmed too: the pass rate under conditions where the true distortion exceeds 0.15 is 0.000.

---

## A.3　Statements that were wrong, and their corrections

What follows is a list of explanations Labretary (the AI research secretary) gave the PI that later proved wrong. **For the sake of scientific integrity, each is recorded together with how the error came to light.**

### 1　"VCAM1 points the opposite way from expectation" (2026-09-27)

- **What we said.** VCAM1, placed as a directional control, came out against expectation in all three seeds.
- **What was actually true.** Only 7, 13 and 15 comparable control genes could be found, below the 20 required by the primary test. All three seeds therefore produced *no value*, not a reversed one. The uncorrected raw values were negative in all three seeds (−2.5 × 10⁻³, −9.9 × 10⁻⁴, −2.8 × 10⁻³), i.e. moving away from the goal exactly as expected.
- **How we found out.** After reporting, Labretary rechecked and found that the evaluation code had been counting missing values as failures. Guardrails were added: never count a missing value as a failure; confirm that a value could be computed before reporting a control's pass or fail.

### 2　"We found one gene that moves cells along the goal axis" (2026-09-28)

- **What we said.** In EXP-20260927-004, SYN3 met the rule in all three seeds out of 10 positive-control candidates; we had found a gene that moves cells along the goal axis.
- **What was actually true.** **RASA1, a control gene expected to have no effect, met the identical rule in all three seeds** (2.76 / 2.91 / 2.79 SD). The measured chance pass rate was 1/82 = 0.0122, not the 0.000 we had estimated in advance. With 10 candidates, 0.12 chance passes are expected; the single observed pass is indistinguishable from chance.
- **How we found out.** Labretary applied the same rule mechanically to the control panel and recomputed. A guardrail was added: apply the rule to the controls and measure chance passes empirically before reporting any pass.
- **A related error in the same experiment.** A script created outside the pre-registered set had entered the primary-metric computation, and we wrote that "the interpretation of the results is unchanged." That too was wrong: removing the centroids produced by that script flips diagnostic M5 from pass to fail (SYN3 3.19 → 2.41 SD).

### 3　"The goal is too close to measure" (2026-09-30)

- **What we said.** The start-to-goal distance is only 1/2.2 of the start-to-myofibroblast distance; the goal is too close, so the signal is intrinsically too small to measure.
- **What was actually true.** Recomputed across the four earlier experiments, the axis-length ratio ranges only from 0.442 to 0.637 and **does not separate success from failure**: the qualified pass (EXP-20260928-001) sits at 0.452 and the failure (EXP-20260929-001) at 0.442. The real cause was the number of positive-control genes (1 / 0 / 7 / 24). The related claim that "redefining the goal made the axis 14% shorter" was also wrong — it compared **two different quantities**, a centroid-to-centroid distance against the length of a direction vector. Compared like for like, the change is only −3% to −9%.
- **How we found out.** Recomputing the quantity across four past experiments (without running any perturbation) in order to set thresholds for the feasibility check showed that it did not discriminate. The check's primary failure condition was changed from the axis-length ratio to the number of positive-control genes (fail below 20).

### 4　"Gate G2 approval is required" (2026-09-28)

- **What we said.** A plan costing 21.8 hours and using 4,616 cells is large enough to require the PI's planning approval (G2).
- **What was actually true.** The G2 trigger conditions are "over 24 hours," "over 10,000 cells," or "a genome-wide screen," and **none applied.** We had passed on a delegate's wording without checking it against the conditions.
- **How we found out.** The PI pointed it out ("isn't anything under 24 hours something you just proceed with?").
- **A second instance (2026-10-01).** A delegate reported that "the cumulative cell count far exceeds 10,000, so G2 is required." In fact the cell count `run_isp.py` inspects is cells *per gene* (200 in this design), not a cumulative total; the structurally identical EXP-20260930-001 had 360,000 cumulative cell-gene evaluations and was flagged as not requiring approval. This time Labretary re-applied the formula and **stopped the unnecessary approval request before it reached the PI.**

### 5　"The housekeeping genes returned to neutrality thanks to the displacement correction" (2026-10-01)

- **What we said.** In EXP-20260930-001, ACTB and GAPDH returned to the neutral range for the first time in four experiments because of the displacement correction introduced in stage 2.
- **What was actually true.** On this experiment's data, **they were inside the range even with no correction at all** (ACTB 0.834 / 0.840 / 0.863; GAPDH 0.753 / 0.740 / 0.792). The reason was not the correction but the rebuilt goal definition and subtraction scheme. The correction genuinely worked on the EXP-20260929-001 design (1.071 → 0.335), but in the new design it was unnecessary and actively degraded other diagnostics.
- **How we found out.** A post-hoc diagnostic laid several correction forms side by side (none, linear, full quadratic, and others). The finding amounted to a proposal to remove a correction Labretary itself had introduced; the PI approved it.

### 6　"Switching to a projection onto the goal direction will make M5 pass" (2026-10-01)

- **What we said.** If the primary metric becomes "project the displacement vector onto the start-to-goal direction" and the displacement correction is dropped, M5 becomes −0.055 / −0.096 / −0.078, inside the 0.10 threshold.
- **What was actually true.** In production (EXP-20261001-001) it came out at **−0.101 / −0.152 / −0.109, failing in all three seeds**; seed 82 was worse than predicted. The prediction had been made on the old data with the old control pool and did not transfer once the controls were redrawn.
- **How we found out.** The production result contradicted it, and the subsequent audit identified the cause. **The embeddings of EXP-20260930-001 and EXP-20261001-001 are bit-identical** (maximum difference 0.0 for centroids and for pre- and post-perturbation vectors alike). The only difference is the membership of the 500 control genes: 203 in common, 59.4% replaced. That alone moved the statistic by 0.05 — essentially the same as its standard error under control-pool redraw (0.045–0.047). **The problem was not the metric; the check's measurement precision was mismatched to its threshold.**

### 7　"Adding controls will make it pass" (2026-10-02)

- **What we said.** Enlarging the control panel beyond 500 genes will tighten the confidence interval and the check will pass.
- **What was actually true.** The correct statement is "**the probability of passing goes up.**" Two things had been missed. (1) The analytic standard error understates the real variability under control-pool redraw by a factor of **1.28** (1.20 / 1.34 / 1.31 by seed), because the formula counts only the uncertainty in the slopes fitted on the controls and ignores that redrawing the pool also moves the standardisation, the quantile scale, and the residuals themselves. (2) Adding controls moves the point estimate too (SD 0.027 for 500 added genes). With both effects included, the pass probability is 0.72 with 500 added, 0.83 at a total of 2,500, and **plateaus at 0.92 even with infinitely many**. The three seeds' deviations correlate at 0.82, so adding seeds adds little independent information.
- **How we found out.** The experiment-planning delegate recomputed and corrected it. The plan now contains an explicit section denying the expectation that "adding more will eventually guarantee a pass."

### 8　Pulling numbers from the post-hoc diagnostic table, we read the wrong row (2026-10-01)

- **What we said.** We reported to the PI, and wrote into the plan template, that adopting "projection onto the direction, without the displacement term" would give M5 values of −0.055 / −0.096 / −0.078.
- **What was actually true.** Those values belong to the row "projection, no displacement term (but with original-rank and input-length corrections)." The method actually specified in context was "projection, no correction at all (intercept only)," whose values are +0.010 / −0.048 / −0.013. **The method we prescribed and the numbers we quoted came from different rows.**
- **How we found out.** A delegate reproduced the table and spotted it. A guardrail was added: transcribe the variant name alongside each number and confirm it matches the prescribed method before reporting.

### 9　The handover bundle shipped with an empty `scripts/` directory (2026-09-30)

- **What we said.** We had packaged the deliverables, the decision history and the next steps into a single zip so that another operator could take over, and had verified its integrity.
- **What was actually true.** `scripts/` inside the zip contained **zero files.** `find -exec cp --parents` fails under the macOS `cp`, but `find` still exits 0, so the fallback we had prepared never fired. Our verification consisted only of the archive integrity test (`unzip -t`); we never reconciled file counts directory by directory.
- **How we found out.** It surfaced after the bundle had reached the PI. A guardrail now requires reconciling per-directory file counts against the source before archiving, and treats a count of zero as an anomaly.

### 10　"There was no effect on the computation itself; it only cost 43 minutes of delay" (2026-09-30)

- **What we said.** The seed 82 → seed 83 handoff died and delayed the start of seed 83 by 43 minutes, but the computation itself was unaffected.
- **What was actually true.** The PI's response was blunt: that is not an acceptable way to think about it. The completion watcher was built to fire "when all three seeds have finished," so **if the handoff dies it never fires at all**. We noticed only because the PI asked about the status; **there was no upper bound on how long the lapse could have gone unnoticed.** The figure of 43 minutes is the loss that happened to be observed, not a bound the system guarantees.
- **How we found out.** The PI pointed it out. From then on, failures are reported in terms of **the upper bound on the time we could have failed to notice** and **why we noticed at all**, and the watcher was rebuilt to detect *lack of progress* rather than completion.

### 11　"CXCL13, CCL19 and CCL21 are not in the Geneformer vocabulary" (2026-09-29)

- **What we said.** These three chemokine genes are absent from the model's token dictionary and therefore cannot be perturbed.
- **What was actually true.** Querying the 20,275-entry dictionary directly showed all six genes present (CXCL13 = 9811, CCL19 = 13095, CCL21 = 7298). **The real reason for the zero count was that the starting state had been defined as "cells with zero counts of these three genes"** — the absence was a logical consequence of the state definition.
- **How we found out.** Labretary queried the dictionary directly. A guardrail was added: before writing "not in the dictionary" or "not in the data," look it up, and distinguish (i) absent from the dictionary, (ii) absent from the data, and (iii) **impossible given how the state is defined**. Confusing (iii) with (i) reports a non-limitation as a limitation of the measurement.

### 12　"Only one gene, COL5A2, overlaps between the positive controls and the control panel" (2026-09-30)

- **What we said.** Checking the 39 new positive controls against the existing 500 control genes found exactly one overlap, COL5A2.
- **What was actually true.** ANKS1B also overlapped. We had checked only against the 42-gene minimal list and **had not examined the additional genes in the 76-gene conservative list used for sensitivity analysis.**
- **How we found out.** The experiment planner ran an exhaustive cross-check and found it. The rule became: check overlaps exhaustively across *every* set that could enter the judgement, sensitivity-analysis sets included.

---

## A.4　Where the PI's judgement changed the course of the work

Decisions by the PI that redirected the research. None of them is a choice Labretary would have made on its own.

### 1　Don't cling to an unmeasurable question — rebuild it into a measurable one (2026-09-28)

When it became clear that three genes (LTBR, RELB, NFKB2) could not be judged for lack of power, the PI **changed the question**: from "do these three genes work?" to "does this pathway's target gene set work?" A three-gene set has a detection floor of 2.00–2.63 SD, barely better than the 2.73 SD of a single gene; the analysis indicated a minimum of 50 genes and a recommended 100–150.

Labretary's default reflex is to re-measure the same question with a better method. The PI declined that and chose instead to rewrite the hypothesis and seek fresh adoption. **From this point on, every judgement in the study was made at the level of gene sets.**

Relatedly, on 2026-09-29, after a refutation, the PI chose "find data where it can be measured and re-pose the question" over "close this line of inquiry" — consistently refusing to **close a direction merely because it could not be measured.**

### 2　Do not return to the substantive question while the positive control is unestablished (2026-09-28)

When EXP-20260927-004 showed that the positive control could not be distinguished from chance, the PI chose to **redo the establishment of the positive control** and deferred re-testing the substantive question (can a pathway move fibroblasts?).

Without that decision we would have gone on measuring the substantive question with an unscaled ruler. As Wall 5 describes, the goal axis did not acquire a scale until three experiments later, in EXP-20260930-001.

### 3　"Before anything else, lay out M1 through M5 again and check whether these tests are actually valid" (2026-10-02)

**This was the decisive turning point.**

Once six hypotheses had been declared inconclusive — the count as it stood on 2026-10-02 — each because a soundness check had failed, the PI wrote:

> Just doing things, with no purpose, and assembling some sort of result each time, is not progress. Start by laying out M1 through M5 again and checking whether these tests are actually valid in the first place.

**Throughout that accumulation of inconclusive verdicts, Labretary had responded each time by fixing the check that failed, and had never once questioned the validity of the checks themselves.** The work looked like progress while circling the same spot.

The audit that followed is Wall 6: a battery that rejects 88% of confounding-free metrics, thresholds with no recorded derivation, two checks correlated at 0.978 measuring the same thing, and a threshold only 2.1 times the statistic's own standard error. On that basis the PI **approved re-deriving the checks from the decision they protect** — accepting that a criterion the team itself had set earlier should be discarded once it is shown to have no basis.

Three general rules were recorded from this. (1) When stuck, audit the premises before adding more work. (2) Derive a check from the decision it is meant to protect; "we used it in earlier experiments" is not a justification. (3) When the same kind of failure recurs three or more times, suspect the framework, not the individual failures.

### 4　Improve the measurement rather than relax the criterion (2026-10-02)

When even the rebuilt check failed — the point estimate of the bias was far inside the budget, but the upper end of the 95% interval overshot by 0.0017 in one seed — the PI chose to **enlarge the control panel, tighten the interval, and measure again.** Relaxing the threshold by 0.0017 would have passed it; that option was not taken.

This is the same kind of response as enlarging the intervention set from 11 to 119 genes to bring the detection floor from 2.00 down to 0.604 SD. **Even knowing the conclusion would be a refutation, invest in earning the right to state it.**

### 5　What held throughout: never move a criterion after seeing the result

From 2026-09-27 onward the PI made no exceptions.

- 2026-09-27: explicitly rejected the option of invalidating a single housekeeping gene post hoc in order to convert the verdict into a refutation.
- 2026-10-01: the primary result missed the absolute condition of 1.36 by only 0.02–0.07; the threshold was not moved.
- 2026-10-02: a check failed by 0.0013; the threshold was not moved.

**When a criterion turns out to be unsuitable, do not issue a conclusion — rebuild the design.** This is why the verdicts in this study can be trusted.

---

## A.5　Operational failures, and the machinery we installed

### 1　A completed run went unnoticed for roughly three hours (2026-09-28)

EXP-20260928-001 finished cleanly at 07:05:31 and the watcher missed it. The monitoring facility in use had a 30-minute ceiling and had to be re-armed each time it expired; after a re-arm at around 04:36, the next opportunity did not arrive for about two hours, leaving a gap in coverage.

**What we installed.** Labretary itself checks status once the estimated finish time has passed. Every long computation is logged in `tasks/long-jobs.md` with what is running, where, until when, and how to verify it — always including a command that can be run by hand.

### 2　Seed handoff depended on someone remembering, and cost 43 minutes (2026-09-30)

With two seeds running sequentially on one machine, a human was expected to launch the second when the first finished. The handoff failed and the machine sat idle for 43 minutes.

**What we installed.** An automated handoff (`scripts/remote/handoff.py`). It starts the next seed only when the previous one satisfies **all three** conditions: the expected gene count reached, zero errors, and the execution session terminated. If any one is missing it **does not start**, writes the reason to the log, and halts — so that it can never silently jump over a stop rule.

### 3　The completion watcher fired only "when all three are done," so a dead handoff meant it never fired at all (2026-09-30)

We noticed the 43-minute delay above only because the PI asked. The watcher was built to signal when all three seeds had finished; if the handoff dies, that condition is never met. **There was no upper bound on the time the lapse could have gone unnoticed.**

**What we installed.** A watcher that detects **absence of progress** rather than completion (`scripts/remote/watch_run.py`). It alarms if the combined completed-gene count across the three seeds **fails to increase for 25 minutes** — catching a dead handoff, an abnormal termination, or a network outage alike.

It earned its keep on 2026-10-02: after the handoff failed at 04:46, it alarmed at 05:12 with "1,434 of 2,151, no increase for 25 minutes," and investigation began at 05:21.

### 4　We wrote the watcher in an old bash and crashed it twice (2026-09-30)

Writing the watcher for the stock macOS `/bin/bash` (version 3.2) produced two defects — no associative arrays, and variable names misparsed when Japanese text immediately follows a variable reference — and it died on startup twice. Because output was being discarded, neither crash could be diagnosed.

Labretary explained this as "the environment is old," and the PI's reply was that on a Mac the shell is zsh, and a mistake in bash is not an excuse. **This machine's default shell is zsh, every computation in this project is in Python, and there was no reason to choose bash 3.2. We chose the tool.**

**What we installed.** The watcher was rewritten in Python and three bash versions were retired. Anything long-running is written in the project's Python environment. Standard error is never discarded; it goes to `logs/*-stderr.log`. Thirty to sixty seconds after launch, we verify both that the process is alive and that its log is growing.

### 5　After 8.6 hours of computation, the evaluation code crashed (2026-10-01)

The perturbation run for EXP-20260930-001 (three seeds, 8.6 hours total) completed normally, and then evaluation stopped on a missing column (`KeyError: 'relpos_med'`). The cause was Labretary's own refactor in stage 2, which failed to carry one diagnostic column across when a function was extracted into a shared module. The module's own self-test passed, but nothing checked what the caller expected. Worse, the smoke run had only verified that perturbation executed; **the evaluation code had never been exercised at all.**

**What we installed.** Running the evaluation code end-to-end on the smoke-run output is now a precondition for starting a production run. When extracting a function, reconcile the required-column list against the self-test.

### 6　A mid-experiment file push broke the link to the model (2026-10-02)

A configuration file needed to launch seed 83 had not reached the machine, so Labretary pushed the missing file with the transfer helper (`isp_remote.sh push`). That helper re-points the model symlink to a fixed location; on that machine the weights live elsewhere, so the link dangled and the launch failed. Recovery took about 39 minutes.

**What we installed.** (a) When one machine runs several seeds, send **all** of that machine's configuration files up front and verify both their existence and their sha256 hashes. (b) **Never run a push in the middle of an experiment**; if something must be sent, verify the link target immediately afterwards. (c) Verify the model weights' sha256 at transfer time and stop if it does not match.

### 7　In the final experiment, this machinery fired twice (2026-10-02–03)

During EXP-20261002-002 the safeguards actually prevented failures.

- **A bad model link at transfer time was caught and the transfer stopped.** When a push would have pointed the model link at the wrong location, the sha256 verification on the weights detected it and halted on the spot. The last time Labretary stepped into the same trap (2026-10-02), nothing noticed until the computation died five hours later.
- **The seed 82 → seed 83 handoff connected automatically, with zero delay.** The previous attempt failed in two stages and cost 43 minutes, and was noticed only because the stall detector alarmed.

Both are cases of the machinery described in A.5 items 2, 3 and 6 working as designed.

---

## A.6　What to take away from this appendix

1. **None of the five inconclusive verdicts is a result that says "there was no effect."** In each, the ruler had no scale, or the test of the scale was itself invalid, so we had no right to interpret the number. They are different in kind from the four refutations, which are answers.
2. **When the same kind of failure recurs three or more times, suspect the framework rather than the individual failures.** We did not see it until the inconclusive count had reached six on 2026-10-02.
3. **Discard a criterion you set yourself once it is shown to have no basis.** Derive the check from the decision it protects, and set the threshold from both "how much would this check distort the conclusion when it sits exactly at threshold" and "can we measure it at that precision at all."
4. **Report a failure in terms of the upper bound on the time it could have gone unnoticed and the reason it was noticed — not in terms of the loss that happened to be observed.**

**Everything reported in this appendix is the output of in silico gene perturbation computations, not an experimentally confirmed biological fact.**
