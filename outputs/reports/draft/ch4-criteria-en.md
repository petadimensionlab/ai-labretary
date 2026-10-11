# Chapter 4: The Judgment Framework — All Decision Criteria

## 4.1 Why Decision Criteria Are Necessary

What this study measures comes down to a single operation. The gene-expression profile of a kidney fibroblast is converted by a foundation model, Geneformer, into a point in a several-hundred-dimensional embedding space (the cell embedding). A gene is then artificially forced to a high expression rank (overexpression), and we measure how far that point moves toward the "target cell state." This is not a laboratory fact; it is a computational simulation called in silico perturbation.

For this measurement to mean anything, several preconditions must hold:

- **The axis must actually exist.** The direction labeled "toward the target state" must be a meaningful direction in the model's embedding space, and this must be confirmed in advance by perturbing known driver genes (positive controls) and checking that they move along that direction.
- **Confounding must be absent, or at least measured and subtracted.** Factors with no biological meaning—how far a gene's expression rank had to be pushed, how many genes were detected in a given cell (input length)—must not be allowed to masquerade as the effect of interest.
- **Circularity must be avoided.** Genes used to define the target state must not be reused as "evidence that the effect occurred."
- **Statistical power (the detection floor) must be sufficient.** Even a real and meaningful effect cannot be detected if the gene set is too small; this is a design-time limitation known in advance, not a measurement error.

These preconditions break down in three main ways: (1) confounding, where a quantity treated as a control turns out to correlate with the outcome, so that bias is mistaken for effect; (2) circularity, where genes used to define the target are also used to measure approach to the target, producing a tautology; and (3) underpowering, where a set is too small to say anything except "inconclusive" or "no effect at this detection floor," regardless of whether a true effect exists.

This project has designed 17 experiments, of which 15 included running the perturbation computation (the remaining two were a pre-feasibility check and a re-evaluation, with no perturbation at all). Of these, 10 hypotheses reached an adjudicated verdict: **5 inconclusive, 4 refuted**, and 1 partially supported. That is, 5 of the 10 adjudicated verdicts ended "inconclusive." Tracing the cause back, in most cases the fault lay not in the measured quantity itself but in **the soundness-check criteria used to validate the measurement.** When the PI directed a retrospective audit of the old checks (M1–M5, discussed below) on 2026-10-02, it turned out that even on a "clean metric"—one that by construction contains no confounding—the full battery of checks passed only 12% of the time. In other words, a perfectly sound metric would fail this battery 88% of the time, and the project had been using exactly this battery across all six inconclusive verdicts **standing as of the audit date, 2026-10-02** (one of which was later resolved to "refuted" after the control pool was enlarged).

The principle this audit finally arrived at is this:

> **Criteria are derived from the judgment they are meant to protect. The value of the primary metric is never used to justify the criteria.**

Before this, neither the choice of checks nor their numeric thresholds had any traceable derivation beyond "this is what the previous experiment used." For some thresholds, no derivation was ever recorded. The rebuilt checks (S-1 through S-4, below) all derive their numbers from a single question fixed in advance: "how much could this set's mean value move from confounding alone, and how much of the judgment threshold we are trying to protect (the set's absolute-value condition of 0.41 standard deviations) would that erosion consume?"

## 4.2 Definition of the Primary Metric

### 4.2.1 Final Form

The primary metric is defined as follows:

> Take the displacement vector of the cell embedding before and after perturbation (post-perturbation point minus pre-perturbation point). Project it onto the unit vector pointing from the centroid of the starting state to the centroid of the target state. Divide the length of that projection by the magnitude of the displacement vector itself (the amount of movement).

This is equivalent to the cosine of the angle between the displacement and the direction toward the target. A larger value means a larger fraction of the perturbation-induced movement points toward the target.

### 4.2.2 Why Not the Difference Between Two Centroid-Approach Scores

An earlier design used the difference between "approach to the target centroid" and "approach to the starting-state centroid" as the primary metric. Empirically, this form had two problems.

**First, the sign changes depending on which state is subtracted.** Perturbations carry a uniform "pulling strength" toward any centroid, regardless of which gene is perturbed, and this strength differs by state (observed: +0.09–0.12 toward the target centroid, +0.18–0.22 toward the starting-state centroid, +0.52 toward the centroid of degenerated medullary fibroblasts). Simply changing what is subtracted—from the average of three states to the starting state alone—flipped the rank correlation between the residual and the movement magnitude from +0.009–+0.049 to −0.122 to −0.148. The subtraction operation itself reverses the sign of the result depending on an arbitrary choice of reference.

**Second, it correlates arithmetically with the neutrality metric.** The difference-form primary metric and the uniform component (a diagnostic quantity meant to capture perturbations that approach every centroid equally) are built from the same two underlying quantities: approach to the target and approach to the starting state. The difference (target − start) and the sum/average (the uniform component) are linear combinations of the same two variables, so they are mathematically guaranteed to correlate. The observed correlation between the two underlying quantities (0.937–0.944) and their variance ratio (1.23–1.25) predict a theoretical correlation of 0.28–0.30, which matched the observed correlation to four decimal places. **This is not confounding, so it cannot be removed by confound correction.** Building a metric that structurally correlates with its own soundness check—in a way no correction can remove—undermines the soundness check itself.

Because of these two defects, the primary metric was changed to the projection form, and the difference form was dropped.

### 4.2.3 Removing Confounding

Two factors with no biological meaning can leak into the primary metric:

- **Original rank** (where the gene sat in the cell's input ranking before the overexpression operation). Because overexpression moves a gene to the top rank, genes that started lower have to move further—an artifact of the mechanics of the operation, not of biology.
- **Input length** (roughly, the number of genes detected and tokenized for that cell). Cells with more detected genes can behave differently under the same operation.

To remove these, the primary metric uses the **residual from a regression fit on control genes** (a large set of genes assumed to have no effect). The regression's predictors are a full quadratic in original rank and input length. For genes under evaluation, the same fitted equation (with coefficients estimated from controls only) is used to compute an expected value, which is subtracted to obtain the residual. One rule is essential here: **a gene being evaluated is always excluded from the control pool used to fit its own regression.** Including the gene in its own reference set could artificially neutralize it through self-reference.

### 4.2.4 Why Movement-Magnitude Correction Is Not Applied by Default

Movement magnitude (‖Δ‖, how far the embedding actually moved) is itself a candidate confounder. Adding it as a regression predictor ("movement correction") is **not applied by default**, because its effect depends on the rest of the design.

In an earlier design (EXP-20260929-001, using the difference-form metric and subtracting the average of three states), adding movement correction pulled the set-average value of housekeeping genes (expected to be neutral, with no reason to approach the target) from 1.071 down to 0.335—exactly the intended correction. Taken alone, this looks like strong evidence that the correction works.

However, in a later design (EXP-20260930-001, using the projection-form metric), **housekeeping genes were already within the neutral range without any movement correction**, and adding it made the soundness checks worse instead (the rank correlation with the uniform component rose from 0.31 to 0.40, and the rank correlation with movement magnitude worsened from −0.07 to −0.12). In other words, an improvement that worked under one design did not generalize to another. **Whenever the design changes—the form of the primary metric, or which states are subtracted—whether a given correction is still needed must be re-measured, not assumed.** Because of this lesson, movement correction is off by default, and any decision to include it must be accompanied by a stated reason and measured evidence.

### 4.2.5 Why the Uniform Component Is Built from States Not Used in the Primary Metric

The uniform component (a diagnostic quantity capturing perturbations that approach every state's centroid equally) is built as the average approach toward centroids of states that are *not* used in the primary metric's projection (here, degenerated fibroblasts and degenerated medullary fibroblasts), rather than from the starting and target states. This avoids the problem described in §4.2.2: building the primary metric and the uniform component from the same two underlying quantities guarantees an arithmetic correlation. If the two shared the same ingredients, their correlation would be a mathematical artifact, not a meaningful soundness signal.

### 4.2.6 Why Cells with Zero Movement Are Excluded

Overexpression perturbation raises a gene's expression rank to the top of the input. If a gene already occupied the top rank before perturbation, the operation changes the input by not a single bit, and the pre- and post-perturbation embeddings are identical. The displacement vector's magnitude (movement) is then exactly zero.

Empirically, 114 of 338,967 cells (0.034%) fell into this category, and **in every one of them, the perturbed gene's original rank was already at the top (rank 1).** This is an inevitable consequence of the perturbation mechanism, not a flaw in the metric. But mixing zero-movement cells into the primary metric either leaves the cosine undefined (division by zero) or, if treated as zero, dilutes the measured effect uniformly. Zero-movement cells are therefore excluded as cases where "no perturbation occurred." It is essential that the same set of excluded cells be used across every version being compared (with or without a given correction, etc.); otherwise, a difference in exclusion criteria would masquerade as a difference between versions.

## 4.3 Soundness Checks S-1 through S-4 (The Final Set)

The four soundness checks (S-1 through S-4) all answer the question "is it valid to interpret the primary metric's value (the set mean)?" **All four must pass before the primary metric's value is used as grounds for any reported conclusion.**

### S-1: Bias Budget

**Definition.** Five candidate confounders are enumerated: movement magnitude (median ‖Δ‖ per gene), original rank (log2(1+rank)), input length (log2(input length)), relative position (rank ÷ input length), and the uniform component (per-gene residual). Using control genes only, the standardized residual z is regressed **jointly** on the quantiles (rank-transformed to 0–1) of all five candidates.

For the set under evaluation (e.g., the 119-gene intervention set), the "bias" contributed by each candidate v is

> bias b_v = (joint-regression partial coefficient fitted on controls) × (the set's mean quantile for v − 0.5)

and the total bias is B = Σ b_v. The set passes if **the upper bound of the 95% confidence interval for |B| (i.e., |B| + 1.96 × standard error)** is at most the budget of 0.10 standard deviations, in all three seeds.

**Why joint regression rather than marginal regression (fitting each candidate separately).** The five candidates are strongly intercorrelated (rank correlations among controls, observed range 0.44–0.88; e.g., movement vs. original rank 0.873–0.878, movement vs. relative position 0.780–0.787, movement vs. uniform component 0.436–0.459, uniform component vs. original rank 0.237–0.264). Summing marginally fitted b_v values counts the same underlying bias multiple times. On synthetic data, **the marginal-regression Σ|b_v| overstated the true bias by a factor of 3.6–4.3.** Worse, on a "clean metric" containing no confounding at all (projection of the 1,152-dimensional embedding onto a random direction), a set whose quantiles were skewed to 0.75 had a **0% pass rate under the marginal-regression version** (versus 51.7% under the joint-regression version). Using marginal regression would repeat the same failure mode seen before: a sound metric gets rejected. Judgment therefore uses joint regression (partial coefficients); the marginal-regression Σ|b_v| is reported only as a "double-counted upper bound" and is never used to decide pass/fail.

**Derivation of the 0.10 threshold**, in five steps, not inherited from past practice:

1. The judgment to be protected is the set's absolute-value condition, **0.41 standard deviations**, fixed by the PI before the results were seen.
2. Bias is required to stay within one-quarter of this judgment: 0.41 ÷ 4 = 0.1025 → rounded to **0.10**.
3. As a fraction of the judgment threshold: 0.10 / 0.41 = **24.4%**.
4. Compared against irreducible noise—the standard deviation of the set-mean z obtained by repeatedly drawing random 119-gene sets from controls alone (observed 0.0802 / 0.0782 / 0.0806 across three seeds)—0.10 / 0.080 ≈ **1.25×**. This keeps systematic bias at roughly the same scale as the chance variability the design already tolerates.
5. Compared against the standard error of the statistic B itself—the median se_B from 1,500 random 119-gene control sets, **0.0073**—0.10 / 0.0073 ≈ **13.7×**. That is, the 0.10 budget is comfortably larger than measurement precision, and so is actually **measurable**.

**Checking the upper end.** Loosening the budget to half the judgment threshold (0.205) would let bias alone move half the judgment. Tightening it to one-eighth (0.051) falls below the 95th percentile of random control-set noise (observed 0.0742), meaning even unbiased sets would start failing. The adopted value, 0.10, lies between these two extremes.

**Why this matters.** Points 4 and 5 holding simultaneously is the crux: the budget must be ≤0.10 to protect the judgment (points 2–4), and 0.10 is 13.7× measurement precision, so it is actually measurable (point 5). The old check, M5 (discussed below), satisfied the analogue of point 5 by only 2.1×—it demanded a precision finer than the measurement could deliver.

**Limitation.** This method only estimates bias explainable by the five enumerated candidates. It cannot establish "no confounding," only "bias from the enumerated candidates is below 0.10 standard deviations." When candidates are noisy proxies for a true, unobserved confounder, recovery of the true bias in simulation ranged from 87% to 95%.

### S-2: Binomial Test on Sham Sets

**Definition.** Two hundred "sham sets" are constructed from control genes, matched to the evaluated set on properties such as expression frequency, position in the input, input length, and cell count. Each is passed through the same pass/fail rule. Using the empirically measured chance pass rate p for a single gene drawn from controls, a binomial test (k ~ Binomial(200, p)) is applied to the number of sham sets that pass, k, and the set is judged acceptable if the one-sided p-value is ≥ 0.05 (i.e., consistent with chance).

**Why "zero of N sets passing" is not used as the requirement.** An earlier phase of the project required "zero of 200 sham sets pass." But when the expected number of chance passes was 1.12, demanding exactly zero meant **a well-designed, unbiased setup would fail this requirement 67% of the time.** A few chance passes are expected and should be evaluated via a binomial test for consistency with chance, not treated as automatic failure.

**Limitation.** For the K=119 intervention set, the chance pass rate is on the order of 5×10⁻⁶, giving an expected pass count among 200 sham sets of only 0.001. At this scale the test has almost no power. Passing this check does not mean "the sham sets failed, therefore the real result is genuine"; this limitation is stated explicitly in the report.

### S-3: Positive Control on the Goal Axis

**Definition.** The set-average value for 71 marker genes characterizing the target state (the activated interstitial state) must exceed the 99th percentile of the control-set distribution, and its absolute value must be at least 0.53 standard deviations, in all three seeds.

**Why this is necessary.** This check asks whether the "direction toward the target" even exists as a meaningful axis in this embedding space. Without passing it, there is no basis for interpreting the primary metric at all. In this project, this gauge **failed to register in four consecutive experiments** — EXP-20260927-004, EXP-20260928-001, EXP-20260928-002, and EXP-20260929-001, with positive-control gene counts of 1, 24, 0, and 7 respectively — before first passing in EXP-20260930-001 (set average 1.002 / 1.093 / 1.067).

**Why the upper bound of bias is subtracted before judging.** The 71 positive-control genes are strongly skewed in their properties (mean quantile among controls: uniform component 0.217, input length 0.334, movement 0.457, original rank 0.545, relative position 0.668). Applying the S-1 bias budget (0.10) directly to such a skewed set would always fail (observed B = +0.126 / +0.129 / +0.150, upper bound 0.238 / 0.239 / 0.258)—but that reflects the fact that the positive control set is inherently skewed, which is a different matter from whether the gauge itself is working. The check therefore uses the positive-control value minus the upper bound of its own bias estimate. Subtracting the upper bound from the previously observed values (0.976 / 1.057 / 1.029) gives 0.718 / 0.818 / 0.771, comfortably above the required 0.53. **Even under the most conservative estimate of bias, the gauge is registering.**

**Derivation of the detection floor, 0.78.** Obtained by simulating the pre-registered decision rule itself: for a set of size K=71, the effect size detectable with 80% probability, incorporating observed between-gene variability (τ = 2.051, see §4.5), is 0.782 standard deviations.

### S-4: Set-Level Neutrality

**Definition.** The set-average value for the four adopted neutral controls must fall within the 0.05–0.95 quantile range of the random control-set distribution in at least two of the three seeds.

**Adoption rules H-1 through H-3 for neutral controls (fixed numerically before seeing results).**

| Rule | Condition | Threshold |
|---|---|---|
| H-1 | Detection rate in both starting and target states | ≥ 0.20 |
| H-2 | Absolute value of the expression bias (log-ratio of detection rates) between states | ≤ 0.25 |
| H-3 | Coefficient of variation of detection rate across donors | ≤ 0.50 |

**Why rules are necessary.** Housekeeping genes (e.g., ACTB, GAPDH, B2M) have traditionally been used as neutral controls under the assumption that they are expressed at a constant level in every cell and therefore have no reason to approach the target state. **This assumption does not hold in single-nucleus RNA sequencing (snRNA-seq).** Observed detection rates for B2M were 0.469 (starting state) and 0.425 (target state); for GAPDH, 0.197 and 0.133. Because snRNA-seq captures nuclear, not cytoplasmic, transcripts, even canonically "constant" genes can be detected unreliably.

**Result of applying the rules.** Applying H-1 through H-3 mechanically to 11 candidates (4 housekeeping genes and 7 neutral-control-panel genes) yielded **only 4 adopted genes: ACTB, CERS4, ANAPC16, and EEF1D.** GAPDH failed all three conditions; B2M failed because its coefficient of variation, 0.5002, exceeded the 0.50 threshold.

**Disclosure of the near-miss.** B2M's coefficient of variation (0.5002) and that of VAPA, a neutral-control-panel candidate (0.5015), both exceed the 0.50 threshold only marginally. **The threshold was not adjusted in response to this near-miss**—doing so would defeat the purpose of fixing the rule in advance. Instead, a version with the threshold relaxed to 0.55 (adding B2M and VAPA, for 6 genes total) was pre-registered as a sensitivity analysis **before results were seen.**

## 4.4 Dropped Checks and Why They Were Dropped

On 2026-10-02, at the PI's direction, the old checks M1 through M5 were audited. The audit used 60 "clean metrics"—metrics guaranteed by construction to be free of confounding and circularity, obtained simply by projecting the displacement vector onto 60 random directions in the 1,152-dimensional embedding space, unrelated to the actual confound-removal pipeline. Applying the old checks M1–M5 to these 60 clean metrics gave the following pass counts:

| Check | Pass rate (of 60) | Finding |
|---|---|---|
| M1 (calibration of the control-gene rank distribution) | 60 / 60 | Does not trigger regardless of injected defects. Carries no information |
| M2 (rank correlation of residual with original position and input length) | — | Undetectable in principle. Since the confound-removal regression's predictors are exactly the full quadratic in original rank and input length, correlations with these quantities are absorbed by definition. Even artificially injecting a ±0.30 defect left the statistic unchanged to four decimal places |
| M3 (rank of the four individual neutral controls) | 43 / 60 | Four genes is too few; the signal is swamped by noise. No derivation for the threshold was ever recorded |
| M4 (rank correlation of residual with the uniform component) | 49 / 60 | No derivation for the threshold was ever recorded |
| **M5 (rank correlation of residual with movement magnitude)** | **8 / 60 (13%)** | The threshold of 0.10 is only **2.1×** the statistic's own standard error (0.045–0.047)—a level of precision the measurement cannot deliver |
| **All five passing together** | **7 / 60 (12%)** | Even a metric completely free of confounding **fails at least one check 88% of the time** |

**The significance of this 88% figure is substantial.** If the probability of failing on any single run is 88%, the probability of failing on six runs is 0.88⁶ ≈ **46%**. In other words, the six "inconclusive" results **standing as of 2026-10-02** fell well within the range expected even if the underlying metric was perfectly sound. This does not mean the metric was flawed and caused the repeated inconclusive outcomes; it means **the check battery itself was constructed in a way that would reject a sound metric.**

**Why each check failed to function:**

- **M2 (original rank, input length) is undetectable in principle**, because the confound-removal regression's formula is exactly that quadratic. The correction absorbs any such defect by construction, so passing this check demonstrates nothing about the absence of confounding.
- **M1 does not trigger at this project's scale** (roughly 500 control genes). No injected defect caused its value to cross the threshold.
- **M3 is swamped by noise with only four genes (neutral controls)**: even a sound metric produces highly variable individual-gene ranks, and four genes provide insufficient stability for a judgment.
- **M5 and the "original rank" component of M2 were measuring substantially the same thing.** Their correlation was **0.978**, driven by the fact that, among control genes, movement magnitude and original rank are themselves strongly correlated (0.873–0.878). Two checks that appeared independent were in fact carrying one piece of information.
- **The thresholds were not calibrated against each other.** Comparing, for a fixed biased set, the bias each check's threshold would tolerate as a fraction of the judgment line (0.41 standard deviations): M5 tolerated 36%, M4 tolerated 104%, and M2 (relative position) tolerated 151%. M4 and M2 (relative position) let through bias large enough to overturn the judgment, while M5 alone was excessively strict.
- **M5's threshold of 0.10 was only 2.1× the statistic's own standard error (0.045–0.047).** The instability this produces is visible directly: between two experiments whose embeddings did not differ by a single bit (i.e., computed from the identical model run), simply swapping 59% of the control-gene roster shifted M5's value from −0.078 to −0.109, **flipping the pass/fail decision.** A check that reverses its verdict purely because of how controls happen to be drawn—with the measured object itself unchanged—cannot be relied on.

These checks are no longer used to decide pass/fail, but their values are retained and reported descriptively in full. "Dropping" a check means removing it from the pass/fail decision, not deleting it from the record.

## 4.5 Designing for Statistical Power

**The detection floor is the smallest true effect size (in units of control-gene standard deviation) that a given decision rule can detect with 80% probability.** It is obtained by directly simulating the pre-registered decision rule itself—random draws from controls, the 99th-percentile condition, the absolute-value condition, and the requirement that all three seeds pass—rather than by any closed-form approximation.

**Why a crude formula (4.5/√K, where K is the set size) is not used in reporting.** Such a formula is useful as a rough planning-stage guide, but it diverges systematically from direct simulation. For a set of size K=11 (the predecessor to the current intervention set), the formula gives 1.36 standard deviations, while direct simulation of the actual decision rule gave **2.00 standard deviations**. For a K=39 set (a subset of positive controls from a separate experiment), the formula gives 0.69, while simulation gave **1.03**. The formula underestimates systematically because it omits two quantities observed in the data: between-gene variability in effect size (observed τ = 2.051) and the correlation between seeds (observed around 0.944). Because of this, the detection floor is always computed by simulating the actual decision rule; the simple formula is confined to early planning and never used in reporting.

**Observed relationship between set size K and detection floor:**

| K | Detection floor (SD) | Source |
|---|---|---|
| 11 | 2.00 | `EXP-20260930-001/work/power_004.json` (set T1) |
| 18 | 1.56 | same source (set T1_T2) |
| 39 | 1.03 | same source (positive_control39); note the same-named set in `EXP-20261001-001/work/power_005.json`, computed against a different control pool, gives 1.06 |
| 71 | 0.78 | `EXP-20261001-001/work/power_005.json` (positive_control71) |
| 119 | 0.604 | same source (INT, the primary judgment) |

Larger K produces a smaller detection floor, making the same true effect easier to detect. The primary judgment's detection floor of 0.604 (K=119) is about **3.3 times smaller** than the detection floor for the K=11 set used in earlier experiments (2.00)—this difference in statistical power is the single largest design change separating this experiment from its predecessors.

**A single gene cannot be reliably assessed.** The detection floor for a single gene (effectively K=1) reaches 2.73 standard deviations, far exceeding even the largest effects observed in this project (roughly 1.0 standard deviation for the positive-control set). This is why the unit of judgment was placed at the level of the gene set rather than the individual gene.

**The pre-feasibility check (`scripts/analysis/axis_feasibility.py`).** Before any perturbation is run, this check determines whether the target axis can be measured at all under the proposed design. The single strongest condition is **a positive-control gene count K ≥ 20**; designs falling below this are not attempted. Empirically, a K=7 positive-control set has a detection floor of 1.71 standard deviations, and **no positive-control set in this project has ever reached that level.** No amount of careful check design can rescue a set too small to detect the effect it is meant to detect. The check returns one of three verdicts—pass, hold-for-review, or fail. A fail blocks the start of the experiment until the design is reworked; a hold-for-review requires the identified weaknesses to be documented in the plan and referred to the PI for judgment.

## 4.6 Discipline of Pre-Registration

Pre-registration means fixing, in a document (`plan.md`, `config.yaml`), the decision criteria, the target genes, and the cell-selection method *before* results are seen. The fixed content is recorded as a sha256 hash (`work/SHA256SUMS-preregistered.txt`) and verified from within the experiment directory. This verification mechanically confirms that not a single byte has changed; it does not rely on memory or self-report.

**Any change made after fixing the pre-registration requires all three of the following:**

1. Record the date, content, reason, and **whether the change was made before or after seeing results** in a change log.
2. Update the pre-registration checklist.
3. Re-compute the sha256 hashes and re-fix until every file matches again.

**Criteria are never changed after seeing results.** Even when a check turns out to be inappropriate, the response is not to draw a conclusion from the flawed check's output, but to redesign the check, **re-fix the redesigned check before seeing results again**, and only then re-run the evaluation. When S-1 through S-4 were rebuilt in this project, the existing perturbation data (embeddings already obtained in EXP-20261001-001) were not re-generated at all; only the definitions of the checks were re-fixed in advance, and the data were then re-evaluated (EXP-20261002-001).

**Re-evaluations performed after seeing results are always disclosed as such.** The report states plainly, "this is not pre-registered analysis performed before running the perturbation, but a re-evaluation of existing data." Throughout this project's reporting, every number is marked as having been fixed before results were seen, or as arising from a post-hoc re-evaluation.
