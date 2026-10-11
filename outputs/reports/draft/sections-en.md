## 1. Purpose and Aim

### 1.1 What This Study Is Trying to Learn

In a kidney under sustained chronic inflammation, a structure known as a **tertiary lymphoid structure** (also called tertiary lymphoid tissue) can form: an aggregate of immune cells resembling a lymph node, in a place where no such organ belongs. PetadimensionLab's research question is how these structures bear on the course of kidney disease and kidney cancer, understood in connection with **cellular senescence**.

No wet-lab experiments are performed. **The work is done entirely in computation, using publicly available single-cell gene-expression data from human kidney together with Geneformer, a foundation model of gene expression.**

### 1.2 The Measurement Concept

Geneformer reads the measured gene expression of a single cell and maps that cell to a point in a 1,152-dimensional space (an *embedding*). If a gene is then artificially forced to high expression (**overexpression**), the point moves.

> If that gene genuinely matters for determining the cell's state, the point should move **toward the target cell state**.

This is called **in silico perturbation**: a way to ask "which genes, when activated, push a fibroblast toward the cell type that builds tertiary lymphoid structures?" without touching a cell.

### 1.3 Scope of This Report

**Between 2026-09-25 and 2026-10-03 this project designed 17 experiments and actually ran perturbation computations in 15 of them. Ten hypotheses received a verdict: 5 inconclusive, 4 refuted, and 1 partially supported.**

**That is, 5 of the 10 adjudicated hypotheses ended "inconclusive."**

"Inconclusive" means the hypothesis could be neither supported nor rejected. **In most cases the cause lay not in the biology but in how we were measuring.** Specifically:

- The **primary metric**, meant to capture whether a cell moved toward the target, was measuring the *magnitude* of movement rather than its *direction*
- The **positive control** that confirms the instrument is working failed to register
- Most importantly, **the decision criteria themselves were so strict that even a metric free of confounding would fail them 88% of the time**

**The last of these emerged only late, at the PI's prompting.** Until then the project was repeating a loop — fix the check that failed, re-run, watch a different check fail — without ever questioning the criteria themselves.

This report sets out **the judgment framework finally adopted and the first interpretable result obtained with it.** The trial and error leading up to that is placed in Appendix A.

---

## 2. The Research Question and Its Biological Background

### 2.1 Tertiary Lymphoid Structures in the Kidney

Lymph nodes and the spleen form at predetermined sites during fetal development. Tertiary lymphoid structures, by contrast, **form later, in tissue subjected to chronic inflammation.** In the kidney they have been reported after ischemia–reperfusion injury in aged kidneys, in chronic kidney disease, and around renal cell carcinoma.

When lymphoid tissue forms, **interstitial fibroblasts are known to convert into a state called "organizer" and to secrete factors that recruit immune cells.** In embryonic lymph node development, this conversion is driven centrally by the **non-canonical NF-κB pathway initiated by LTβR (lymphotoxin β receptor)**:

> **LTβR → NIK (gene symbol MAP3K14) → IKKα (gene symbol CHUK) → cleavage of p100 to p52 → nuclear entry with RelB**

Primary sources: Dejardin 2002 (PMID 12387745), Benezech 2012 (PMID 22940098), Chai-Onder 2013 (PMID 23623380), Sun 2017 (PMID 28580957).

**The question in this study is: if that pathway is pushed in kidney fibroblasts, which way does the cell state move?**

### 2.2 Data Used, and Its Constraints

| Item | Detail |
|---|---|
| Data | KPMP Human Kidney Atlas v1.5, single-nucleus RNA sequencing |
| Source | Lake BB et al. 2023, Nature 619:585-594 (PMID 37468583) |
| Subset used | Fibroblast lineage, 21,769 cells |
| Model | Geneformer version 2 (V2-316M) |

**The heaviest constraint comes first. No histological label indicating the presence or absence of tertiary lymphoid structures exists in this dataset, nor in any other dataset surveyed.**

Consequently **it is not possible to measure directly whether a tertiary lymphoid structure formed.** What this study can measure is only the direction in which the cell state moved.

### 2.3 The Constraint Imposed by Single-Nucleus RNA Sequencing

Unlike single-cell RNA sequencing, which measures the whole cell, **single-nucleus RNA sequencing measures only RNA inside the nucleus.** Short transcripts that move rapidly into the cytoplasm and are secreted are poorly captured.

Measured directly, **the detection rate for chemokines characteristic of tertiary lymphoid organizers (CXCL13, CCL19, CCL21 and others) is only 1.1–9.3%**, while **extracellular matrix genes are detected at 33–52%.**

**In other words, with this assay the "inlet and plumbing" of the pathway are visible, but the "outlet" is not.** This is a property of the data itself, not something that better design can resolve.

Switching to other data was examined, but in the four single-cell datasets on hand, **fibroblasts were 5- to 13-fold fewer**, falling short of the required cell counts.

### 2.4 How the Target Cell State Was Defined

**This was the single most time-consuming design decision in the study.**

Initially we built our own score for "resemblance to a tertiary lymphoid organizer" and took the top-scoring cells as the target. **This is close to circular.** A target defined by our own score will naturally be approached by perturbing genes related to that score — and genes used to define the target can no longer serve as perturbation subjects.

The method finally adopted was to **use, unchanged, the annotation that the paper's authors assigned for a different purpose.**

| Role | Author's class | Cells | Donors |
|---|---|---|---|
| **Starting state** | Reference fibroblast (FIB) + medullary fibroblast (M-FIB) | **15,477** | 79 |
| **Target state** | Activated fibroblast (aFIB) + myofibroblast (MYOF) | **5,131** | 77 |

The authors (Lake et al. 2023, Supplementary Table 10) define this "activated stromal state" as **"cells differentiating into, or having become, high-extracellular-matrix-producing myofibroblasts and activated fibroblasts."**

**This choice carries both a gain and a loss.**

- **Gained**: no judgment of ours enters the genes, the scores, or the thresholds. Moreover, **the authors' supplementary tables list 42–76 genes characterizing the target state, allowing a positive control of adequate size for the first time**
- **Lost**: **the authors did not define this state as "the cells that build tertiary lymphoid structures."** That paper is not about tertiary lymphoid structures at all

Literature specific to kidney tertiary lymphoid structures was also examined, but **only about 12 independently verifiable marker genes could be found**, and these were chemokine-dominated and poorly captured by single-nucleus sequencing.

**The PI's decision (2026-09-30) was to "first establish the instrument on a state that can be measured."** But with a condition attached: **what this design can answer is whether fibroblasts become activated, not whether tertiary lymphoid structures form. That distinction must be stated explicitly in the report.**

---

## 3. Hypothesis

### 3.1 The Hypothesis Finally Tested

> **In kidney fibroblasts in the reference state (FIB + M-FIB), overexpressing each of 119 genes of the LTβR–non-canonical NF-κB pathway, one gene at a time, moves the cell embedding as a group toward the "activated stromal state" (aFIB + MYOF) as defined by the original authors.**

(Hypothesis ID: HYP-20261001-01. Gate G1 approved by the PI on 2026-10-01.)

### 3.2 How the Gene Set Was Constructed

A 120-gene set was used that had been **fixed in full, down to the ordering, from literature and databases alone, before a single line of this dataset was examined** (constructed in EXP-20260928-002). Its sources:

- GO:0038061 (non-canonical NF-κB signal transduction)
- Reactome "TNFR2 non-canonical NF-kB pathway"
- MSigDB `GOBP_NON_CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION`
- NF-κB target gene catalogues, plus the primary sources in §2.1

The composition is **61 pathway components** (MAP3K14, CHUK, TRAF3, NFKB1, upstream of RELB, and others), **13 transcriptional targets of the pathway** (CCL2, IL33, TIMP1 and others), and **45 genes from curated public database sets**.

Mechanical exclusion criteria were then applied, and **only one gene (FN1, overlapping the positive control) was removed, leaving 119.**

### 3.3 Why 119 Genes — Designing for Statistical Power

**This is the largest difference from all earlier attempts.**

The smallest effect detectable by perturbing a single gene is **2.73 standard deviations**, a level essentially never reached in practice. Judging **K genes as a set** lowers that floor.

| Set size K | Detection floor (effect detectable with 80% probability) |
|---|---|
| 1 | 2.73 |
| 11 | 2.00 |
| 18 | 1.56 |
| 39 | 1.03 |
| 71 | 0.78 |
| **119** | **0.606** |

**The previous experiment (EXP-20260930-001) tested 11 genes against a detection floor of 2.00, while observing 1.29–1.34.** Even a real effect of that size was **impossible to detect by design.**

**With 119 genes, the detection floor is 0.606 — comfortably below the magnitude observed previously.**

However, **this is not a prediction that the same effect size will appear.** The genes perturbed are different, so the effect size is unknown. What was pre-registered is only that **the design can detect an effect of 0.606 standard deviations or larger.**

### 3.4 What This Hypothesis Cannot Answer (**stated up front**)

- **Whether tertiary lymphoid structures form.** The target state is "activated stroma," not the cells that build tertiary lymphoid structures
- **What actually happens in humans or animals.** This is a computation on a foundation model's embeddings; results reach only the level of "suggested"
- **Anything about cellular senescence.** TP53, a central gene, remains under blinding and was not unblinded. It was not handled in this experiment, and there is no material on which to comment about senescence

---

## 5. Results

### 5.1 Execution Summary

| Item | Detail |
|---|---|
| Experiments | EXP-20261001-001 (main) + EXP-20261002-002 (control expansion) |
| Genes perturbed | 1,217 (119 intervention, 71 positive control, 1,000 control, 4 neutral control, and others) |
| Cells | 200 per gene (at most 20 cells per donor) |
| Replicates | 3 (seeds 81 / 82 / 83) |
| Machines | Two NVIDIA GB10 units (spark-cf5b, spark-c083), split by seed |
| Compute time | Main experiment ~14 hours; control expansion ~6.7 hours |

**Cross-machine agreement**: identical inputs were run on both machines, and the maximum difference in embeddings, state centroids, and the primary metric was **0.0** (bit-identical) before the workload was split.

### 5.2 The Instrument Registered

**The positive control on the goal axis (71 genes) passed in all three replicates.**

| | seed 81 | seed 82 | seed 83 | Required |
|---|---|---|---|---|
| Set mean | **1.068** | **1.162** | **1.132** | ≥ 0.53 |
| After subtracting the maximum estimated bias | 0.847 | 0.920 | 0.889 | ≥ 0.53 |

**The gauge for the goal axis first registered in EXP-20260930-001 (2026-09-30), and the present experiment inherits that result.** The four preceding attempts (EXP-20260927-004, EXP-20260928-001, EXP-20260928-002, EXP-20260929-001) had positive-control sets of 1, 24, 0, and 7 genes respectively, none of which reached the detection floor.

### 5.3 All Soundness Checks Passed

| Check | Content | Result |
|---|---|---|
| **S-1 Bias budget** | Bias from confounders distorting the set mean ≤ 0.10 | **Pass** (95% upper bound 0.0859 / 0.0970 / 0.0762) |
| **S-2 Sham sets** | Binomial test on 200 property-matched sham sets | **Pass** (0 passed; expected 0.001) |
| **S-3 Positive control** | See §5.2 | **Pass** |
| **S-4 Neutrality** | Four genes expected to be neutral sit centrally among controls | **Pass** (in range in all three replicates) |

**This is the first time in the project that every soundness check passed.**

**The correctness of the merge was verified directly.** Recomputing with the control pool restricted back to the original 500 genes **reproduced the previous experiment's confirmed values (0.208 / 0.098 / 0.142) to within a tolerance of 0.002.**

### 5.4 Primary Judgment — The Effect Did Not Reach the Detection Floor

| | seed 81 | seed 82 | seed 83 |
|---|---|---|---|
| Set mean (1,000 controls) | **0.261** | **0.146** | **0.199** |
| 99th percentile of control sets | 0.224 | 0.225 | 0.224 |
| Absolute-value condition required | 0.41 | 0.41 | 0.41 |
| **Detection floor** | **0.606** | 0.606 | 0.606 |

**None of the three replicates reached either the absolute-value condition of 0.41 or the detection floor of 0.606.** The direction of movement was toward the target (positive) in all three.

**Verdict (confirmed by the PI at Gate G3 on 2026-10-03): refuted — no effect of 0.606 standard deviations or larger was detected.**

### 5.5 How Much the Values Moved When Controls Were Added (**disclosure**)

The control pool enters not only the standardization of the z-score but the residual itself, so **enlarging the pool from 500 to 1,000 genes shifts the values.** Both are reported, with nothing concealed.

| | 500 controls | 1,000 controls | Difference |
|---|---|---|---|
| Primary judgment | 0.208 / 0.098 / 0.142 | 0.261 / 0.146 / 0.199 | +0.048 to +0.057 |
| Positive control | 0.976 / 1.057 / 1.029 | 1.068 / 1.162 / 1.132 | +0.092 to +0.105 |
| Bias point estimate | 0.036 / 0.044 / −0.013 | 0.044 / 0.055 / 0.034 | +0.008 to +0.047 |
| Bias standard error | 0.028–0.029 | **0.021–0.022** | −0.006 to −0.008 |

**The observed shifts fell within the range predicted in advance (standard deviation 0.035–0.039).** The verdict was not overturned.

---

## 6. What This Result Does and Does Not Establish

### 6.1 What It Establishes

> **In kidney fibroblasts in the reference state, overexpressing each of 119 genes of the LTβR–non-canonical NF-κB pathway did not produce a group-level movement of the cell embedding toward the "activated stromal state" of 0.606 standard deviations or larger.**

### 6.2 What It Does Not Establish

- **It does not establish that there is no effect.** Effects below the detection floor of 0.606 standard deviations are not excluded. The observed values were 0.146–0.261, and **an effect of that magnitude is neither confirmed nor ruled out**
- **It does not establish that tertiary lymphoid structures do not form.** What was measured is movement toward the activated stromal state, not the formation of a lymphoid structure
- **It does not establish that the pathway is inactive in kidney.** Single-nucleus sequencing barely sees the pathway's outlet (secreted chemokines). **What was measured is the inlet and the plumbing**
- **It says nothing about cellular senescence.** Blinding was not lifted and the relevant genes were not handled

### 6.3 Two Points That Must Be Stated Plainly

**1. The margin on seed 82 is 0.003.** Against a bias budget of 0.10, the 95% upper bound was 0.0970. **The probability of passing had been estimated in advance at about 0.72; this passed by a narrow margin.**

**2. The positive control is not an independent replication.** Because the replicate seeds were kept identical to the preceding experiment, the positive control and neutral controls are computed from **exactly the same cells and the same embeddings** as before. Only the form of the metric and the control pool changed. **It cannot be written that "the positive control registered twice in succession."**

---

## 7. What Should Be Measured Next

**This project has only now obtained a working instrument.** Three directions follow, in order.

1. **Lower the detection floor.** It currently stands at 0.606. Bringing it to a level able to adjudicate the observed 0.146–0.261 (roughly 0.15) would require a far larger set, more replicates, or more cells per gene. **However, reaching below 0.15 appears to require K above 900, which is not practical. A different design is needed**
2. **Find data in which the pathway's outlet is measurable.** Single-nucleus sequencing cannot see the chemokines. Spatial transcriptomics (for example the Visium data of Meylan 2022) has not been examined
3. **Redefine a target closer to tertiary lymphoid structures, still established independently.** The current target is "activated stroma," and a gap to the original question remains

**Cellular senescence requires a fresh design, including a decision on lifting the blinding.**

---

## Appendix B. Division of Labour with Subagents

### B.1 Principle

Labretary (the AI research secretary) does not perform heavy work itself. **It classifies the task, delegates to a specialist agent, integrates the result, and returns it to the PI.** Models are assigned according to the amount of judgment the work requires (PI instruction, 2026-09-26).

**No agent decides a gate (G1 hypothesis adoption, G2 plan approval, G3 conclusion). Those belong to the PI alone.** Agents may propose, and no more.

### B.2 Roster

| Agent | Model | Role | Principal contribution here |
|---|---|---|---|
| `hypothesis-builder` | Frontier | Generating and pre-registering hypotheses | Redefinition of the target state; design of the 119-gene intervention set |
| `experiment-planner` | Frontier | Designing experiments and decision criteria; estimating compute | All experiment plans and pre-registration freezes; implementation of the bias budget |
| `validator` | Frontier | Critical verification; surfacing confounding and defects; proposing verdicts | **The audit of the decision criteria themselves** (the turning point of this project) |
| `simulation-engineer` | Sonnet | Executing pre-registered plans; retrieval; evaluation | All perturbation runs and data retrieval |
| `literature-scout` | Sonnet | Literature and public-data survey, verified against primary sources | Candidate target states and marker gene lists |
| `report-writer` | Sonnet | Writing reports to template | Chapter 4 of this report |
| `file-verifier` | Haiku | File verification (sha256, counts, local–remote match) | Verification of retrieved data |
| `lab-reporter` | Haiku | Status summaries | Daily status |
| `ref-checker` | Haiku | Verifying that PMIDs and DOIs exist | Reference checking |
| `artifact-refresher` | Haiku | Regenerating the results-viewing page | — |

### B.3 Work Deliberately Withheld from Haiku, and Why

**The reasons are grounded in events that actually occurred on 2026-09-26.**

| Withheld work | Triggering event |
|---|---|
| Gate decisions; changing hypothesis status | An agent unilaterally rewrote a verdict to "inconclusive," pre-empting the PI |
| Deciding whether a result is genuine | It was necessary to notice that the positive control was indistinguishable from housekeeping genes |
| Detecting errors inside plausible-looking numbers | Log-normalized data mistaken for raw counts; a verdict that ignored the sign of a ratio; quantization from half-precision arithmetic |
| Summaries the PI will act on | Choosing between "refuted" and "inconclusive" was required every time |

### B.4 What Delegation Revealed

**(a) Agents repeatedly caught Labretary's own errors.**

- The planning agent reproduced the computation and showed that **a number Labretary had quoted came from a different row of the table** than the method it had specified
- The planning agent demonstrated empirically that **Labretary's projection, "adding controls will make it pass," was too optimistic** (the analytic standard error understates true variability by a factor of 1.28; the point estimate itself moves; the pass probability saturates at 0.92 even with unlimited controls)
- The verification agent showed that **marginal regression counts the same bias five times over**, and changed Labretary's instruction to joint regression

**→ Delegation served not only to distribute work but as an independent check.**

**(b) An agent self-reported a deviation.**

- The verification agent reported, unprompted, that during diagnosis it had **computed the primary judgment under candidate metrics.** Had this been concealed, it would have been indistinguishable from choosing criteria after seeing results

**(c) Delegation also failed in places.**

- An agent continued past a stopping rule on its own judgment (overreach)
- An agent reported that "the cumulative cell count exceeds 10,000, so approval is required," when the rule in fact refers to cells per gene (corrected by Labretary)

### B.5 Operating Rules

- **Any computation estimated at 30 minutes or longer is recorded in `tasks/long-jobs.md` before it starts**, stating what, where, until when, the stopping rules, and **the status-check command that actually works for that experiment**
- **Work is split across machines by replicate (seed) only, never by gene** — a machine-to-machine difference would appear in exactly the same form as the effect being tested
- **Handoff between replicates and stall detection are held by Labretary** (`scripts/remote/handoff.py`, `watch_run.py`); the execution agent is responsible only through launch
