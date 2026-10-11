# Technical Changes — Shared Scripts and Execution Environment

Date: 2026-10-03
Period covered: 2026-09-25 through 2026-10-03
Repository: `/Users/petadimensionlab/workspace/research/ai-labretary`

---

## 1. Scope

### In scope

- Changes to shared scripts under `scripts/`, each with symptom, cause, fix, verification, and date
- Design intent of newly added shared modules, and the contract callers must honor
- Known pitfalls that have **not** been fixed in code, with the procedural workaround
- Environment-level failures (macOS, shell, GPU, numeric precision)
- Statistical and numerical implementation fixes
- Operating procedures that have been frozen as rules
- Open problems

### Out of scope

- **No research judgments.** Hypothesis acceptance, experiment pass/fail, metric interpretation, and biological conclusions are not covered here. Those live in `hypotheses/tracker.md`, each experiment's `plan.md` / `validation.md`, and the research report.
- Debate over whether a threshold or decision rule is *scientifically* sound. This document covers only how such rules are implemented, where they broke, and how they were repaired.
- Per-experiment results. Numbers appearing here are limited to **measurements that motivated or validated a code change**.

### On line numbers

Every line number below was read from the actual file on 2026-10-03. `scripts/` is not under version control, so later edits will shift these.

---

## 2. Changes to shared scripts

### 2.1 `scripts/geneformer/run_isp.py` (1,116 lines)

#### (1) Gate G2 threshold raised from 2 GPU-hours to 24

| | |
|---|---|
| **Symptom** | Under the old 2-hour threshold, nearly every production run required approval, blocking work that could have proceeded automatically. |
| **Cause** | The threshold was set well below the real scale of a single job (several to a dozen-plus hours). Experiments grew larger than the original design assumed. |
| **Fix** | `run_isp.py:415`, with the directive date retained in the comment:<br>`out["needs_G2_approval"] = hours > 24 or n_cells > 10000 or g == "all"  # GPU 24 h（PI 2026-09-27）` |
| **Verification** | Run `run_isp.py --dry-run --config <config>` and read `needs_G2_approval` from the output JSON. No computation is performed. |
| **Date** | 2026-09-27 (PI directive) |

**What `n_cells` means — cells *per gene*, not the cumulative total.**

`run_isp.py:387`:

```python
n_cells = sel.get("max_cells_per_gene") or cfg["data"].get("max_ncells") or 10000
```

`sel` is the `cell_selection` section, and `max_cells_per_gene` caps the starting cells drawn **for one perturbation target gene**. So `n_cells > 10000` does *not* ask whether the cumulative cell count exceeds 10,000.

- Concrete case: EXP-20260930-001 totalled 360,000 cell×gene forwards, yet `max_cells_per_gene` was at or below 10,000, so `needs_G2_approval` was `false`. A delegate misread this on 2026-10-01 and raised an unnecessary approval request.
- **Rule**: never compute G2 eligibility by hand. Read `needs_G2_approval` from `--dry-run`.

**The only mechanical enforcement** is `run_isp.py:988-989`:

```python
if os.environ.get("LABRETARY_MODE") == "night" and est["needs_G2_approval"] and not args.approved:
    sys.exit("夜間モード: G2 承認が必要なジョブです（--approved が必要）")
```

Daytime runs are not blocked. Outside night mode, the only thing stopping an unapproved production run is the operating rule in `CLAUDE.md` §4.

---

#### (2) `--smoke` did not lower `min_cells_per_gene`, so every gene was skipped and the run "succeeded"

| | |
|---|---|
| **Symptom** | A smoke run exited 0 and logged completion, yet **not a single perturbation had executed**. In EXP-20260928-001 all 41 genes were skipped and downstream evaluation failed with "stats CSV not found." |
| **Cause** | The old `--smoke` lowered only `max_cells_per_gene` to 50, leaving `min_cells_per_gene` (the floor on cells actually drawn) at its production value, e.g. 120. Every gene then hit the skip branch in `select_cells()` at `run_isp.py:668-669`: "drew 50 < min_cells_per_gene → skip". A skip emits a warning and does not change the exit code. |
| **Fix** | `run_isp.py:1016-1022` now lowers the floor to at most 5 under `--smoke` **and logs that it did so**:<br>`sel_s["min_cells_per_gene"] = min(min_before, 5)` |
| **Verification** | After a smoke run, confirm `<results_dir>/**/emb/{orig,pert}.npy` exist with `ls`. Also count `WARN ... → skip` lines in the log to confirm how many genes actually ran. |
| **Date** | 2026-09-28 |

**Related pitfall (not fixed in code).** `run_isp.py:978`:

```python
base = cfg.get("output", {}).get("results_dir") or ("results-smoke" if args.smoke else "results")
```

If a config sets `output.results_dir` explicitly, **`--smoke` writes into the same directory as production**. Combined with the resume logic (skip a gene if its stats already exist), a 50-cell smoke output can silently contaminate a production run. Always give smoke runs their own `results_dir`.

---

#### (3) `--smoke` forces 2 controls and 3 targets, so evaluation code cannot be exercised on its output

| | |
|---|---|
| **Symptom** | The requirement "run the evaluation code end to end before production" cannot be satisfied using `--smoke` output. Evaluation either fails to run or runs without exercising the confound-removal path. |
| **Cause** | Regardless of what the config says, `--smoke` forces two reductions:<br>• `run_isp.py:1012` — `cfg["controls"]["negative_random_n"] = min(2, cfg.get("controls", {}).get("negative_random_n", 0))`, capping random negative controls at **2**<br>• `run_isp.py:1063` — `targets = g[:3] if args.smoke else g`, cutting targets to the **first 3**<br>**Two controls cannot support leave-self-out confound removal.** `movement_correction.gene_stats()` fits the regression on the control pool with the target gene's own contribution subtracted out; with two controls, one remains. |
| **Fix** | No code change. **The procedure changed.** Smoke runs now use a dedicated config file (`work/config-smoke.yaml`, with explicit floor, cap, and output directory) run **without** the `--smoke` flag. The smoke set must contain at least one gene of every role used in the decision (primary target, positive control, negative control, uniform-component reference), and roughly **30 control genes**. |
| **Verification** | Run the production evaluation script to completion on the smoke output and confirm the decision JSON is produced without exception. For EXP-20261002-002, evaluation over 723 genes / 506 controls / 1 seed completed in 75 seconds. |
| **Date** | 2026-10-01 (discovered and codified) |

---

#### (4) New: `state_emb_dataset_with_rowid()` and `rows_file_path()`

**`state_emb_dataset_with_rowid()` (`run_isp.py:630-646`)**

| | |
|---|---|
| **Symptom** | `precheck_ci.py` tried to join `cell_embs.csv` to the input dataset by row position and halted. In EXP-20260927-002 only 39 of 2,311 cells happened to line up. |
| **Cause** | `EmbExtractor.extract_embs` reorders rows internally via `downsample_and_sort` (stable sort, descending input length) for batching efficiency. The `(state, donor)` pair cannot disambiguate multiple cells from the same donor. |
| **Fix** | Materialize a dataset carrying the pre-sort row index as a `row_id` column and include `row_id` in `emb_label` (called from `precheck()` at `run_isp.py:787`; the `EmbExtractor` call passes `emb_label=[key, dkey, "row_id"]`). `cell_embs.csv` then carries `row_id` and the join is by ID. **No additional forward passes** — the dataset gains a column, nothing more. |
| **Verification** | Confirm the joined row count equals the input cell count and that `row_id` has no duplicates. |
| **Date** | 2026-09-27 |

**`rows_file_path()` (`run_isp.py:571-574`)**

| | |
|---|---|
| **Purpose** | Let an experiment's analysis code pre-register the exact cells used for state centroids as a **JSON list of row indices** (e.g. depth-matched centroids). |
| **Implementation** | `state_embs.rows_file` may contain `{seed}`; `rows_file_path(cfg, seed)` resolves it via `str(rf).format(seed=seed)` and returns `None` when unset. |
| **Consumers** | `emb_name()` embeds the first 8 hex of the rows file's sha256 into the dataset name, so **a dataset built under different selection rules is never silently reused**. `state_emb_dataset()` (`run_isp.py:588-627`) selects exactly those rows. |
| **Built-in checks** | Three hard stops, all `SystemExit`: rows missing for a state (`run_isp.py:606`), a listed row whose state label does not match (`609`), duplicate row indices (`612`). |
| **Caveat** | `rows_file` is **not** transferred by `isp_remote.sh push` — see §2.5. |
| **Date** | 2026-09-27 |

---

### 2.2 `scripts/analysis/movement_correction.py` (new, 260 lines)

#### (1) Movement-magnitude correction

**What it does.** Adds the log of perturbation displacement, `ld = log2(‖Δ‖)`, to the quadratic confound-removal regression.

| Variant | Regressors | Design-matrix columns |
|---|---|---|
| Previous | full quadratic in `lp` (log original rank) and `ll` (log input length) | 6 |
| This module | full quadratic in `lp`, `ll`, **`ld`** | 10 |

`design()` (`movement_correction.py:66-87`) returns 6 columns when `ld is None` and 10 otherwise. **The leading columns are identical across variants**, so coefficients can be compared directly. A single `use_dnorm` argument switches between them (`_cols()`: `movement_correction.py:89-99`; `gene_stats()`: `104-150`).

With `use_dnorm=True`, `_cols()` raises on two conditions: `KeyError` if the `dnorm` column is absent, and `ValueError` if any cell has displacement ≤ 0 (no logarithm).

**Measurement that motivated this** (EXP-20260929-001, validation.md §2.5(a)): median ‖Δ‖ for housekeeping genes sat at the 99.8th percentile of 500 controls for ACTB, 98.6th for GAPDH, 97.8th for PPIA, and the primary metric retained a +0.171 rank correlation with displacement.

**Important limitation, stated as fact.** The correction is not universally beneficial. On EXP-20260930-001 the uncorrected metric was already within range and the correction *degraded* two diagnostic statistics. **Which variant is used for the decision must be fixed in the preregistration before results are seen.** The module computes both side by side by design.

| | |
|---|---|
| **Verification** | `.venv/bin/python scripts/analysis/movement_correction.py --self-test` (`movement_correction.py:183-252`). On synthetic data it checks (a) the fast normal-equations path matches a naive reference to < 1e-9, (b) adding displacement cuts the residual confound by more than half, (c) on confound-free data the two variants agree (correlation > 0.95). |
| **Date** | 2026-09-30 (implementation); 2026-10-01 (`REQUIRED_COLUMNS`) |

---

#### (2) The dropped `relpos_med` column that halted evaluation with `KeyError`

| | |
|---|---|
| **Symptom** | **After 8.6 hours of perturbation compute had finished**, evaluation died with `KeyError: 'relpos_med'`. |
| **Cause** | When `gene_stats()` was extracted from the experiment-local `evaluate_003.py` into the shared `movement_correction.py`, one diagnostic column — `relpos_med` — was not carried over. **The module's own self-test passed**, because it never checked the columns the caller expected. |
| **Fix** | Two parts.<br>(a) Column restored at `movement_correction.py:146`, with an inline note that it is a diagnostic aggregate unrelated to any threshold or gene set.<br>(b) `REQUIRED_COLUMNS` declared at `movement_correction.py:153-154` and enforced at the top of the self-test (`208-211`):<br>```REQUIRED_COLUMNS = ("gene", "col", "D_raw", "q", "n_elig", "n_cells", "lp_med", "ll_med", "relpos_med")``` |
| **Verification** | Self-test item (0), "required columns present", must pass. Immediately after any extraction, diff the new module's output columns against the original implementation mechanically. |
| **Date** | 2026-10-01 |

**Applied elsewhere.** The same guard was added to `bias_budget.py:57` (`REQUIRED_COLUMNS = ("gene", "z")`), enforced at the top of `bias_budget()` (`bias_budget.py:139-141`).

---

### 2.3 `scripts/analysis/axis_feasibility.py` (new, 558 lines)

#### (1) Pre-commitment check: three-valued verdict and exit codes

**Purpose.** Decide, **without running a single perturbation**, whether an experiment's goal axis can be measured at all. Four experiments (EXP-20260927-004, EXP-20260928-001, EXP-20260928-002, EXP-20260929-001) stalled on the same cause — no working scale on the goal axis — and EXP-20260929-001 burned 11.3 machine-hours before that became apparent. Every quantity needed for the verdict is computable without perturbation.

**Quantities measured:**

| # | Quantity | Input |
|---|---|---|
| 1 | Positive-control gene count `K_positive_control` | gene list only |
| 2 | Marker detection rate `marker_detection_median` | expression data |
| 3 | Marker-set coherence `marker_corr_median` | expression data |
| 4 | Per-marker skew `skew_*` | expression data |
| 5 | Axis length ratio `axis_length_ratio` | state centroids |

**Three-valued verdict** (`judge()`: `axis_feasibility.py:391-480`; decision at `472`):

```python
verdict = "不合格" if blocks else ("判定保留" if holds else "合格")
```

i.e. *fail* if any blocking condition fires, *hold* if any weakness fires, otherwise *pass*.

**Exit codes** (`axis_feasibility.py:543-545`):

```python
if a.no_exit_code:
    return 0
return {"合格": 0, "判定保留": 1, "不合格": 2}[res["verdict"]]
```

Anything that is not a verdict — a config error, a bad input, an unhandled exception — exits **3** (`axis_feasibility.py:548-558`). Callers can therefore distinguish pass (0), hold (1), fail (2), and *check did not run* (3). **Do not collapse 3 into 2**: it means the check was not performed, not that the design failed.

**Thresholds** (`DEFAULT_THRESHOLDS`, `axis_feasibility.py:90-105`) were calibrated against the four prior experiments and can be overridden from the config's `thresholds:` section (`axis_feasibility.py:492`). The strongest blocking condition is `block_K_min = 20` — fewer than 20 usable positive-control genes means do not start. `block_axis_ratio_min = 0.35` is only a backstop: the axis length ratio did not separate the four experiments (0.442–0.637 overlapping, with the one partially successful experiment at 0.452), so it is never used as a standalone stop.

| | |
|---|---|
| **Verification** | `.venv/bin/python scripts/analysis/axis_feasibility.py --config <config yaml> [--out <json>]`. All inputs are **read-only**; the only write is the `--out` file. Include the output JSON in the preregistration hash. |
| **Date** | 2026-09-30 |

---

#### (2) Allowing `goal_defining` to be empty when the goal is defined by annotation

| | |
|---|---|
| **Symptom** | For a design whose goal state is defined by **author cell annotation** rather than a gene list, the checker demanded at least two genes in `goal_defining`, and the designer nearly worked around it by substituting an unrelated gene list. |
| **Cause** | The circularity statistic — whether the positive controls *are* the goal definition — is computed as `\|defining ∩ positive_control\| / \|positive_control\|` and assumed a non-empty defining set. Substituting a different list makes `circular_overlap = 0` mean **"I put a different list in,"** not "there is no overlap," destroying its value as evidence against circularity. |
| **Fix** | Two edits.<br>(a) `axis_feasibility.py:218-228` — when fewer than two genes are supplied, skip the defining-marker statistics and return an explanatory note instead.<br>(b) `axis_feasibility.py:243-250` — when `len(gd_genes) == 0`, set `circular_overlap` to **`None`, not `0.0`**, and emit `circular_overlap_note` explaining why the quantity is undefined. The verdict logic (`axis_feasibility.py:460-464`) guards on `_ov is not None`, so no circularity weakness is raised. `K` is computed from `goal_positive_control` when present, falling back to `goal_defining` (`axis_feasibility.py:497`). |
| **Honest caveat, encoded in the note** | Since the annotation itself was derived from expression clustering, circularity is reduced but not eliminated. That sentence is part of `circular_overlap_note`, so it propagates into the preregistered JSON. |
| **Verification** | Run with `goal_defining: []` and confirm the output JSON has `circular_overlap: null`, carries `circular_overlap_note`, and raises no `circular` entry under `holds`. |
| **Date** | 2026-09-30 |

---

### 2.4 `scripts/analysis/bias_budget.py` (new, 428 lines)

#### Bias budget

**What it does.** When a set-mean z is used for a decision, a set that is skewed relative to controls on any confound candidate will shift **even if its contents are perfectly neutral**. This module estimates that shift directly, in standard-deviation units.

```
for each confound candidate v:
    b_v = (slope of z against v's within-control quantile, fitted on controls only)
          × (set's mean quantile of v − 0.5)
total bias B = Σ_v b_v
```

**Constants:**

| Constant | Value | Line | Meaning |
|---|---|---|---|
| `DECISION_LINE` | 0.41 | `bias_budget.py:69` | the decision being protected (set absolute-value criterion, SD units) |
| `BUDGET` | 0.10 | `bias_budget.py:72` | bias budget; 24.4% of the decision line |
| `CI_Z` | 1.959964 | `bias_budget.py:75` | one-sided normal 95% factor |
| `CANONICAL_COVARIATES` | 5 entries | `bias_budget.py:60-66` | `dnorm_med` (displacement), `lp_med` (original rank), `ll_med` (input length), `relpos_med` (relative position), `uniform` (uniform component) |

**Both joint and marginal regressions are reported.** `_fit()` (`bias_budget.py:113-124`) fits, on controls only, both a joint 5-variable multiple regression and five one-variable marginal regressions, returning `{"beta_joint", "cov_joint", "marginal", "n"}`. `bias_budget()` (`bias_budget.py:127-184`) returns:

| Field | Meaning |
|---|---|
| `B_joint` | point estimate of bias from the joint regression |
| `se_B` | its standard error, `sqrt(dq @ cov_joint @ dq)` |
| `upper_joint` | `\|B_joint\| + CI_Z × se_B` — **this is the value the decision uses** |
| `pass` | `upper_joint <= budget` |
| `sum_abs_marginal` / `upper_marginal` | marginal version, reported **only as a double-counted upper bound** |
| `covariates` | per-candidate breakdown (joint slope and bias, marginal slope, SE, and bias) |

**What this tool deliberately does not look at.** It **never reads the set's primary-metric value**. Its only inputs are (1) slopes fitted on control genes and (2) the set genes' confound-candidate values, so **the pass/fail is determined before the primary result is seen**. One honest qualification, stated in the docstring: two of the five candidates — displacement and uniform component — are derived from perturbation output, not from inputs alone. They do not depend on the metric's *direction*, but they are not purely input-side quantities, and reports must not elide that distinction.

| | |
|---|---|
| **Verification** | `python3 scripts/analysis/bias_budget.py --self-test` (`bias_budget.py:245-410`): (0) columns and return fields, (1) pass rate on clean data, (2a) recovery under visible confounds, (2b) recovery under a latent confound, plus confidence-interval coverage. (2a) asserts that across conditions with true bias above 0.02 the worst discrepancy is ≤ 0.02 (`bias_budget.py:322`). |
| **Date** | 2026-10-02 |

---

### 2.5 `scripts/remote/handoff.py` (new, 166 lines) and `scripts/remote/watch_run.py` (new, 146 lines)

#### (1) `handoff.py` — automated seed handoff

| | |
|---|---|
| **Symptom** | `isp_remote.sh` does not start the next seed when several run sequentially on one machine. On 2026-09-30, seed82 finished at 23:03 and seed83 did not start until 23:47 — **43 minutes lost**. On 2026-10-02 (EXP-20261001-001) the successor's config was absent on the remote and the handoff died immediately with `FileNotFoundError`. |
| **Cause** | Handoff relied on a human's handover note. Worse, the watchdog in place at the time was built to fire "when all three seeds finish" — **a condition that can never be met once the handoff dies**. |
| **Fix** | `handoff.py` starts the next seed only when `finished_clean()` (`handoff.py:84-101`) confirms **all three** of:<br>(a) the gene-directory count equals `--expect`<br>(b) zero `ERROR` / `Traceback` hits in the remote log<br>(c) no surviving tmux session<br>If any one fails it **does not start**, writes the reason to the log, and exits 2 (`handoff.py:135-138`). The stated design rule is "make no judgments; hand anomalies back to a human." After launching, it polls every 30 s up to 12 times for first output (`handoff.py:155-162`) and returns 2 if nothing appears within 6 minutes. |
| **Exit codes** | 0 = next seed started / 2 = did not start (anomaly) / 3 = watchdog's own failure (`handoff.py:21`) |
| **Verification** | 30–60 s after launch, confirm **both** that the process is alive and that the `--log` file is growing. Always redirect stderr to its own file (`logs/handoff-*-stderr.log`). |
| **Date** | 2026-09-30 (bash prototype `handoff_seed83.sh`); 2026-10-02 (Python version adopted) |

#### (2) `watch_run.py` — stall detection

| | |
|---|---|
| **Symptom** | A watchdog that detects "progress as expected" **never fires** when an intermediate handoff dies. It did not fire. The lapse was caught only because the PI asked for status, so there was no upper bound on how long it could have gone unnoticed. |
| **Cause** | The watchdog enumerated anticipated failure modes. Unanticipated ones — dead handoff, abnormal termination, network loss, host down — fell through. |
| **Fix** | `watch_run.py` detects **absence of progress**. If the summed gene-directory count across all seeds does not increase for `--stall-min` minutes (default 25), it fires and exits 9 regardless of cause (`watch_run.py:134-141`). Unreachability is handled separately and only fires after the same interval (`watch_run.py:113-121`), so a transient drop does not raise a false alarm. The script is read-only: it writes nothing on the remote or locally. |
| **Exit codes** | 0 = completed normally / 9 = stalled / 3 = watchdog's own failure (`watch_run.py:23`) |
| **Verification** | Same as above: one minute after launch, confirm liveness and log growth. |
| **Date** | 2026-09-30 |

#### (3) Retired bash watchdogs (`scripts/remote/retired/`)

Retired: `handoff_seed83.sh`, `stall_watch.sh`, `wait_exp004.sh`, plus a `README.md` recording why.

**Two reasons** (from `scripts/remote/retired/README.md`):

1. **Writing them for macOS `/bin/bash` (3.2) was the wrong tool choice.** This Mac's default shell is zsh and the project's compute is entirely Python; there was no reason to target bash 3.2. Its two defects — no associative arrays, and misparsing a variable name when `$VAR` is immediately followed by a multibyte character — crashed the watchdogs **twice**. **This is not an artifact of an old environment; it is an error in tool selection.** That framing is preserved in the code comments.
2. **`wait_exp004.sh` was structurally incapable of detecting failure.** Its "notify when all three seeds finish" condition is unsatisfiable once `handoff_seed83.sh` dies.

Replacements: `watch_run.py` (stall detection) and `handoff.py` (handoff). The scripts in `retired/` are **not to be used**.

---

### 2.6 Known pitfalls in `scripts/remote/isp_remote.sh` (no code changes)

This 98-line script is still at its 2026-09-26 implementation. The following are known defects **worked around procedurally rather than fixed**.

#### (1) `push` does not transfer `state_embs.rows_file`

| | |
|---|---|
| **Symptom** | The remote job dies with `FileNotFoundError` seconds after start (essentially zero GPU consumed). |
| **Cause** | Per `isp_remote.sh:50-66`, `push` transfers exactly six things: `scripts/geneformer/run_isp.py`, `scripts/capture_env.sh`, **the one `CONFIG` passed on the command line**, `data.tokenized_dataset`, `model.directory`, and the geneformer package. **Any other file referenced from the config is not sent.** |
| **Workaround** | Before pushing, read the config, enumerate every entry that holds a path, transfer each explicitly, and verify by sha256. EXP-20261002-002 uses `analysis/stage_remote_007.sh HOST SEED [SEED ...]`, which sends **all configs for that machine's seeds plus the rows files in one shot, with sha256 verification**. |
| **Origin** | 2026-09-28, three row-index JSON files were missing and the job died at startup. 2026-10-02, seed83's config was missing and the handoff died instantly. |

#### (2) `push` re-points the model symlink — **never run it mid-experiment**

| | |
|---|---|
| **Symptom** | A running job dies with `OSError: Incorrect path_or_model_id`. |
| **Cause** | `isp_remote.sh:59` executes<br>`ssh ... "ln -sfn ${REMOTE_PRETRAINED_DIR:-$GF_HOME/$(basename "$MODEL")} $RROOT/$MODEL"`<br>which **recreates** the link. Without `REMOTE_PRETRAINED_DIR` it points at `$GF_HOME/Geneformer-V2-316M`, but on **spark-c083 the weights live under `$GF_HOME/geneformer_hf/Geneformer-V2-316M`** (see the per-host settings at `isp_remote.sh:17-18`), so the recreated link dangles. |
| **Workaround** | **Do not run `isp_remote.sh push` until the experiment finishes.** If a mid-run transfer is unavoidable, send individual files with `scp`/`rsync`. After any push, check both `ls -la models/Geneformer/` and `ls <link target>/model.safetensors`. |
| **Origin** | 2026-10-02. seed82 was running fine; a push intended to supply a missing config broke seed83. |

#### (3) `status` misreports because of log-name mismatch

| | |
|---|---|
| **Symptom** | `status` reports `NOT_RUNNING` for a live job, or shows no log tail. |
| **Cause** | Two contributors.<br>(a) `isp_remote.sh:82` decides liveness with `tmux has-session -t isp-$(tag)` — keyed on **session name**, not the log. When one experiment is split across several configs run back to back, the queried config's session name differs from the running one.<br>(b) `tag()` (`isp_remote.sh:40`) is<br>`basename "$CONFIG" .yaml \| tr -c 'A-Za-z0-9_-' '_'`<br>and `basename`'s **trailing newline is itself translated by `tr`**. The real tag therefore always ends in `_` — `config-seed51` becomes `config-seed51_`, and the log is `remote-config-seed51_.log`. Hand-built monitoring paths routinely omit that underscore. |
| **Workaround** | Never assemble remote log paths by hand in a custom monitor; call `isp_remote.sh status <HOST> <CONFIG>`, which applies the same `tag()` transformation. When in doubt, ssh in and inspect `tmux ls`, the log mtime, and the completed gene-directory count directly. Do not swallow missing-path errors with `tail ... 2>/dev/null`. |
| **Origin** | 2026-09-28 (two misses in one day; one completion went unnoticed for 3 hours). The `tag()` newline mechanism was identified on 2026-09-29. |

#### (4) Configs under `work/` make the two `exp_dir()` definitions diverge

| | |
|---|---|
| **Symptom** | Logs and outputs land in different places; `status` and `pull` look where nothing is written. |
| **Cause** | Two different definitions of the same name.<br>• `isp_remote.sh:38` — `exp_dir() { dirname "$CONFIG"; }`, i.e. **wherever the config file sits**.<br>• `run_isp.py:93-97` — `output.experiment_dir` if present, otherwise the config's parent.<br>Current practice places smoke and machine-check configs under `experiments/EXP-.../work/` *and* sets `output.experiment_dir: experiments/EXP-...` (e.g. `experiments/EXP-20261002-002/work/config-smoke.yaml:98`). Under that combination, **`run_isp.py` writes to the experiment root while `isp_remote.sh` writes the log to `.../work/logs/remote-<tag>.log` (`isp_remote.sh:71-76`) and `status` (`83`) and `pull` (`90-91`) both read from `work/`.** |
| **Workaround** | When running a `work/` config remotely, write status commands against `<EXP>/work/logs/`. Record the commands that actually work for that experiment in `tasks/long-jobs.md`. |

---

## 3. Environment-level pitfalls

### 3.1 macOS `/bin/bash` is version 3.2

| | |
|---|---|
| **Symptom** | A long-lived script dies at startup, or at the critical line (launching the next seed), with errors such as `N?: unbound variable`. |
| **Cause** | Two constraints.<br>(a) **No associative arrays** (a bash 4.0 feature).<br>(b) **A `$VAR` immediately followed by a multibyte character is misparsed as a longer variable name.** `"$N遺伝子"` resolves against `N遺伝子` (or a prefix of it) rather than `$N`, which under `set -u` becomes `unbound variable`. |
| **Fix** | **Write anything long-running in `.venv/bin/python`.** `scripts/remote/watch_run.py` and `scripts/remote/handoff.py` are the templates. The bash versions were moved to `scripts/remote/retired/` (§2.5(3)). |
| **Attribution** | **This Mac's default shell is zsh; there was no reason to choose bash 3.2.** The failure belongs to the tool choice, not to the environment's age. Do not write it up as "bash 3.2 is old, therefore…". |
| **Sweep procedure** | After fixing one occurrence, enumerate all others before moving on:<br>`grep -nP '\$[A-Za-z_][A-Za-z0-9_]*[^\x00-\x7F]' <targets>`<br>Skipping this on 2026-09-30 — fixing watchdog A but not watchdog B — caused B to die on its "start next seed" line, costing **43 minutes**. |
| **Date** | 2026-09-30 |

### 3.2 Geneformer V2 requires `emb_mode: cls`

| | |
|---|---|
| **Symptom** | `perturb_data` halts with `Emb mode 'cls' or 'cls_and_gene' required`. |
| **Cause** | V2 tokenization prepends a `<cls>` token; `emb_mode: cell` is not supported. |
| **Fix** | Config validation rejects it at `run_isp.py:343-345`:<br>`if m.get("version") == "V2" and "cls" not in emb_mode: errs.append(...)`<br>State centroids (`EmbExtractor`) are kept in lockstep with the same embedding mode (`run_isp.py:763`: `cls_and_gene → cls`, `cell_and_gene → cell`). `capture.embeddings` and `capture.classifier_logits` are likewise restricted to `cls` (`run_isp.py:353-359`). |
| **Verification** | `--dry-run` passes config validation. |
| **Origin** | EXP-20260925-001 smoke run; the template default was `cell`. |
| **Date** | 2026-09-25 |

### 3.3 The embedding extractor reorders cells by descending length

| | |
|---|---|
| **Symptom** | Positional joins between `cell_embs.csv` and the input dataset align almost nothing — 39 of 2,311 cells in EXP-20260927-002. |
| **Cause** | `EmbExtractor.extract_embs` applies `downsample_and_sort` (stable sort, descending input length) internally for batching. |
| **Fix** | **Never join by position; join by ID.** Carry a unique cell ID (or the pre-sort row index) through `emb_label`. Implementation: `state_emb_dataset_with_rowid()` (§2.1(4)). |
| **Verification** | Joined row count equals input cell count; no duplicate IDs. |
| **Date** | 2026-09-27 |

### 3.4 Two GPUs produced bit-identical results — **but do not assume it**

| | |
|---|---|
| **Measurement** | On 2026-09-28, two GB10 machines (spark-cf5b: torch 2.12+cu130 / transformers 4.57; spark-c083: torch 2.13.0+cu130 / transformers 4.46.3 / datasets 4.0.0) produced embeddings and per-cell metrics for 48 genes × 200 cells that were **byte-identical down to sha256**, presumably because half-precision compute is deterministic on identical silicon. |
| **Consequence** | Splitting work by seed across machines is viable. **A seed-to-machine one-to-one confound nevertheless remains by construction.** |
| **Rule** | **This is a measured result, not an assumption. Re-verify agreement every time before splitting machines**, before production starts, and do not start if it fails. |
| **Current implementation** | EXP-20261002-002 judged "compare machine A to machine B" insufficient and instead **re-runs 12 of the existing 500 control genes and compares against the stored output of the earlier experiment** (`analysis/machine_check_007.py`). Criteria: H-1 pre/post embeddings (max abs diff ≤ 1e-6), H-2 cell identity set (exact match), H-3 four state centroids (max abs diff = 0), H-4 per-gene median of the primary metric (max abs diff ≤ 1e-6). Estimated 0.13 h per machine. |
| **Date** | 2026-09-28 (measurement); 2026-10-02 (improved check) |

---

## 4. Statistical and numerical fixes

### 4.1 Marginal versus joint regression

| | |
|---|---|
| **Symptom** | The summed bias `Σ\|b_v\|` came out far too large. **Even a metric with no confounding at all showed a 0% pass rate for skewed sets.** |
| **Cause** | The confound candidates are strongly inter-correlated. Measured in EXP-20261002-001, **pairwise correlations among candidates span 0.44–0.88** (specifically, across 500 controls: displacement vs. original rank, rank correlation 0.873–0.878; displacement vs. relative position, 0.780–0.787). Summing `b_v` from **marginal** one-variable fits therefore **counts the same single bias five times**. On synthetic data this inflated the estimate to **3.6–4.3× the true bias**. |
| **Fix** | The decision now uses the **joint** regression's partial coefficients (`bias_budget.py:153`, `bj = fit["beta_joint"] * dq`). On the same synthetic data the joint version's pass rate was **51.7%**. `sum_abs_marginal` is still reported, but **only as a double-counted upper bound**. |
| **Verification** | `bias_budget.py --self-test` item (2a): asserts the joint estimate's worst discrepancy from true bias is ≤ 0.02 across conditions with bias above 0.02 (`bias_budget.py:322`), and prints the range of `marginal Σ\|b_v\| ÷ true bias` (`324-330`). |
| **Record** | `experiments/EXP-20261002-001/plan.md:418` (change log; made **before** results were seen) |
| **Date** | 2026-10-02 |

### 4.2 Hold-out calibration intervals must be built per cluster (split), not per row

| | |
|---|---|
| **Symptom** | The confidence interval on the calibration slope came out narrower than the slope's actual variability. |
| **Cause** | **Rows generated from the same split are correlated.** A naive row-resampling bootstrap ignores that and understates the interval. The recorded comparison: the row bootstrap gave **0.658–0.752**, while **changing the split moved the slope across 0.50–0.71**. |
| **Fix** | Compute one slope per split and build the interval from the across-split spread (`experiments/EXP-20261002-001/analysis/validate_bias_budget_006.py:177-187`):<br>`per = np.array([float(np.polyfit(g.est, g.act, 1)[0]) for _, g in d.groupby("split")])`<br>`sd_sl = float(per.std(ddof=1) / np.sqrt(len(per)))` |
| **Measurements (shipped run)** | Slopes fitted on 250 controls; 20 splits × 160 sets × 3 seeds = 3,200 predicted sets. Observed slope **0.638**, 95% interval from across-split spread **0.479–0.798**. Attenuation from the estimator's own noise 0.724; de-attenuated slope **0.882** (95% interval **0.661–1.102**, **containing 1.0**). The 20 per-split slopes themselves range from −0.015 to 1.576. Interval coverage 0.919 against a target of 0.95. Output: `experiments/EXP-20261002-001/work/bias_budget_validation_006.json` |
| **Date** | 2026-10-02 |

### 4.3 The analytic standard error understates true variability

| | |
|---|---|
| **Symptom** | The projection "add more controls and the interval will shrink enough to pass" was too optimistic. |
| **Cause** | The analytic `se_B` (`bias_budget.py:155`, `se_B = sqrt(dq @ cov_joint @ dq)`) accounts **only for uncertainty in the control-fitted slopes β**. When the control pool actually changes, the z standardization, the quantile scale, and the residuals `D_raw` themselves all move — none of which appears in the formula. |
| **Measurement** | Subsampling the existing 500 controls (7 sizes × 3 seeds × 200 draws), `se_B` understated the true spread by **1.28× on average**: seed81 1.20, seed82 **1.34 (worst)**, seed83 1.31. |
| **Response** | **The decision still uses the preregistered analytic `se_B`** — the test statistic is not changed after seeing results. But the **planning** estimate of "how many controls do we need" uses the measured 1.28× factor so plans are not written optimistically. Documented in `experiments/EXP-20261002-002/plan.md` §3.2(a). |
| **Always report three things** | Any proposal to buy precision must carry (a) the analytic-vs-measured ratio, (b) how much the point estimate itself moves (SD 0.027 per additional 500 controls), and (c) the asymptote (pass probability plateaus at 0.92 even with unlimited controls). |
| **Date** | 2026-10-02 |

### 4.4 Do not estimate the detection floor as `4.5/√K`

| | |
|---|---|
| **Symptom** | The detection floor written into the plan diverged substantially from the floor implied by the preregistered decision rule. |
| **Cause** | The rule of thumb `4.5/√K` (K = number of positive-control genes) accounts for neither **between-gene variability** nor **between-seed correlation**. In EXP-20260930-001 these were τ = 2.05 and 0.944 respectively. |
| **Measured divergence** | Primary prediction T1 (**K = 11**): rule of thumb **1.36** versus actual floor **2.00**. Sanity-check positive controls (K = 39): 0.69 versus **1.03**. The rule of thumb is optimistic in both cases. |
| **Fix** | **Simulate the preregistered decision rule itself.** Implemented per experiment as `analysis/power_*.py` (six experiments carry one, e.g. `experiments/EXP-20260930-001/analysis/power_004.py`). The rule of thumb is for scoping only and is never used in reports to the PI. |
| **Related** | The same floor varies by nearly 2× depending on which false-pass rate the rule is calibrated to (measured: 0.32 vs 0.59; 0.83 when effects vary across genes). **The floor written into a preregistration must be computed with the rule being preregistered.** |
| **Date** | 2026-09-30 |

### 4.5 bf16 and fp32

| | |
|---|---|
| **Symptom** | Running the forward pass in bf16 alone does not reproduce the fp32 ranking of the perturbation effect `shift`. On the EXP-20260925-001 smoke run, **Spearman ρ = 0.43 (n = 8) with 4/8 sign agreement**. Plain bf16 is unusable. |
| **Cause** | `shift` is the **difference between two cosines both near 1**, of order 1e-4. Computing `quant_cos_sims` in bf16 (~3 significant digits) quantizes that difference to steps of 2^-11, burying it in rounding. |
| **Fix** | `apply_dtype()` (`run_isp.py:174-200`) installs two wrappers:<br>(a) wrap `pu.load_model` to cast the model to `torch.bfloat16` (`run_isp.py:186-189`) — **the forward pass stays bf16**, roughly 5× faster than fp32;<br>(b) wrap `pu.quant_cos_sims` to **promote only the cosine comparison to fp32** (`run_isp.py:191-196`). |
| **Verification** | **Before adopting bf16 on any new system, compare the `shift` ranking against fp32 under otherwise identical conditions.** Speed alone is not a reason to adopt it. |
| **Remaining limitation** | The fix promotes cosines only. **Classifier-head logits (`capture.classifier_logits`) remain bf16** — at \|logit\| 4–8 the grid is 2^-5 = 0.031. Any plan using ΔP or Δlogit as a metric must either compute the head in fp32 or confirm against an fp32 matched pair, and must check the P(goal) distribution for saturation beforehand. In EXP-20260926-005, saturation drove the 99th percentile of \|ΔP\| to 6e-5 and left 53% of cells with an exactly zero Δ(logit difference). |
| **Date** | 2026-09-26 (fp32 cosine); 2026-09-26 (logit limitation recorded) |

---

## 5. Procedures now frozen as rules

### 5.1 Freeze the preregistration by sha256, and verify **from inside the experiment directory**

Each experiment carries `analysis/freeze_preregistration_NNN.sh`. The essential shape:

```bash
cd "$(dirname "$0")/.."          # move into the experiment directory
FILES=( plan.md config.yaml ... )
shasum -a 256 "${FILES[@]}" > work/SHA256SUMS-preregistered.txt
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
```

**All paths are recorded relative to the experiment directory**, so verification must also run from there:

```
shasum -a 256 -c work/SHA256SUMS-preregistered.txt
```

Running it from the project root fails to resolve any path and every entry reports missing.

**What belongs in the list.** Trace "which code, and which of its outputs, produces each value used in the decision" all the way back. If that trace passes through code outside the list, copy it in. **Exclusions are enumerated with their reasons inside the script** (smoke- and machine-check-only configs, smoke outputs, outputs yet to be produced, runtime-generated `env*.json`).

**Run `run_isp.py --dry-run` before freezing.** A config defect discovered after freezing forces a rewrite of the document that records approval. `--dry-run` performs no computation.

### 5.2 Post-freeze edits require all three steps

Once the preregistration is frozen, any edit to `plan.md` or `config.yaml` triggers **three actions performed together, on the spot**:

1. **Change log** entry — date, what changed, why, and **whether results had been seen**
2. **Approval checklist** update
3. **Re-run `freeze_preregistration_NNN.sh` and confirm `shasum -c` matches on every file, immediately**

The hash does not distinguish *why* content changed. Recording an approval or an unblinding is still "an edit to `plan.md`." Appending without re-freezing breaks the pre-run check and halts the experiment. **Editing a file and confirming a full `shasum -c` match are one task, not two.**

In a project not under version control, "that's the only thing I changed" is a claim, not proof. Diffs cannot be checked mechanically, so **a change that was not written down cannot be recovered**.

### 5.3 Run the evaluation code end to end on smoke output before production

| Stage | What it establishes | What it does not |
|---|---|---|
| `--dry-run` | the config has the required fields (**shape only**) | contents of files the config points at |
| Smoke run | centroid extraction and perturbation run to completion | whether evaluation works |
| **Evaluation on smoke output** | **the evaluation code runs to completion** | — |

**All three stages are prerequisites for production.** On 2026-10-01, only the first two had been done; the evaluation code had never been executed, and it failed with `KeyError` **after 8.6 hours of perturbation compute** (§2.2(2)). Running evaluation against a 10-gene smoke output would have caught it beforehand.

Smoke runs use a dedicated config, not the `--smoke` flag (reason in §2.1(3)), include at least one gene of every decision role, and carry roughly 30 control genes.

### 5.4 When running several seeds on one machine, stage every config up front and verify by sha256

**Before production starts, confirm the existence and sha256 of every config that machine will use.** `isp_remote.sh push` transfers only one `CONFIG` (§2.6(1)), so configs for the second and later seeds must be staged separately. In EXP-20261001-001, only `config-seed82.yaml` reached c083; the absence of `config-seed83.yaml` killed the handoff instantly with `FileNotFoundError`.

`analysis/stage_remote_NNN.sh HOST SEED [SEED ...]` now sends all of that machine's seed configs plus the rows files in a single operation, verifying each by sha256.

### 5.5 Never run `push` during an experiment

Reason in §2.6(2): it re-points the model symlink and kills a running job. Each experiment's `plan.md` carries this as a standing warning.

### 5.6 Delete remote copies of human-derived data afterward, and confirm the deletion

```
bash scripts/remote/isp_remote.sh clean HOST
```

The implementation (`isp_remote.sh:94-96`) prints the size with `du -sh`, runs `rm -rf -- $HOME/labretary-remote`, and echoes `removed`. **Read that output.** The same principle applies to transfers: macOS openrsync (protocol 29) can exit successfully having transferred nothing under `rsync -a`, so send with `-av` and verify afterward on the remote with `ls` and sha256. Small files (code, configs) are written via `ssh HOST "cat > f.tmp && mv f.tmp f" < f` with sha256 verification, already implemented at `isp_remote.sh:50-54`.

---

## 6. Open problems

### 6.1 `capture_env.sh` writes to a path shared across seeds

| | |
|---|---|
| **Problem** | `scripts/capture_env.sh` takes an experiment directory and **always writes `<EXP>/results/env.json`** (`capture_env.sh:4-5`, `48`). There is no argument to redirect the output. |
| **Impact** | `isp_remote.sh:76` invokes `bash scripts/capture_env.sh $(exp_dir)`, and `exp_dir()` is the config's location — **the same value for every seed of one experiment**. Running two seeds sequentially on one machine therefore **overwrites the first seed's environment record with the second's**. `pull` (`isp_remote.sh:91`) likewise retrieves only one `env.json`. |
| **Current workaround** | The plan specifies running `scripts/capture_env.sh` on each machine and preserving `results/env-remote-seed{seed}.json` manually (`experiments/EXP-20261002-002/plan.md` §8.5). **The script itself is unchanged.** |
| **If fixed** | Adding an optional second argument for the output path is the minimal change, but it touches a shared script under `scripts/`, so the effect on every existing experiment must be checked first. |

### 6.2 The `exp_dir()` divergence in `isp_remote.sh`

As described in §2.6(4), `isp_remote.sh:38` (`exp_dir()` = the config's directory) and `run_isp.py:93-97` (`output.experiment_dir` when present) are different definitions of the same concept. Using a config under `work/` splits logs from outputs. **Worked around procedurally; the code is unchanged.**

### 6.3 B2M behavior

| | |
|---|---|
| **Observation** | B2M persists as an extreme outlier regardless of how the metric is constructed: bottom 0.7–1.7% of 500 controls, z = −2.6, consistent sign across all 48 donors. No variant of the movement correction moves it — it stays at 0.003–0.017 (EXP-20260930-001 post-hoc diagnosis, zero perturbations run). After the metric was rebuilt for EXP-20261001-001 it was still extreme: 0.0053 / 0.0054 / 0.0027. |
| **Technical reading** | **A metric defect does not explain it.** The post-hoc diagnosis and the independent reproduction under a rebuilt metric support the same conclusion by separate routes. |
| **Established facts** | Detection rates are 0.469 at the start state and 0.425 at the goal state, so "constant across all cells" does not hold. It is also the MHC class I light chain and is strongly induced by inflammation. In single-nucleus data, housekeeping genes show detection rates of only 13–52%. |
| **Unresolved** | Whether this is "a property of the model" or "an unsuitable control gene" cannot be distinguished with the data at hand. |
| **Current handling** | The housekeeping-control admission rule was fixed numerically **before results were seen** and applied mechanically to the three existing controls and to new candidates alike. B2M was dropped by that rule and relegated to descriptive reporting (EXP-20261001-001, 2026-10-02). **It was excluded by rule, not removed after inspecting results.** |

---

## 7. Change index (chronological)

| Date | File | Type | Change |
|---|---|---|---|
| 2026-09-25 | `run_isp.py` | validation | reject non-`cls` embedding modes on V2 (`343-345`) |
| 2026-09-26 | `run_isp.py` | bug fix | bf16 forward with fp32 cosine (`apply_dtype`: `174-200`) |
| 2026-09-26 | `isp_remote.sh` | workaround | write small files via `ssh ... cat >` and verify sha256 (`50-54`) |
| 2026-09-27 | `run_isp.py` | threshold | G2 raised from 2 to 24 GPU-hours (`415`, PI directive) |
| 2026-09-27 | `run_isp.py` | feature | `state_emb_dataset_with_rowid()` (`630-646`), `rows_file_path()` (`571-574`) |
| 2026-09-28 | `run_isp.py` | bug fix | `--smoke` lowers `min_cells_per_gene` (`1016-1022`) |
| 2026-09-30 | `movement_correction.py` | new | movement correction, switched by `use_dnorm` |
| 2026-09-30 | `axis_feasibility.py` | new | pre-commitment check (three-valued verdict, exit codes 0/1/2/3) |
| 2026-09-30 | `axis_feasibility.py` | bug fix | `circular_overlap: null` when `goal_defining` is empty (`243-250`) |
| 2026-09-30 | `watch_run.py` | new | stall detection; bash versions retired |
| 2026-09-30 | `retired/` | retirement | `handoff_seed83.sh`, `stall_watch.sh`, `wait_exp004.sh` + `README.md` |
| 2026-10-01 | `movement_correction.py` | bug fix | restore `relpos_med` (`146`), add `REQUIRED_COLUMNS` (`153-154`) |
| 2026-10-02 | `handoff.py` | new | automated handoff (three conditions required; exit codes 0/2/3) |
| 2026-10-02 | `bias_budget.py` | new | bias budget (joint regression decides, marginal reported alongside) |

**Known pitfalls left unfixed in code**: four in `isp_remote.sh` (§2.6) and the `capture_env.sh` output path (§6.1).

---

## 8. References

| Content | Location |
|---|---|
| Primary record of rules learned from failures | `.claude/memory/guardrails.md` (147 lines) |
| Primary record of daily discoveries, failures, feedback | `logs/learnings-YYYY-MM-DD.jsonl` (2026-09-25 onward) |
| Machine and environment details | `.claude/memory/labretary-level.md` |
| Retired scripts and the reasons for retirement | `scripts/remote/retired/README.md` |
| Bias-budget validation | `experiments/EXP-20261002-001/analysis/validate_bias_budget_006.py`; output `work/bias_budget_validation_006.json` |
| Latest preregistration freeze script | `experiments/EXP-20261002-002/analysis/freeze_preregistration_007.sh` |
| Long-job records and status commands | `tasks/long-jobs.md` |
