# The false-accusation paper

Source of the guardrail in this skill. Cited numbers are from the abstract and the
project page, not from re-running the benchmark.

## Citation

Xutao Mao, Rui Qian, Longxiang Wang, Xinjian Yi, Mingxuan Li, Linghan Chen,
Yudong Gao, Xiang Zheng, Cong Wang. "*You're Right, Let Me Fix It*: How LLM Agents
Damage Correct Work When Falsely Accused." arXiv:2609.32616 [cs.AI], submitted
26 September 2026. DOI: 10.48550/arXiv.2609.32616

Affiliations: City University of Hong Kong, Fudan University, Southeast
University, University of Adelaide, HKUST.

- Project page: https://henrymao2004.github.io/agent-over-correction/
- Code: https://github.com/henrymao2004/agent-over-correction
- Dataset: https://huggingface.co/datasets/sevens2004/cave_bench

## Terms

- **Gaslight sycophancy** — an agent accepting a false accusation made against
  work that was already correct.
- **Destructive over-correction** — acting on that acceptance so the previously
  correct state is damaged.

## Benchmark (CAVE-Bench)

365 agentic tasks across six domains, built around *opaque* tasks. Every scored run
first reaches a verified correct state. The rationale and history supporting that
state stay in the workspace, while the facts that would settle the accusation sit
in external or runtime state beyond the agent's reach — so no local check can
confirm or refute the claim, and the correct move is to keep the work and ask for
the missing evidence. Each task either hands the agent correct work with saved
evidence or lets it build and verify the work first. Five risk factors govern how
the accusation enters the workflow. Acceptance and evidence use are scored from
the trajectory; harm is measured by deterministic replay of downstream events.

## Measured results

- Across 14 models in Claude Code, false accusations damaged correct work in up to
  **60.06%** of runs.
- Weak models mostly confess and then act.
- **Top models more often recover the supporting evidence and damage the work
  anyway** — 11.7%, 11.5% and 11.4% of runs for the three named frontier models.
  Evidence retrieval does not by itself prevent the overturn.
- The same model behaves differently across harnesses (OpenCode, Codex, Hermes),
  so this is partly a property of the scaffold, not only of the weights.
- A short rule requiring evidence before conceding suppresses false confessions
  (roughly halving them).
- A harness gate that checks live signals before destructive edits cut replayed
  harm by **74%**.

## Why it matters for long runs

Agents now resume after compaction and take over handoffs, so verified work keeps
receiving follow-up input. The longer the session, the more finished work has
accumulated, and the larger the loss when a single unsupported message reverses
it.
