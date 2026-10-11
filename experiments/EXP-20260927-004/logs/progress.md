# EXP-20260927-004 進行ログ

## 2026-09-27 (experiment-planner)
- 20:2x 作業開始。EXP-20260927-003 の validation.md / plan.md / config.yaml / evaluate_fib.py、guardrails、labretary-level、tracker を読了。
- テンプレートから experiments/EXP-20260927-004/ を作成（plan.md, config.yaml, analysis/, work/, logs/, results/）。
- literature-scout に「目標の軸での陽性対照」の文献調査を委託（Q1 線維芽細胞→三次リンパ組織オーガナイザー様への転換を駆動する遺伝子、Q2 代替の確立した単一遺伝子転換、Q3 オーガナイザー様のマーカー）。
- 既存データ（EXP-20260927-003 results, 84 遺伝子 × 3 seed）で計算できる指標候補の事後評価を準備。**主対象 LTBR/RELB/NFKB2 と CDKN1A は読み込み直後に除外**し、対照遺伝子だけで評価する。
- 利用可能なセントロイド: FIB_org（検出された遺伝子数の分布をそろえたもの, goal）, FIB_neg_TLShigh（start）, FIB_neg_TLSlow（alt）, MYOF_neg（状態の対照）, FIB_org_unmatched（感度分析）。細胞ごとの元・摂動後の CLS 埋め込み（orig.npy / pert.npy）も保存済み → 方向の射影指標も計算できる。
- 指標候補の事後評価を完了（`analysis/metric_candidates.py` → `work/metric_candidates.{json,csv}`）。12 候補。主対象 4 遺伝子を除いた 80 遺伝子 × 3 seed（22,614 行）で計算。
  - 現行 A0 の再現を確認（一様成分との ρ 0.974〜0.985、置き換え 0.89〜0.92、交絡 ll −0.489）＝ validation.md と一致。
  - 最良は **A11 = 目標 − (出発点・対照の状態・筋線維芽細胞の平均)** と **A12 = A11 / 移動量**。一様成分との ρ 0.106、置き換え 0.25〜0.26、交絡 ≤0.244/0.220、遺伝子固有の分散割合 0.83/0.82（A0 は 0.32）。
  - 筋線維芽細胞の軸の陽性対照は A12 で ACTA2 7.8〜9.1 標準偏差、TAGLN 3.4〜4.2（A0 は 4.0〜4.8 / 1.2〜1.4）。
- literature-scout 完了（`outputs/literature/TLS-organizer-drivers-2026-09-27.md`）。**Q1 は該当なし**: 線維芽細胞で単一遺伝子を過剰発現して三次リンパ組織オーガナイザー様へ転換させた一次文献は見つからなかった。LTβR はアゴニスト刺激（受容体シグナル）での十分性（PMID 22940098）とノックアウトでの必要性（PMID 23623380）のみ。しかも LTβR 経路は本仮説の主対象そのもので、陽性対照に使うと循環する。
  → 代替（マーカー操作型の陽性対照）を設計する。

## 2026-09-27 (experiment-planner / 再開 2 回目)
- 前担当は API 上限で停止。再開時点の実地確認（`find` と `ls -lT`）:
  - 完了: `analysis/metric_candidates.py` + `work/metric_candidates.{json,csv}`（21:51）、`analysis/select_genes_004.py`、
    `work/{positive_controls_004.csv, null_pool_004.json, state_rows-seed3{1,2,3}.json, unmatched_goal_rows-seed3{1,2,3}.json,
    depth_match_summary_004.json, precheck_coverage_004.json, tokenized.dataset}`（21:56）。
  - 未完了: `work/genes_004.json` が無い。`select_genes_004.py` は最後の
    `assert nc.min() >= MIN_CELLS_PER_GENE` で落ちている（**CA8 は 3 seed とも 111 細胞 < 下限 120**）。
  - `plan.md` と `config.yaml` はテンプレートのまま（未着手）。`analysis/evaluate_004.py` も無い。
- `precheck_coverage_004.json` の追加の不適合（主検定の下限 20 を下回る「使える対照遺伝子」）:
  IL7 17/19/17、CA8 9/15/15、HPSE2 14/10/14、FNTB 19（seed32）、ZNF587B 19（seed33）。
  いずれも元の位置（高発現側）が対照遺伝子の分布から外れていることによる。
- 再開後にやること: (1) CA8 の差し替えと対照遺伝子の補充 → `genes_004.json`、(2) 実データからの検出力見積もり、
  (3) `analysis/evaluate_004.py`（EXP-003 の副検定 (d) の不具合を引き継がない）、(4) plan.md / config.yaml / sha256。

### 完了（再開 2 回目・同日）
- `select_genes_004.py` を修正して完走。変更点:
  - data 層の陽性対照を固定リストから**機械的な規則**（`DATA_TIER_RULE`）に変更。実行可能性（抽出で 120 細胞以上、
    補充後に使える対照遺伝子 20 以上）を選定条件に組み込んだ。→ **CA8 は 111 細胞で不採用**（`positive_controls_rejected_004.json`）、
    HPSE2 は順位 10 位で上位 5 に入らず。data 層は FKBP5 / SLC2A3 / CCN1 / SYN3 / STEAP4。
  - 使える対照遺伝子が 20 未満の**判定対象の遺伝子**に対照遺伝子を自動補充する処理を追加（IL7 に 9 個）。
    対照遺伝子自身が 20 未満の場合は較正（calibration。計測分野の「校正」(JIS Z 8103) とは別語）の欠測として扱う（guardrails: 欠測を不合格と数えない）。
  - 結果: 摂動する遺伝子 **106**（陽性対照 10 / 筋線維芽細胞の陽性対照 2 / 恒常発現遺伝子 3 / 中立対照パネル 7 / 対照 84）。
    判定対象 22 遺伝子はすべて使える対照遺伝子 20 以上。較正に使える対照×seed は 248/252。
    延べ 63,570 細胞×遺伝子、総トークン 1.248e8。
- `analysis/power_004.py` → `work/power_004.json`。**実データのばらつきのみ**（合成データ不使用）。
  - **重要**: A12 では対照のばらつきの 81.5% が遺伝子固有の再現成分（A0 は 32.1%）。順位だけの規則では
    偶然合格率が下がらない（1 遺伝子 0.104、候補 10 個で 0.665）。→ **2.5 標準偏差の条件 ＋ 3 seed すべて**を採用（0/48）。
  - 最小検出可能効果量（200 細胞・規則込み・検出力 80%）= **+9.01e-5 = 2.97 標準偏差**。参考: ACTA2 の筋線維芽細胞の軸での効果は 7.76/7.82/9.09 標準偏差。
  - 細胞数 100→200 で対照の標準偏差は 5% しか縮まない（3.19e-5 → 3.04e-5）。検出力は細胞数では上がらない。
- `analysis/evaluate_004.py` を新規作成。EXP-003 の副検定 (d) の不具合は引き継がず、`compare_same_scale()` で
  尺度違いの比較を停止させる（自己検査済み: 同尺度は通過、異尺度は停止）。
- `plan.md` と `config.yaml` を作成。Ensembl ID 22 件を `gene_stats_fib.csv` と `input.h5ad` の `gene_name` で照合（全件一致）。
- `work/SHA256SUMS-preregistered.txt` に **17 ファイルを固定し、`shasum -c` で照合（全件 OK）**。
- learnings 4 件を `logs/learnings-2026-09-27.jsonl` へ（json.dumps で追記し、全 23 行を読み直して検証。壊れた行 0）。
- **摂動は実行していない。** G2 は不要（GB10 で約 4.0 時間 < 上限 24 時間）。

## 2026-09-27 22:41〜22:50 (simulation-engineer)
- 前提確認: `shasum -a 256 -c work/SHA256SUMS-preregistered.txt`（実験ディレクトリ内で実行）→ 17 ファイル全 OK。
- `scripts/capture_env.sh` でローカルの `results/env.json` を保存（Mac, MPS, torch 2.14.0, geneformer 0.1.0）。
- リモート計算機 2 号（coi-next1@100.94.57.17, spark-c083, GB10）へ push（サンドボックス外で実行）:
  - `config.yaml` / `run_isp.py` / `capture_env.sh` / `genes_004.json` / `null_pool_004.json` /
    `state_rows-seed{31,32,33}.json` を sha256 照合つきで送付。`tokenized.dataset` は rsync 後にリモートで
    `ls -la` して local と size/mtime が一致することを確認。事前学習モデルの重みも sha256 一致を確認。
- **`--dry-run` で停止を検出。config.yaml に必須の `stats:` セクションが無い**（`run_isp.py` の `validate()` が
  `セクション 'stats' がありません` で `sys.exit(1)`。実行前のため計算・回収は一切していない）。
  - `EXP-20260927-003/config.yaml` には `stats: {mode: goal_state_shift, fdr_threshold: 0.05}` があるが、
    本実験の `config.yaml`（sha256 固定済み）には無い。`run_stats()`（各遺伝子の ISP 後に必須で呼ばれる）が
    `cfg["stats"]["mode"]` を参照するため、`stats` が無いと validate() の時点で必ず止まる。
  - `evaluate_004.py` はこの `stats` サブディレクトリの出力を一切読んでいない（`cells.csv` と
    `genes/<ENSG>/emb/{orig,pert}.npy`、`state_emb_vectors.npz`、`extra_centroids-seed*.npz` のみ使用）。
    したがって値そのもの（`mode: goal_state_shift, fdr_threshold: 0.05`＝EXP-003 と同一）は M1〜M5 の
    判定基準・遺伝子リスト・細胞選抜のいずれにも影響しない、パイプラインの技術的な必須項目に見える。
  - ただし config.yaml は sha256 固定済みの事前登録ファイルであり、変更するには plan.md の変更ログ・
    承認チェックリスト・ハッシュの再固定が要る（guardrails 2026-09-27）。**判定基準や設計を変える権限は
    simulation-engineer にないため、ここで実行を止めて報告する。** 自己判断で追記・再ハッシュ・再実行はしていない。
  - リモートの複製（tokenized dataset 等、ヒト由来データ由来）は `isp_remote.sh clean` で削除し、
    `ls` 相当（ssh 経由の存在確認）で消えたことを確認した。
- **状態: 停止（config.yaml の不備待ち）。摂動は 0 遺伝子×0 seed。GPU 計算は未着手。**

## 2026-09-27 22:58〜22:59 (simulation-engineer / 再開)
- Labretary から config.yaml 修正済みの連絡。**自分でも独立に再確認**（信用せず検証）:
  - 実験ディレクトリ内で `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` → 17/17 OK。
  - `config.yaml` を読み、`stats: {mode: goal_state_shift, fdr_threshold: 0.05}` の追記を確認。
  - python で遺伝子集合を機械的に検算: `genes_to_perturb`(10) + `positive`(2) + `negative_fixed`(3) +
    `extra_fixed` パネル(7) + `negative_sampler.genes`(84) = **106、重複 0、`blocked_genes` との重なり 0**。
  - ローカルで `run_isp.py --config ... --dry-run` → `forward_passes: 63600`（= 106×200×3, 計画と一致）、
    `dry-run OK`。内蔵見積もり `gpu_hours_estimate: 0.12` / `measured_hours_total: 2.33` は
    guardrails 訂正に従い**採らない**（実測モデルの 4.0 時間を上限判断に使う）。
- リモート計算機 2 号（coi-next1@100.94.57.17, spark-c083）へ再 push（config.yaml, run_isp.py,
  capture_env.sh, state_rows-seed{31,32,33}.json を sha256 照合つきで送付。tokenized.dataset とモデル重みは
  クリーン後の再送要否を ls で確認し、既存ファイルのサイズ・mtime が一致していたため再送不要と判断）。
  - genes_004.json / null_pool_004.json はもう config.yaml から参照されない（遺伝子リストを直接展開したため）ので送っていない。
- リモートでも `--dry-run` を実行 → `forward_passes: 63600`, `dry-run OK`（`WARN compute.device=mps` は
  想定どおり。`isp_remote.sh run` は `LABRETARY_DEVICE=cuda` で上書きする）。
- **`isp_remote.sh run` で tmux 起動。実開始 2026-09-27 22:59:30 JST。**
  中断規則: seed 31 完了時点で 2.0 h 超（実開始から 01:00 頃）なら停止・報告。全体 6.0 h 超見込みでも同様。

## 2026-09-27 22:59:35〜22:59:36 (simulation-engineer) — 3 回目の停止（実行中に自動クラッシュ）
- `isp_remote.sh status` を実行したところ `NOT_RUNNING`。ログを見ると、起動 6 秒後（22:59:35〜36）に
  **`run-31-L-1` の状態のセントロイド抽出（`state_emb_dataset()`）で `KeyError: 'rows'` が発生し、tmux セッションごと終了していた**。
  - 原因: `config.yaml` の `state_embs.rows_file` は `work/state_rows-seed{31,32,33}.json` を指すが、
    この 3 ファイル（sha256 固定済み・事前登録）の中身は `{"FIB_org": [...], "FIB_neg_TLSlow": [...], "FIB_neg_TLShigh": [...]}`
    という**フラットな辞書**である。一方 `run_isp.py` の `state_emb_dataset()` は
    `json.loads(rf.read_text())["rows"]` として読むので、**トップレベルに `rows` キーで包んだ形**
    （EXP-20260927-003 の `work/state_rows-seed21.json` は `{"seed":21,"depth_match":...,"rows":{...},...}` の形）
    を要求している。**3 seed とも同じ形式不一致で、1 つも state セントロイドを作れない。**
  - 実害の確認: GPU 上の実計算（forward pass）はゼロ。`results/env.json` と `results/device.json` は
    書かれたが、`results/genes/` 等の摂動結果は一切生成されていない。tmux は自動終了済み（再起動していない）。
  - **この 3 ファイルも sha256 固定済みの事前登録ファイルであり、simulation-engineer の権限では修正しない。**
    stats セクションのときと同じ理由（設計側の変更ログ・チェックリスト・再ハッシュが必要）で、
    実行を止めて報告する。自己判断でファイルの形式を書き換えたり、run_isp.py 側を変えたりしていない。
  - リモートの複製（tokenized dataset 等）は `isp_remote.sh clean` で削除し、`ssh` での存在確認（`ls` が
    「そのようなファイルやディレクトリはありません」）で消えたことを確認した。
- **状態: 停止（work/state_rows-seed{31,32,33}.json の形式不一致待ち）。摂動は 0 遺伝子×0 seed。GPU 計算は実質ゼロ（6 秒でクラッシュ）。**

### 固定後の修正 2 件目（experiment-planner, 2026-09-27 23:1x）
- 実行担当の dry-run で 2 件目の不備。`perturbation.genes_from` / `genes_to_perturb: from_file` /
  `controls.negative_sampler.genes_from` という外部ファイル参照の書き方は `run_isp.py` が対応していなかった。
- **方式 (a)（設定に直接展開）を採用**。共有スクリプトは触っていない。`genes_004.json` / `null_pool_004.json` から
  機械生成して config.yaml に展開（転記ミスを避けるため手書きしていない）。役割の対応:
  `genes_to_perturb` = 陽性対照候補 10（run_isp の role は target）、`controls.positive` = ACTA2/TAGLN、
  `negative_fixed` + `fixed_role: housekeeping` = ACTB/GAPDH/PPIA、`extra_fixed` = 中立対照パネル 7、
  `negative_sampler{method: explicit, genes}` = 対照 84。
- 併せて `controls.negative_random_n: 84` を明示（explicit のとき dry-run の内蔵見積もりが対照遺伝子を数えるようになる。
  guardrails 2026-09-27 の「dry-run が explicit の null を数えない」は、この項目が未設定だったことが原因）、
  `compute.device` を `cuda` → `mps`（EXP-20260927-003 と同じ書き方。実行は isp_remote.sh が cuda で上書き）。
- **遺伝子の集合は不変を確認**: 役割ごとに件数と Ensembl ID の集合を照合 →
  positive_goal 10 / positive_myof 2 / housekeeping 3 / negctrl_panel 7 / negative 84 = 106 が完全一致、
  重複なし、回さない遺伝子（LTBR/RELB/NFKB2/CDKN1A/盲検 3/ラベル 3）の混入なし。
  判定基準・細胞選抜・状態の定義・検出力の見積もりは変更なし。
- **`--dry-run` が警告なしで通ることを確認**（forward_passes 63,600 = 106 遺伝子 × 200 細胞 × 3 seed）。
- plan.md §8 に dry-run の内蔵見積もりの扱いを追記、§12 変更ログに 1 行、§14 チェックリストに 3 行。
  `work/SHA256SUMS-preregistered.txt` を再固定し、**`shasum -c` で 17/17 OK**。固定後の dry-run でも 17/17 一致。
- **摂動は実行していない。GPU 時間の消費ゼロ。**

## 2026-09-27 23:05〜23:33 (simulation-engineer / 4 回目の再開 — smoke → 追加のセントロイド → 本番へ)
- Labretary から `work/state_rows-seed{31,32,33}.json` を `{"seed":..,"note":..,"rows":{...}}` の形に
  包み直した連絡。**自分でも独立に再確認**:
  - 実験ディレクトリ内で `shasum -a 256 -c work/SHA256SUMS-preregistered.txt` → 17/17 OK。
  - 3 ファイルとも `list(d.keys()) == ['seed','note','rows']` を確認し、`rows` の中身の細胞数が
    Labretary の報告どおり一致することを検算（FIB_org 614/600/603, FIB_neg_TLShigh 614/600/603,
    FIB_neg_TLSlow 677/663/655）。
  - `.venv/bin/python scripts/geneformer/run_isp.py --config experiments/EXP-20260927-004/config.yaml --dry-run`
    → `dry-run OK`（forward_passes 63600, 変わらず）。
- **smoke（本番前の小さな試し実行）を実施**。事前登録の config.yaml 本体は変更せず、複製に対してのみ
  上書きした一時 config を `/private/tmp/.../scratchpad/exp004-smoke/config.yaml` に作成
  （python で yaml を読み書き。上書き内容: `output.experiment_dir` を scratchpad 配下に、
  `robustness.seeds: [31]`、`perturbation.genes_to_perturb` を先頭 2 遺伝子(VCAM1, ICAM1)、
  `cell_selection.max_cells_per_gene: 10`(`max_cells_per_donor: 10`, `min_cells_per_gene: 5`)、
  対照は `negative_sampler.genes` の先頭 2 個のみに縮小、`positive`/`negative_fixed`/`extra_fixed` は空に）。
  - `--dry-run` で先に確認（`forward_passes: 40`, `dry-run OK`。ポジティブ対照 0 件・seed 1 個の WARN のみ）。
  - 本実行: `.venv/bin/python scripts/geneformer/run_isp.py --config <smoke config> `
    （TOKENIZERS_PARALLELISM=false）。**exit code 0、約 8 分半で完了。**
  - 確認できたこと:
    - **状態のセントロイドの抽出（今回まで 3 回落ちていた箇所）が通った**:
      `results/run-31-L-1/state_emb_vectors.npz` に `FIB_neg_TLShigh`/`FIB_org`/`FIB_neg_TLSlow`
      が 1152 次元ベクトルで保存されていることを確認。
    - 摂動が 4 遺伝子（target 2 + negative 2）とも最後まで通り、
      `results/run-31-L-1/genes/<ENSG>/emb/{orig,pert}.npy` が生成されることを確認。
    - `results/cells.csv`（41 行 = 40 細胞 + ヘッダ）、`results/summary.csv` を確認。
  - **smoke の出力（`/private/tmp/.../scratchpad/exp004-smoke/`, 約 84M）は確認後に `rm -rf` で削除済み。**
    本番の `experiments/EXP-20260927-004/results/` には触れていない。
- **smoke とは別に発見した欠落: `work/extra_centroids-seed{31,32,33}.npz`（`evaluate_004.py` が
  `m_primary`/`m_myof` などの計算で必須に読む）が未生成だった。** `config.yaml` の `extra_centroids`
  セクションは `run_isp.py` からは一切参照されない（grep で該当なし）。EXP-20260927-003 の
  `analysis/extra_centroids.py`（config を引数に取る汎用スクリプト。`select_genes_fib.py` の
  `balanced()`/`depth_match()`（`DEPTH_TOL=0.10`）を使うだけで、遺伝子選定など EXP-003 固有の判断は
  含まない）を EXP-20260927-004 の config.yaml に向けて実行し、3 seed 分を生成した
  （run_isp.py 同様、シミュレーション実行に使う既存インフラとして扱った。新規の分析コードは書いていない）:
  - `.venv/bin/python experiments/EXP-20260927-003/analysis/extra_centroids.py --config experiments/EXP-20260927-004/config.yaml --seed {31,32,33}`
  - **事故と修復**: seed 32 について、手動でバックグラウンド起動したプロセスと `run_in_background` の
    プロセスが二重に起動してしまった（同じ出力パスへ競合）。`ps aux` で検出し、後発の重複プロセスを
    `kill -9` した後、途中生成物（`work/extra-centroid-cells-seed32.dataset`,
    `work/extra-centroid-embs-seed32`）を削除してから単独プロセスでやり直した。
    最終生成物（`extra_centroids-seed32.{npz,json}`）は競合のなかった単独実行の出力のみ。
  - 結果（判定には使わない記述値。§7 頑健性チェック用）:
    | seed | MYOF_neg 細胞数/提供者数 | FIB_org_unmatched 細胞数/提供者数 | dim |
    |---|---|---|---|
    | 31 | 615 / 32 | 711 / 65 | 1152 |
    | 32 | 615 / 30 | 711 / 65 | 1152 |
    | 33 | 597 / 34 | 711 / 65 | 1152 |
  - 中間生成物 `work/extra-centroid-{cells,embs}-seed{31,32,33}*`（計 約80M）はキャッシュとして残す
    （EXP-003 と同じ扱い。人由来データそのものではなく、プロジェクト work/ 内のローカルキャッシュ）。

## 2026-09-27 23:34〜23:35 (simulation-engineer) — 本番再開
- リモート計算機 2 号へ config.yaml・run_isp.py・capture_env.sh・state_rows-seed{31,32,33}.json を
  sha256 照合つきで再送。tokenized.dataset とモデル重みは既存のもの（サイズ・mtime 一致）を再利用。
- リモートで --dry-run → forward_passes: 63600, dry-run OK。
- **23:34:53 に isp_remote.sh run を実行したが exit=127 で即失敗**（REMOTE_GF_HOME の export を
  同じ Bash 呼び出しに含め忘れ、GF_HOME が既定値 ~/workspace/Geneformer に戻っていたため
  .venv/bin/python が見つからなかった。simulation-engineer 側の操作ミスで、config やパイプラインの
  不備ではない。GPU 計算なし、0.1 秒で終了）。
- REMOTE_GF_HOME を同じコマンドに含めて再実行。**実開始 2026-09-27 23:35:40 JST（ログの === start は
  23:35:45）。** 20 秒後に status で RUNNING と === run-31-L-1 === を確認（今回は 6 秒で落ちていない）。
- 監視用のバックグラウンドループを起動（120 秒間隔で status 相当を ssh 経由で確認。
  seed 31→32 の切り替わりを検出した時点で経過が 2.0 h を超えていれば tmux kill-session で停止、
  全体が 6.0 h を超えたら同様に停止。ジョブが自然終了/クラッシュしても検出して抜ける）。
  開始 epoch（リモート基準）= 1790519745。

## 2026-09-28 01:37〜 (simulation-engineer) — 完了・回収・評価
- **ジョブは 23:35:40〜01:37:21 JST（実測 2.03 h）で exit=0（正常終了）。**
  seed 別の所要時間: seed31 23:35:45→00:16:23（0.68 h）、seed32 00:16:23→00:56:52（0.67 h）、
  seed33 00:56:52→01:37:20（0.67 h）。**中断規則（seed31 完了時点で 2.0 h 超／全体 6 h 超）には
  抵触していない。** 監視ループ（120 秒間隔で ssh 経由の tmux 生死確認）は log パスの変数展開ミス
  （ローカルの `~` 展開でリモートパスとして無効になっていた。`ALIVE` 判定自体は影響を受けず機能）で
  seed 単位の進捗行を拾えなかったが、ジョブ終了は正しく検出した。
- `isp_remote.sh pull` で回収（約 601 MB）。**ファイル数の照合**:
  - `results/summary.csv` 319 行（ヘッダ + 318 = 106 遺伝子 × 3 seed、計画と一致）
  - `results/cells.csv` 63,571 行（ヘッダ + 63,570 = 計画の延べ細胞×遺伝子と完全一致）
  - 3 seed とも `run-{31,32,33}-L-1/genes/` 配下に 106 遺伝子ディレクトリ、
    各 `emb/{orig,pert}.npy` が 106×3=318 件とも存在することを `find | wc -l` で確認
  - `cells.csv` の遺伝子集合と `blocked_genes`（LTBR/RELB/NFKB2/CDKN1A/TP53/SPP1/LCN2）の重なりを
    python で検算 → **重なり 0 件**。ユニーク遺伝子数 106。
- `isp_remote.sh clean` でリモートの複製（3.3G）を削除し、`ssh` での存在確認（ディレクトリ不存在）で
  消えたことを確認した。
- **評価**: `analysis/evaluate_004.py` を実行（変更なし・そのまま）。
  `results/eval/criteria_004.json`, `results/eval/gene_seed_004.csv` を生成。
  - **M1（較正）〜M4（方向性）: すべて合格。**
  - **M5（目標の軸での陽性対照, 最重要）: 合格。候補 10 のうち SYN3（データ層）の 1 つが
    3 seed すべてで q≥0.90 かつ D_raw>0 かつ 2.5 標準偏差以上を満たした
    （sd_units 3.19/3.46/3.62、q 1.0/1.0/1.0）。**
  - 偶然合格率の実測: 較正に使えた対照遺伝子 82 個（84 中 2 個は 3 seed 中どこかで
    使える対照遺伝子 <20 のため欠測、不合格には数えていない）で同じ M5 規則を当てはめたところ
    1/82 = 0.0122 が偶然合格。
  - 欠測（q が計算できなかった null 遺伝子×seed）: 4 件（ZBTB7C が 3 seed とも、MICALL2 が seed32 のみ）。
    判定対象 10 遺伝子・筋線維芽細胞陽性対照 2・恒常発現遺伝子 3・中立対照パネル 7 には欠測なし。
- **状態: 完了。** 出力: `experiments/EXP-20260927-004/results/`（summary.csv, cells.csv, run-{31,32,33}-L-1/,
  env.json, env-remote.json, device.json）、`experiments/EXP-20260927-004/results/eval/{criteria_004.json,
  gene_seed_004.csv}`。**判定の確定は PI（G3）。仮説の status は変更していない。**
