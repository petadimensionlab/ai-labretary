.PHONY: help ghostty ghostty-splits status night night-dry demo setup
# make setup で作られた .venv があればそれを使う（上書き: make demo PY=...）
PY ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)
help:            ## このヘルプ
	@grep -E '^[a-z-]+:.*##' Makefile | sed 's/:.*##/ —/'
ghostty:         ## Ghostty ダッシュボードを開く（ウィンドウ方式）
	bash ./labretary-ghostty.sh
ghostty-splits:  ## Ghostty ダッシュボード（1 ウィンドウ 4 分割, macOS）
	bash ./labretary-ghostty.sh --splits
status:          ## 研究状況
	bash scripts/lab_status.sh
night:           ## 夜間シフトを今すぐ実行
	bash ./labretary-nightshift.sh
night-dry:       ## 夜間シフトで実行されるプロンプトを表示
	DRY=1 bash ./labretary-nightshift.sh
demo:            ## GPU なしで検証・図のパイプラインを試す（合成データ）
	$(PY) scripts/geneformer/run_isp.py --config experiments/EXP-00000000-000-demo/config.yaml --mock
	$(PY) scripts/analysis/validate_isp.py --exp experiments/EXP-00000000-000-demo --results results-mock
	$(PY) scripts/analysis/plot_isp.py --exp experiments/EXP-00000000-000-demo --results results-mock
setup:           ## 既存の Geneformer (../Geneformer) を接続
	bash scripts/setup_geneformer.sh
