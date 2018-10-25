.DEFAULT_GOAL := analysis

requirements:
	pip install -r requirements.txt

process:
	looper run metadata/project_config.yaml

summarize:
	looper summarize metadata/project_config.yaml

mklog:
	mkdir -p log

merge:
	merge_signal -j -a cell_line,ip,condition metadata/project_config.yaml

peaks:
	call_peaks -j -c metadata/comparison_table.csv metadata/project_config.yaml

chip_analysis: summarize peaks
	ngs_analysis --data-type ChIP-seq metadata/project_config.yaml
chip_analysis_job: summarize peaks mklog
	sbatch -p longq --time 8-00:00:00 -c 12 --mem 80000 -J mthfd1 -o log/$(shell date +"%Y%m%d-%H%M%S").mthfd1.chipseq.log \
	--wrap "ngs_analysis --data-type ChIP-seq metadata/project_config.yaml"

analysis: merge
	python -u src/call_peaks.py
	sh src/peak_overlap.sh
	python -u src/peak_overlap.plot.py

	sh src/make_bigwigs.sh
	sh src/heatmaps.sh

	Rscript src/diffBind.R
	python -u src/diffBind.plot.py
	python -u src/diffBind.great_enrichments.plot.py


all: requirements process merge peaks analysis

.PHONY: requirements process summarize mklog analysis all
