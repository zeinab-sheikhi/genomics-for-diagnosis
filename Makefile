.PHONY: setup-envs dry-run run generate-dag-graph generate-rule-graph generate-report

setup-envs:
	@echo "Setting up environments..."
	@chmod +x install.sh
	@./install.sh

dry-run:
	@snakemake -s workflow/Snakefile --dryrun

run:
	@snakemake -s workflow/Snakefile --use-conda --cores 4

generate-dag-graph:
	@snakemake -s workflow/Snakefile --dag | dot -Tpng > workflow/reports/dag.png

generate-rule-graph:
	@snakemake -s workflow/Snakefile --rulegraph | dot -Tpng > workflow/reports/rulegraph.png

generate-report:
	@snakemake -s workflow/Snakefile --report workflow/reports/report.html

help:
	@echo "setup-envs - Create all conda environments"
	@echo "dry-run    - Run the pipeline with dry run"
	@echo "run        - Run the genomic pipeline for diagnosis"
	@echo "generate-dag-graph - Generate DAG graph"
	@echo "generate-rule-graph - Generate rule graph(relationships between rules)"
	@echo "generate-report - Generate HTML report"