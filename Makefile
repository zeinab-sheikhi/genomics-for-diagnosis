.PHONY: setup-envs

setup-envs:
	@echo "Setting up environments..."
	@setup_envs.sh

dry-run:
	@snakemake --dryrun 

run:
	@snakemake --use-conda --cores 4

help:
	@echo "setup-envs - Create all conda environments"
	@echo "dry-run    - Run the pipeline with dry run"
	@echo "run        - Run the genomic pipeline for diagnosis"