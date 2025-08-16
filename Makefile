.PHONY: setup-envs

setup-envs:
	@echo "Setting up environments..."
	@workflow/scripts/setup_envs.sh

help:
	@echo "setup-envs - Create all conda environments"
	@echo "run        - Run the genomic pipeline for diagnosis"