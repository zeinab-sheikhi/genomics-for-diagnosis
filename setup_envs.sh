#!/bin/bash
# setup_envs.sh

echo "Creating genomics environments..."

envs=("qc" "alignment" "svcalling")

for env in "${envs[@]}"; do
    if conda env list | grep -q "^${env} "; then
        echo "Environment '$env' already exists - skipping"
    else
        echo "Creating environment: $env"
        conda env create -f "workflow/envs/${env}.yaml"
    fi
done

echo "All environments ready!"