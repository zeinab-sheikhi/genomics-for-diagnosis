#!/bin/bash
# setup_envs.sh
echo "Creating genomics environments..."
conda env create -f workflow/envs/qc.yaml
conda env create -f workflow/envs/alignment.yaml  
conda env create -f workflow/envs/svcalling.yaml
conda env create -f workflow/envs/annotation.yaml
conda env create -f workflow/envs/processing.yaml
echo "All environments ready!"