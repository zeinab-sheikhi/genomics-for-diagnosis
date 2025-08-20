#!/usr/bin/env bash

set -xe
MY_ENV=genomics

eval "$(conda shell.bash hook)"
printenv

conda install -n base -c conda-forge mamba
conda activate base 
conda list

mamba env create -n $MY_ENV -f environment.yaml
conda activate $MY_ENV
conda list
