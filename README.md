# Genomics Pipeline for Structural Variant Diagnostics

A Snakemake-based workflow for processing whole genome sequencing data to identify structural variants for clinical diagnostics. This pipeline takes patient FASTQ files and reference genome, then produces a clinical report of structural variants in CSV format.

## Installation
1. Download Miniconda3 installer and install conda
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh
```
2. Clone this repository:
   ```bash
   git clone <repository-url>
   cd genomics-for-diagnosis
   ```

3. Run the installation script and activate conda environment:
   ```bash
   chmod +x install.sh
   source install.sh
   ```

## Configuration

**Configure the workflow.**

-   **config files**:
    -   [`config.yaml`](/config/config.yaml) - main workflow configuration including sample names, file paths, tool settings, and resource allocation


-   **input files**:
    -   **Required data to be added by user**:
        - paired-end FASTQ files (place in `workflow/data/` directory):
          - `workflow/data/HG002_R1_wgs_chr21.fastq.gz` 
          - `workflow/data/HG002_R2_wgs_chr21.fastq.gz`
    -   **Provided reference data**:
        - reference genome: `workflow/data/fasta/chr21.fa`

-   **output files**:
    -   quality control reports in `workflow/reports/`
    -   **final clinical report**: `workflow/reports/HG002_structural_variants.csv`


**Note**: The example dataset uses chromosome 21 data from sample HG002. Place your paired-end FASTQ files in the `workflow/data/` directory. Update the file paths in `config/config.yaml` if using different sample names or file locations.
