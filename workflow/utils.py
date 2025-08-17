from pathlib import Path

from validator import load_and_validate_config


config_path = Path(__file__).parent.parent / "config" / "config.yaml"
config = load_and_validate_config(config_path)


def get_env(env_name: str) -> str:
    """Get conda environment file path."""
    env_path = Path(__file__).parent / "envs" / env_name
    if not env_path.exists():
        raise FileNotFoundError(f"Conda environment file not found: {env_path}")
    return str(env_path)


def get_fasta() -> str:
    """Get reference FASTA file path with comprehensive validation."""
    fasta_path = Path(config.reference.fa)
    
    if not fasta_path.exists():
        raise FileNotFoundError(f"Reference FASTA file not found: {fasta_path}")
    
    valid_extensions = ('.fa', '.fasta', '.fna')
    if not fasta_path.suffix.lower() in valid_extensions:
        raise ValueError(f"Invalid FASTA file extension. Expected {valid_extensions}, got: {fasta_path.suffix}")
    
    if fasta_path.stat().st_size == 0:
        raise ValueError(f"Reference FASTA file is empty: {fasta_path}")
    
    return str(fasta_path)    


def get_bwa_idx(check_exists: bool = False) -> list[str]:
    """Get BWA index file paths.
    
    Args:
        check_exists: If True, validates that files exist (for inputs).
                     If False, just returns paths (for outputs).
    """
    fasta_path = get_fasta()
    
    bwa_extensions = config.file_exts.bwa_idx
    index_files = [fasta_path + ext for ext in bwa_extensions]
    
    if check_exists:
        for index_file in index_files:
            if not Path(index_file).exists():
                raise FileNotFoundError(f"BWA index file not found: {index_file}")
    
    return index_files


def get_fasta_idx(check_exists: bool = False) -> str:
    """Get samtools faidx output file path (.fai).
    
    Args:
        check_exists: If True, validates that file exists (for inputs).
                     If False, just returns path (for outputs).
    """
    fasta_path = get_fasta()
    
    fai_ext = config.file_exts.fasta_idx
    fai_file = fasta_path + fai_ext
    
    if check_exists:
        if not Path(fai_file).exists():
            raise FileNotFoundError(f"FASTA index file not found: {fai_file}")
    
    return fai_file


def get_fastq_r1(check_exists: bool = True) -> str:
    """Get R1 FASTQ file path with validation."""
    r1_path = Path(config.fastq.r1)
    
    if check_exists:
        if not r1_path.exists():
            raise FileNotFoundError(f"R1 FASTQ file not found: {r1_path}")
        
        if not str(r1_path).endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
            raise ValueError(f"Invalid R1 FASTQ file extension: {r1_path.suffix}")
        
        if r1_path.stat().st_size == 0:
            raise ValueError(f"R1 FASTQ file is empty: {r1_path}")
    
    return str(r1_path)


def get_fastq_r2(check_exists: bool = True) -> str:
    """Get R2 FASTQ file path with validation."""
    r2_path = Path(config.fastq.r2)
    
    if check_exists:
        if not r2_path.exists():
            raise FileNotFoundError(f"R2 FASTQ file not found: {r2_path}")
        
        if not str(r2_path).endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
            raise ValueError(f"Invalid R2 FASTQ file extension: {r2_path.suffix}")
        
        if r2_path.stat().st_size == 0:
            raise ValueError(f"R2 FASTQ file is empty: {r2_path}")
    
    return str(r2_path)


def get_qc_outputs(sample: str | None = None) -> dict:
    """Get FastQC output file paths."""
    if sample is None:
        sample = config.sample
    
    qc_outdir = config.tools.fastqc.outdir
    
    return {
        "r1_html": f"{qc_outdir}/{sample}_R1_wgs_chr21_fastqc.html",
        "r2_html": f"{qc_outdir}/{sample}_R2_wgs_chr21_fastqc.html",
    }


def get_alignment_outputs(sample: str | None = None, check_exists: bool = False) -> dict:
    """Get alignment output file paths."""
    if sample is None:
        sample = config.sample
    
    paths = {
        "raw_bam": f"workflow/data/bam/{sample}_raw.bam",
        "sorted_bam": f"workflow/data/bam/{sample}_sorted.bam", 
        "bam_index": f"workflow/data/bam/{sample}_sorted.bam.bai"
    }
    if check_exists:
        if not Path(paths["sorted_bam"]).exists():
            raise FileNotFoundError(f"Sorted BAM file not found: {paths['sorted_bam']}")
        if not Path(paths["bam_index"]).exists():
            raise FileNotFoundError(f"BAM index file not found: {paths['bam_index']}")
    
    return paths
    

def get_bwt() -> str:
    """Get BWA .bwt index file."""
    fasta_path = get_fasta()
    return fasta_path + ".bwt"


def get_variant_outputs(sample: str | None = None, check_exists: bool = False) -> dict:
    """Get structural variant calling output file paths."""
    if sample is None:
        sample = config.sample
    
    delly_outdir = config.tools.delly.outdir
    
    paths = {
        "bcf": f"{delly_outdir}/{sample}_svs.bcf",
        "vcf": f"{delly_outdir}/{sample}_svs.vcf.gz", 
        "vcf_index": f"{delly_outdir}/{sample}_svs.vcf.gz.tbi",
        "csv": f"workflow/reports/{sample}_structural_variants.csv" 
    }
    
    if check_exists:
        for path in paths.values():
            if not Path(path).exists():
                raise FileNotFoundError(f"Variant file not found: {path}")
    
    return paths


def make_all_outputs() -> list[str]:
    """Make all output paths for the workflow."""
    outfiles = []
    sample = config.sample

    # Reference files
    outfiles.extend(get_bwa_idx(check_exists=False))  # BWA index files
    outfiles.append(get_fasta_idx(check_exists=False))  # FASTA index

    # QC files
    qc_outputs = get_qc_outputs(sample)
    outfiles.extend(qc_outputs.values())

    # Alignment files 
    alignment_outputs = get_alignment_outputs(sample, check_exists=False)
    outfiles.extend([
        alignment_outputs["sorted_bam"],
        alignment_outputs["bam_index"]
    ])

    # Final variant output (CSV file)
    variant_outputs = get_variant_outputs(sample, check_exists=False)
    outfiles.append(variant_outputs["csv"])

    return outfiles
