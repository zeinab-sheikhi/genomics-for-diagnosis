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


def get_fasta_fai(check_exists: bool = False) -> str:
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


def make_all_outputs() -> list[str]:
    """Make all output paths for the workflow."""
    outfiles = []
    sample = config.sample

    outfiles.extend(get_bwa_idx(check_exists=False))  # BWA index files
    outfiles.append(get_fasta_fai(check_exists=False))  # FASTA index

    return outfiles
