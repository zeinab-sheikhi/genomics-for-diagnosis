"""YAML validator for config.yaml file."""
import yaml
import sys
from typing import List, Optional
from pydantic import BaseModel, Field, field_validator


class FileExtensions(BaseModel):
    """File extensions configuration."""
    fasta: str = ".fasta"
    bwa_idx: List[str] = [".bwt", ".amb", ".ann", ".pac", ".sa"]
    fasta_idx: str = ".fai"
    bam: str = ".bam"
    bam_idx: str = ".bam.bai"
    bcf: str = ".bcf"
    vcf: str = ".vcf"


class Reference(BaseModel):
    """Reference genome configuration."""
    fa: str = Field(..., description="Path to reference FASTA file")
    bed: str = Field(..., description="Path to reference BED file")
    
    @field_validator('fa')
    @classmethod
    def validate_fasta_path(cls, v):
        if not v.endswith(('.fasta', '.fa', '.fna')):
            raise ValueError('Reference file must have .fasta, .fa, or .fna extension')
        return v.strip()
    
    @field_validator('bed')
    @classmethod
    def validate_bed_path(cls, v):
        if not v.endswith(('.bed')):
            raise ValueError('Reference BED file must have .bed extension')
        return v.strip()


class FastqFiles(BaseModel):
    """FASTQ input files configuration."""
    r1: str = Field(..., description="Path to R1 FASTQ file")
    r2: str = Field(..., description="Path to R2 FASTQ file")
    
    @field_validator('r1', 'r2')
    @classmethod
    def validate_fastq_extension(cls, v):
        if not v.endswith(('.fastq', '.fastq.gz', '.fq', '.fq.gz')):
            raise ValueError('FASTQ files must have .fastq, .fastq.gz, .fq, or .fq.gz extension')
        return v


class ToolResource(BaseModel):
    """Base tool resource configuration."""
    threads: int = Field(..., ge=1, description="Number of threads")
    memory: Optional[int] = Field(None, ge=1, description="Memory in MB")
    outdir: Optional[str] = Field(None, description="Output directory")


class BWA(ToolResource):
    """BWA alignment configuration."""
    pass


class Samtools(ToolResource):
    """Samtools configuration."""
    pass


class FastQC(ToolResource):
    """FastQC configuration."""
    pass


class PostProc(ToolResource):
    """SV post-processing configuration."""
    pass


class Delly(ToolResource):
    """Delly structural variant caller configuration."""
    sv_types: List[str] = Field(
        default=["BND", "DEL", "DUP", "INS", "INV"],
        description="Structural variant types to call"
    )
    
    @field_validator('sv_types')
    @classmethod
    def validate_sv_types(cls, v):
        valid_types = {"BND", "DEL", "DUP", "INS", "INV"}
        for sv_type in v:
            if sv_type not in valid_types:
                raise ValueError(f"Invalid SV type: {sv_type}. Must be one of {valid_types}")
        return v


class Tools(BaseModel):
    """All tools configuration."""
    bwa: BWA
    samtools: Samtools
    fastqc: FastQC
    delly: Delly
    postproc: PostProc


class FinalOutput(BaseModel):
    """Final output configuration."""
    filename: str = Field(..., description="Output filename")
    columns: List[str] = Field(..., description="Output CSV columns")
    results_dir: str = Field(default="workflow/results", description="Results directory")
    
    @field_validator('filename')
    @classmethod
    def validate_filename(cls, v):
        if not v.endswith('.csv'):
            raise ValueError('Output filename must end with .csv')
        return v


class GenomicsConfig(BaseModel):
    """Main genomics workflow configuration."""
    sample: str = Field(..., description="Sample identifier")
    reference: Reference
    fastq: FastqFiles
    file_exts: FileExtensions
    tools: Tools
    final_output: FinalOutput
    
    @field_validator('sample')
    @classmethod
    def validate_sample_name(cls, v):
        if not v or not v.strip():
            raise ValueError('Sample name cannot be empty')
        return v.strip()


def load_and_validate_config(config_file: str) -> GenomicsConfig:
    """Load and validate configuration file."""
    try:
        with open(config_file, 'r') as f:
            config_data = yaml.safe_load(f)
        
        # Validate using Pydantic
        config = GenomicsConfig(**config_data)
        return config
        
    except FileNotFoundError:
        raise FileNotFoundError(f"Configuration file not found: {config_file}")
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML syntax: {e}")
    except Exception as e:
        raise ValueError(f"Configuration validation failed: {e}")


def main():
    """Main validation function."""
    if len(sys.argv) != 2:
        print("Usage: python yaml_validator.py <config.yaml>")
        sys.exit(1)
    
    config_file = sys.argv[1]
    
    try:
        _ = load_and_validate_config(config_file)
        print("✅ Configuration is valid!")
        sys.exit(0)
        
    except Exception as e:
        print("❌ Configuration validation failed!")
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
