rule bwa_index:
    input:
        get_fasta(),
    output:
        get_bwa_idx(),
    resources:
        mem_mb=config.tools.bwa.memory,
    threads: config.tools.bwa.threads
    conda:
        get_env("wgs.yaml")
    shell:
        "bwa index {input:q}"


rule samtools_faidx:
    input:
        get_fasta(),
    output:
        get_fasta_idx(),
    resources:
        mem_mb=config.tools.samtools.memory,
    threads: config.tools.samtools.threads
    conda:
        get_env("wgs.yaml")
    shell:
        "samtools faidx {input:q}"
