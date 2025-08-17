rule bwa_index:
    input: 
        get_fasta()
    output:
        get_bwa_idx()
    conda: get_env("wgs.yaml")
    shell:  "bwa index {input:q}"

rule samtools_faidx:
    input: 
        get_fasta()
    output:
        get_fasta_fai()
    conda: get_env("wgs.yaml")
    shell: "samtools faidx {input:q}"