rule bwa_align:
    input:
        r1=config["fastq"]["r1"],
        r2=config["fastq"]["r2"],
        ref=config["reference"]["fa"]
        bwt=config["reference"]["fa"] + ".bwt"
    output:
        temp("results/alignment/HG002_raw.bam")
    conda: "../env/alignment.yaml"
    threads: config["threads"]["bam_align"]
    shell: """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} | \
        samtools view -b -o {output} - 2>> {log}
    """

rule sort_bam:
    input:
    output:
    conda: 
    threads:
    log:
    shell: ""

rule index_bam:
    input:
    output:
    conda:
    log:
    shell: ""

