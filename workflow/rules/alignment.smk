rule bwa_align:
    input:
        r1=config["fastq"]["r1"],
        r2=config["fastq"]["r2"],
        ref=config["reference"]["fa"],
        bwt=config["reference"]["fa"] + ".bwt"
    output:
        temp("workflow/results/alignment/HG002_raw.bam")
    conda: "../envs/alignment.yaml"
    threads: config["threads"]["bwa_align"]
    log: "logs/alignment/bwa_align.log"
    shell: """
        bwa mem -t {threads} {input.ref} {input.r1} {input.r2} 2> {log} | \
        samtools view -b -o {output} - 2>> {log}
    """

rule sort_bam:
    input: "workflow/results/alignment/HG002_raw.bam"
    output: "workflow/results/alignment/HG002_sorted.bam"
    conda: "../envs/alignment.yaml"
    threads: config["threads"]["sort_bam"]
    log: "logs/alignment/sort_bam.log"
    shell: "samtools sort -@ {threads} -o {output} {input} 2> {log}"

rule index_bam:
    input: "workflow/results/alignment/HG002_sorted.bam"
    output: "workflow/results/alignment/HG002_sorted.bam.bai"
    conda: "../envs/alignment.yaml"
    threads: config["threads"]["index_bam"]
    log: "logs/alignment/index_bam.log"
    shell: "samtools index {input} {output} 2> {log}"
