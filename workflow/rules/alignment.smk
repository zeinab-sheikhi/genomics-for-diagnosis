rule bwa_align:
    input:
        r1=r1,
        r2=r2,
        ref=REF,
        bwt=f"{REF}.bwt"
    output:
        temp(f"{RESULTS}/alignment/{{sample}}_raw.bam")
    conda: ENV("alignment.yaml")
    threads: config["threads"]["bwa_align"]
    log: f"{LOGS}/alignment/{{sample}}_bwa_align.log"
    shell: """
        bwa mem -t {threads} {input.ref:q} {input.r1:q} {input.r2:q} 2> {log} | \
        samtools view -b -o {output:q} - 2>> {log}
    """

rule sort_bam:
    input: f"{RESULTS}/alignment/{{sample}}_raw.bam"
    output: f"{RESULTS}/alignment/{{sample}}_sorted.bam"
    conda: ENV("alignment.yaml")
    threads: config["threads"]["sort_bam"]
    log: f"{LOGS}/alignment/{{sample}}_sort_bam.log"
    shell: "samtools sort -@ {threads} -o {output:q} {input:q} 2> {log}"

rule index_bam:
    input: f"{RESULTS}/alignment/{{sample}}_sorted.bam"
    output: f"{RESULTS}/alignment/{{sample}}_sorted.bam.bai"
    conda: ENV("alignment.yaml")
    threads: config["threads"]["index_bam"]
    log: f"{LOGS}/alignment/{{sample}}_index_bam.log"
    shell: "samtools index {input:q} {output:q} 2> {log}"
