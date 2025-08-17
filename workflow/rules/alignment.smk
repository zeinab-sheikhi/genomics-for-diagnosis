rule bwa_align:
    input:
        r1=get_fastq_r1(),
        r2=get_fastq_r2(),
        ref=get_fasta(),
        bwt=get_bwt()
    output:
        temp(get_alignment_outputs()["raw_bam"])
    conda: get_env("wgs.yaml")
    threads: config.tools.bwa.threads
    shell: """
        mkdir -p workflow/data/bam && \
        bwa mem -t {threads} {input.ref:q} {input.r1:q} {input.r2:q} | \
        samtools view -b -o {output:q} - 
    """

rule sort_bam:
    input: get_alignment_outputs()["raw_bam"]
    output: get_alignment_outputs()["sorted_bam"]
    conda: get_env("wgs.yaml")
    threads: config.tools.samtools.threads
    shell: "samtools sort -@ {threads} -o {output:q} {input:q}"

rule index_bam:
    input: get_alignment_outputs()["sorted_bam"]
    output: get_alignment_outputs()["bam_index"]
    conda: get_env("wgs.yaml")
    threads: config.tools.samtools.threads
    shell: "samtools index {input:q} {output:q}"
