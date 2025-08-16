rule bwa_index:
    input: REF
    output:
        expand("{fa}.{ext}", fa=REF, ext=["bwt", "sa", "pac", "ann", "amb"])
    conda: ENV("alignment.yaml")
    shell:  "bwa index {input:q}"

rule samtools_faidx:
    input: REF
    output: f"{REF}.fai"
    conda: ENV("alignment.yaml")
    log: f"{LOGS}/samtools_faidx.log"
    shell: "samtools faidx {input:q} 2> {log}"