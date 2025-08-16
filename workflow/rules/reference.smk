rule bwa_index:
    input: config["reference"]["fa"]
    output:
        config["reference"]["fa"] + ".bwt",
        config["reference"]["fa"] + ".sa",
        config["reference"]["fa"] + ".pac",
        config["reference"]["fa"] + ".ann",
        config["reference"]["fa"] + ".amb"
    conda: "../envs/alignment.yaml"
    shell:  "bwa index {input}"

rule samtools_faidx:
    input: config["reference"]["fa"]
    output: config["reference"]["fa"] + ".fai"
    conda: "../envs/alignment.yaml"
    log: "logs/samtools_faidx.log"
    shell: "samtools faidx {input} 2> {log}"