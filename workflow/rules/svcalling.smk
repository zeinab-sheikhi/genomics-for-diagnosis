rule delly_call:
    input: 
        bam="workflow/results/alignment/HG002_sorted.bam",
        bai="workflow/results/alignment/HG002_sorted.bam.bai",
        ref=config["reference"]["fa"],
        fai=config["reference"]["fa"] + ".fai"
    output: 
        bcf="workflow/results/variants/HG002_svs.bcf"
    conda: "../envs/svcalling.yaml"
    threads: config["threads"]["delly_call"]
    log: "logs/delly_call.log"
    shell: """
        delly call \
        -g {input.ref} \
        -o {output.bcf} \
        {input.bam} 2> {log}
    """

rule bcf_to_vcf:
    input: "workflow/results/variants/HG002_svs.bcf"
    output: 
        vcf="workflow/results/variants/HG002_svs.vcf.gz",
        tbi="workflow/results/variants/HG002_svs.vcf.gz.tbi"
    conda: "../envs/svcalling.yaml"
    log: "logs/bcf_to_vcf.log"
    shell: """
        bcftools convert -O z -o {output.vcf} {input} 2> {log}
        tabix -p vcf {output.vcf} 2>> {log}
    """

rule vcf_to_csv:
    input: 
        vcf="workflow/results/variants/HG002_svs.vcf.gz",
        tbi="workflow/results/variants/HG002_svs.vcf.gz.tbi"
    output: "workflow/results/structural_variants.csv"
    conda: "../envs/svcalling.yaml"
    log: "logs/vcf_to_csv.log"
    shell: "bash workflow/scripts/vcf_to_csv.sh {input.vcf} {output} {log}"