rule delly_call:
    input: 
        bam=f"{RESULTS}/alignment/{{sample}}_sorted.bam",
        bai=f"{RESULTS}/alignment/{{sample}}_sorted.bam.bai",
        ref=REF,
        fai=f"{REF}.fai"
    output: 
        bcf=f"{RESULTS}/variants/{{sample}}_svs.bcf"
    conda: ENV("svcalling.yaml")
    threads: config["threads"]["delly_call"]
    log: f"{LOGS}/variants/{{sample}}_delly_call.log"
    shell: """
        delly call \
        -g {input.ref:q} \
        -o {output.bcf:q} \
        {input.bam:q} 2> {log}
    """

rule bcf_to_vcf:
    input: f"{RESULTS}/variants/{{sample}}_svs.bcf"
    output: 
        vcf=f"{RESULTS}/variants/{{sample}}_svs.vcf.gz",
        tbi=f"{RESULTS}/variants/{{sample}}_svs.vcf.gz.tbi"
    conda: ENV("svcalling.yaml")
    log: f"{LOGS}/variants/{{sample}}_bcf_to_vcf.log"
    shell: 
        "bcftools convert -O z -o {output.vcf:q} {input:q} 2> {log} ; "
        "tabix -p vcf {output.vcf:q} 2>> {log}"

rule vcf_to_csv:
    input: 
        vcf=f"{RESULTS}/variants/{{sample}}_svs.vcf.gz",
        tbi=f"{RESULTS}/variants/{{sample}}_svs.vcf.gz.tbi"
    output: f"{RESULTS}/{{sample}}_structural_variants.csv"
    conda: ENV("svcalling.yaml")
    log: f"{LOGS}/variants/{{sample}}_vcf_to_csv.log"
    shell: "bash workflow/scripts/vcf_to_csv.sh {input.vcf:q} {output:q} {log}"