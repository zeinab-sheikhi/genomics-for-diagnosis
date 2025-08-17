rule delly_call:
    input: 
        bam=get_alignment_outputs()["sorted_bam"],
        bai=get_alignment_outputs()["bam_index"],
        ref=get_fasta(),
        fai=get_fasta_idx()
    output: 
        bcf=get_variant_outputs()["bcf"]
    conda: get_env("wgs.yaml")
    threads: config.tools.delly.threads
    shell: """
        mkdir -p {config.tools.delly.outdir} && \
        delly call \
        -g {input.ref:q} \
        -o {output.bcf:q} \
        {input.bam:q}
    """

rule bcf_to_vcf:
    input: get_variant_outputs()["bcf"]
    output: 
        vcf=get_variant_outputs()["vcf"],
        tbi=get_variant_outputs()["vcf_index"]
    conda: get_env("wgs.yaml")
    shell: 
        "bcftools convert -O z -o {output.vcf:q} {input:q} && "
        "tabix -p vcf {output.vcf:q}"

rule vcf_to_csv:
    input: 
        vcf=get_variant_outputs()["vcf"],
        tbi=get_variant_outputs()["vcf_index"]
    output: get_variant_outputs()["csv"]
    conda: get_env("wgs.yaml")
    shell: "bash workflow/scripts/vcf2csv.sh {input.vcf:q} {output:q}"