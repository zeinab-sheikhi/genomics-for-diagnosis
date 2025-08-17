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

rule vcf_to_bed:
    input:
        vcf=get_variant_outputs()["vcf"],
        tbi=get_variant_outputs()["vcf_index"]
    output: get_variant_outputs()["bed"]
    conda: get_env("wgs.yaml")
    shell: "bash workflow/scripts/vcf2bed.sh {input.vcf:q} {output:q}"

rule sv_gene_annotation:
    input:
        sv_bed=get_variant_outputs()["bed"],
        genes_bed=get_genes_bed()
    output: get_variant_outputs()["annotated"]
    conda: get_env("wgs.yaml")
    shell: "bash workflow/scripts/annotate.sh {input.sv_bed:q} {input.genes_bed:q} {output:q}"
