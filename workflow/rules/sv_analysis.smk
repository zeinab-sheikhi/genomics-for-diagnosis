rule sv_analysis:
    input: 
        vcf=get_sv_inputs()["vcf"],
        tbi=get_sv_inputs()["vcf_index"],
        csv=get_sv_inputs()["csv"]
    output:
        genome_plot=get_analysis_outputs()["genome_plot"],
        size_dist=get_analysis_outputs()["size_dist"],
        quality_plot=get_analysis_outputs()["quality_plot"],
    conda: get_env("postproc.yaml")
    params:
        sample=config.sample,
        outdir=config.tools.postproc.outdir
    threads: config.tools.postproc.threads
    shell: """
        mkdir -p {params.outdir} && \
        python workflow/scripts/sv_analysis.py \
            --vcf {input.vcf:q} \
            --csv {input.csv:q} \
            --sample {params.sample} \
            --outdir {params.outdir}
    """
