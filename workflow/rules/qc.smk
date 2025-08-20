rule fastqc:
    input: 
        r1=get_fastq_r1(),
        r2=get_fastq_r2()
    output:
        r1_html=get_qc_outputs()["r1_html"],
        r2_html=get_qc_outputs()["r2_html"]
    conda: get_env("wgs.yaml")
    resources: 
        mem_mb=config.tools.fastqc.memory,  
    threads: config.tools.fastqc.threads
    shell: 
        "mkdir -p workflow/reports/ && "
        "fastqc {input.r1:q} {input.r2:q} -o workflow/reports/"