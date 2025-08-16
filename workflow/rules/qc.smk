rule fastqc:
    input: 
        r1=config["fastq"]["r1"], 
        r2=config["fastq"]["r2"]
    output:
        "results/qc/HG002_R1_wgs_chr21_fastqc.html",    
        "results/qc/HG002_R2_wgs_chr21_fastqc.html"
    conda: "../envs/qc.yaml"
    shell: "fastqc {input.r1} {input.r2} -o results/qc/"