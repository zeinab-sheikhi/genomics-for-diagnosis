rule fastqc:
    input: r1=r1, r2=r2
    output:
        R1_html = f"{RESULTS}/qc/{{sample}}_R1_wgs_chr21_fastqc.html",
        R2_html = f"{RESULTS}/qc/{{sample}}_R2_wgs_chr21_fastqc.html"
    conda: ENV("qc.yaml")
    shell: "fastqc {input.r1:q} {input.r2:q} -o {RESULTS}/qc/"