rule multiqc:
    input:
        expand("qc/fastqc/{sample}", sample=samples),
    output:
        directory("qc/multiqc"),
    log:
        ".log/multiqc.log",
    benchmark:
        ".log/multiqc.bm"
    conda:
        config["conda"]["multiqc"]
    shell:
        "multiqc {input} --outdir {output} 2> {log}"
