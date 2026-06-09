rule fastqc_raw:
    input:
        rules.create_symlinks.output.fq1,
        rules.create_symlinks.output.fq2,
    output:
        directory("qc/fastqc/{sample}"),
    log:
        ".log/fastqc/{sample}.log",
    benchmark:
        ".log/fastqc/{sample}.bm"
    conda:
        config["conda"]["fastqc"]
    threads: config["threads"]["low"]
    shell:
        "mkdir -p {output} && fastqc -t {threads} {input} -o {output} &> {log}"
