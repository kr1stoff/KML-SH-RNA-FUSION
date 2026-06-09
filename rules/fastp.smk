rule fastp:
    input:
        fq1=rules.create_symlinks.output.fq1,
        fq2=rules.create_symlinks.output.fq2,
    output:
        fq1="trim/{sample}_1.fastq.gz",
        fq2="trim/{sample}_2.fastq.gz",
        html="trim/{sample}.html",
        json="trim/{sample}.json",
    log:
        ".log/fastp/{sample}.log",
    benchmark:
        ".log/fastp/{sample}.bm"
    conda:
        config["conda"]["fastp"]
    threads: config["threads"]["medium"]
    params:
        extra=config["fastp"]["extra"],
    shell:
        """
        fastp -w {threads} {params.extra} \
            -i {input.fq1} -I {input.fq2} \
            -o {output.fq1} -O {output.fq2} \
            -h {output.html} -j {output.json} 2> {log}
        """
