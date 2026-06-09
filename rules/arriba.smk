rule arriba:
    input:
        bam=rules.samtools_index.input,
        bai=rules.samtools_index.output,
        ref=config["database"]["reference"],
        gtf=config["database"]["gtf"],
        blacklist=config["database"]["arriba_blacklist"],
        known_fusions=config["database"]["arriba_known_fusions"],
        protein_domains=config["database"]["arriba_protein_domains"],
    output:
        fusions="arriba/{sample}_fusions.tsv",
        discarded="arriba/{sample}_fusions_discarded.tsv",
    log:
        ".log/arriba/{sample}.log",
    benchmark:
        ".log/arriba/{sample}.bm"
    conda:
        config["conda"]["arriba"]
    threads: config["threads"]["medium"]
    shell:
        """
        arriba -x {input.bam} \
            -a {input.ref} \
            -g {input.gtf} \
            -o {output.fusions} \
            -O {output.discarded} \
            -b {input.blacklist} \
            -k {input.known_fusions} \
            -t {input.known_fusions} \
            -p {input.protein_domains} \
            -T 5 &> {log}
        """


rule draw_fusions:
    input:
        fusions=rules.arriba.output.fusions,
        gtf=config["database"]["gtf"],
        cytobands=config["database"]["arriba_cytobands"],
        protein_domains=config["database"]["arriba_protein_domains"],
    output:
        pdf="arriba/{sample}_fusions.pdf",
    log:
        ".log/arriba/{sample}.draw_fusions.log",
    benchmark:
        ".log/arriba/{sample}.draw_fusions.bm"
    conda:
        config["conda"]["arriba"]
    shell:
        """
        draw_fusions.R --annotation={input.gtf} \
            --fusions={input.fusions} \
            --output={output.pdf} \
            --cytobands={input.cytobands} \
            --proteinDomains={input.protein_domains} &> {log}
        """
