rule star_fusion:
    input:
        chim_junction=rules.star_align.output.chim_junction,
    output:
        abridged="star_fusion/{sample}/star-fusion.fusion_predictions.abridged.tsv",
        full="star_fusion/{sample}/star-fusion.fusion_predictions.tsv",
    log:
        ".log/star_fusion/{sample}.log",
    benchmark:
        ".log/star_fusion/{sample}.bm"
    conda:
        config["conda"]["star_fusion"]
    threads: config["threads"]["high"]
    params:
        ctat_lib=config["database"]["ctat_lib"],
    shell:
        """
        STAR-Fusion --genome_lib_dir {params.ctat_lib} \
            -J {input.chim_junction} \
            --output_dir star_fusion/{wildcards.sample} \
            --CPU {threads} \
            --FusionInspector validate \
            --examine_coding_effect \
            --denovo_reconstruct \
            &> {log}
        """
