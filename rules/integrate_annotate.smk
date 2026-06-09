rule integrate_fusions:
    input:
        sf_abridged=rules.star_fusion.output.abridged,
        arr_fusions=rules.arriba.output.fusions,
    output:
        filtered="annotation/{sample}_filtered_fusions.tsv",
        all_merged="annotation/{sample}_all_merged_fusions.tsv",
    log:
        ".log/integrate/{sample}.log",
    benchmark:
        ".log/integrate/{sample}.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/integrate_fusions.py"


rule annotate_fusions:
    input:
        filtered=rules.integrate_fusions.output.filtered,
    output:
        report_tsv="final/{sample}_fusion_report.tsv",
        report_xlsx="final/{sample}_fusion_report.xlsx",
        priority="final/{sample}_priority_fusions.tsv",
    log:
        ".log/annotate/{sample}.log",
    benchmark:
        ".log/annotate/{sample}.bm"
    conda:
        config["conda"]["python"]
    script:
        "../scripts/annotate_fusions.py"
