import pandas as pd
from snakemake.utils import validate


configfile: workflow.source_path("config.yaml")


validate(config, "config.schema.yaml")

shell.executable("/bin/bash")
shell.prefix("set +eu; ")

col_names = ["sample", "fq1", "fq2"]
samples_df = pd.read_table(config["samples_tsv"], header=None, names=col_names)
samples = samples_df["sample"].tolist()


rule all:
    input:
        # todo: star-fusion 无结果会报错, 后面按需再调试
        # expand("annotation/{sample}_filtered_fusions.tsv", sample=samples),
        # expand("final/{sample}_fusion_report.tsv", sample=samples),
        # expand("final/{sample}_fusion_report.xlsx", sample=samples),
        # expand("final/{sample}_priority_fusions.tsv", sample=samples),
        "qc/multiqc",
        expand("star_fusion/{sample}/star-fusion.fusion_predictions.abridged.tsv", sample=samples),
        expand("arriba/{sample}_fusions.pdf", sample=samples),


include: "rules/common.smk"
include: "rules/fastqc.smk"
include: "rules/fastp.smk"
include: "rules/star_align.smk"
include: "rules/star_fusion.smk"
include: "rules/arriba.smk"
include: "rules/integrate_annotate.smk"
include: "rules/multiqc.smk"
