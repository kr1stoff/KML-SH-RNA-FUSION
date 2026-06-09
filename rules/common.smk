rule create_symlinks:
    input:
        fq1=lambda w: samples_df.loc[samples_df["sample"] == w.sample, "fq1"].iloc[0],
        fq2=lambda w: samples_df.loc[samples_df["sample"] == w.sample, "fq2"].iloc[0],
    output:
        fq1="rawdata/{sample}_R1.fastq.gz",
        fq2="rawdata/{sample}_R2.fastq.gz",
    log:
        ".log/create_symlinks/{sample}.log",
    benchmark:
        ".log/create_symlinks/{sample}.bm"
    run:
        import os

        os.makedirs("rawdata", exist_ok=True)
        os.symlink(input.fq1, output.fq1)
        os.symlink(input.fq2, output.fq2)
