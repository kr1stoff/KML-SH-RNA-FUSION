rule star_align:
    input:
        fq1=rules.fastp.output.fq1,
        fq2=rules.fastp.output.fq2,
    output:
        bam="star/{sample}_Aligned.sortedByCoord.out.bam",
        chim_junction="star/{sample}_Chimeric.out.junction",
    log:
        ".log/star/{sample}.log",
    benchmark:
        ".log/star/{sample}.bm"
    conda:
        config["conda"]["star"]
    threads: config["threads"]["high"]
    params:
        genome=config["database"]["star_index"],
    shell:
        """
        STAR --runThreadN {threads} \
            --genomeDir {params.genome} \
            --readFilesIn {input.fq1} {input.fq2} \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD \
            --outBAMsortingThreadN {threads} \
            --outFileNamePrefix star/{wildcards.sample}_ \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 8 \
            --chimOutJunctionFormat 1 \
            --chimOutType Junctions WithinBAM \
            --chimNonchimScoreDropMin 10 \
            --chimMultimapScoreRange 3 \
            --chimMultimapNmax 20 \
            --chimScoreJunctionNonGTAG -4 \
            --chimScoreSeparation 1 \
            --alignSJstitchMismatchNmax 5 -1 5 5 \
            --alignIntronMax 500000 \
            --alignMatesGapMax 1000000 \
            --alignSJDBoverhangMin 3 \
            --outSJtype Standard \
            --twopassMode Basic \
            --peOverlapNbasesMin 12 \
            --peOverlapMMp 0.1 \
            --outFilterMultimapNmax 20 \
            --outFilterMismatchNmax 10 \
            --outFilterScoreMinOverLread 0 \
            --outFilterMatchNminOverLread 0 \
            --outFilterMatchNmin 0 \
            --quantMode GeneCounts \
            &> {log}
        """


rule samtools_index:
    input:
        rules.star_align.output.bam,
    output:
        "star/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        ".log/samtools/{sample}.index.log",
    benchmark:
        ".log/samtools/{sample}.bm"
    conda:
        config["conda"]["samtools"]
    threads: config["threads"]["low"]
    shell:
        "samtools index -@ {threads} {input} &> {log}"
