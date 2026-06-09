SAMPLE="INDEX-11"
FQ1="/data/mengxf/Project/KML260605_SH-RNA-FUSION/1M/INDEX-11_S12_R1_001.fastq.gz"
FQ2="/data/mengxf/Project/KML260605_SH-RNA-FUSION/1M/INDEX-11_S12_R2_001.fastq.gz"
THREADS=16

# trim fq
mamba -n basic2 run fastp \
    -w ${THREADS} \
    --qualified_quality_phred 20 \
    --length_required 36 \
    --cut_front \
    --cut_tail \
    --cut_mean_quality 20 \
    --detect_adapter_for_pe \
    --trim_poly_g \
    -i $FQ1 \
    -I $FQ2 \
    -o ${SAMPLE}.trim.1.fq.gz \
    -O ${SAMPLE}.trim.2.fq.gz \
    -h ${SAMPLE}.html \
    -j ${SAMPLE}.json

# align
# ! 参数有问题, snakemake 需要改, --chimNonchimScoreDropMin --chimOutType Junctions WithinBAM
mamba -n basic2 run STAR \
    --runThreadN ${THREADS} \
    --genomeDir "/data/mengxf/Database/STAR/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa.star.idx" \
    --readFilesIn ${SAMPLE}.trim.1.fq.gz ${SAMPLE}.trim.2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outBAMsortingThreadN ${THREADS} \
    --outFileNamePrefix "${SAMPLE}_" \
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
    --quantMode GeneCounts

samtools index INDEX-11_Aligned.sortedByCoord.out.bam

# fusion detection
mamba -n basic2 run STAR-Fusion \
    --genome_lib_dir /data/mengxf/Database/STAR/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir \
    -J INDEX-11_Chimeric.out.junction \
    --output_dir star_fusion \
    --CPU ${THREADS} \
    --FusionInspector validate \
    --examine_coding_effect \
    --denovo_reconstruct

mamba -n fusion run arriba \
    -x INDEX-11_Aligned.sortedByCoord.out.bam \
    -a /data/mengxf/Database/STAR/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_genome.fa \
    -g /data/mengxf/Database/STAR/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf \
    -o ${SAMPLE}_fusions.tsv \
    -O ${SAMPLE}_fusions_discarded.tsv \
    -b /home/mengxf/miniforge3/envs/fusion/var/lib/arriba/blacklist_hg38_GRCh38_v2.5.1.tsv.gz \
    -k /home/mengxf/miniforge3/envs/fusion/var/lib/arriba/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
    -t /home/mengxf/miniforge3/envs/fusion/var/lib/arriba/known_fusions_hg38_GRCh38_v2.5.1.tsv.gz \
    -p /home/mengxf/miniforge3/envs/fusion/var/lib/arriba/protein_domains_hg38_GRCh38_v2.5.1.gff3 \
    -T 5

# ! draw_fusions.R 不是 draw_fusions, 参数也要改
mamba -n fusion run draw_fusions.R \
    --annotation=/data/mengxf/Database/STAR/GRCh38_gencode_v44_CTAT_lib_Oct292023.plug-n-play/ctat_genome_lib_build_dir/ref_annot.gtf \
    --fusions=${SAMPLE}_fusions.tsv \
    --output=${SAMPLE}_fusions.pdf \
    --cytobands=/home/mengxf/miniforge3/envs/fusion/var/lib/arriba/cytobands_hg38_GRCh38_v2.5.1.tsv \
    --proteinDomains=/home/mengxf/miniforge3/envs/fusion/var/lib/arriba/protein_domains_hg38_GRCh38_v2.5.1.gff3
