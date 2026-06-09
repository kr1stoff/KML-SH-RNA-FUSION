# KML-SH-RNA-FUSION

RNA-seq 融合基因分析流程。FASTQ → QC → STAR-Fusion + Arriba → 注释 & 整合。

## 目录结构

```text
KML-SH-RNA-FUSION/
├── Snakefile                    # 主入口
├── config.yaml                  # 配置
├── config.schema.yaml           # schema 校验
├── rules/
│   ├── common.smk               # create_symlinks
│   ├── fastqc.smk               # FastQC 原始数据质控
│   ├── fastp.smk                # fastp 质控修剪
│   ├── star_align.smk           # STAR 比对 + samtools index
│   ├── star_fusion.smk          # STAR-Fusion 融合检测
│   ├── arriba.smk               # Arriba 融合检测 + draw_fusions 可视化
│   ├── integrate_annotate.smk   # 整合 + 注释规则
│   └── multiqc.smk              # MultiQC 汇总
├── scripts/
│   ├── integrate_fusions.py     # STAR-Fusion + Arriba 合并过滤
│   └── annotate_fusions.py      # 临床注释 + 优先级评分
└── tests/
    └── example.tsv              # 样本输入模板
```

## 流程

```text
rawdata (create_symlinks)
  └─ fastqc
     └─ fastp (trimming)
        └─ STAR (alignment + chimeric detection)
           ├─ STAR-Fusion (fusion detection)
           ├─ Arriba (fusion detection + PDF visualization)
           └─ samtools index
              ├─ integrate_fusions (merge ＋ filter)
              │  └─ annotate_fusions (clinical annotation ＋ scoring)
              └─ multiqc (aggregate QC report)
```

| 步骤 | 规则 | 输入 | 输出 |
|------|------|------|------|
| 1 | `create_symlinks` | input.tsv (fq1, fq2) | rawdata/*.fastq.gz |
| 2 | `fastqc_raw` | raw FQ | qc/fastqc/{sample}/ |
| 3 | `fastp` | raw FQ | trim/{sample}_{1,2}.fastq.gz |
| 4 | `star_align` | trimmed FQ | star/{sample}_Aligned.bam + Chimeric.junction |
| 5 | `samtools_index` | STAR BAM | star/{sample}.bam.bai |
| 6 | `star_fusion` | chimeric junction | star_fusion/{sample}/fusion_predictions*.tsv |
| 7 | `arriba` + `draw_fusions` | STAR BAM | arriba/{sample}_fusions.tsv + .pdf |
| 8 | `integrate_fusions` | SF + Arriba results | annotation/{sample}_*_fusions.tsv |
| 9 | `annotate_fusions` | filtered fusions | final/{sample}_fusion_report.{tsv,xlsx} + priority_fusions.tsv |
| 10 | `multiqc` | fastqc + STAR log | qc/multiqc/{sample}_multiqc.html |

## 使用

```bash
# 1. 修改 config.yaml 中的样本表和数据库路径
# 2. 运行
snakemake --cores 32 --use-conda --rerun-incomplete --scheduler greedy \
   --config samples_tsv=$PWD/tests/example.tsv \
   --directory /path/to/output/directory
```

### 输入格式

`tests/example.tsv`（三列，制表符分隔）：

```text
sample fq1 fq2
SAMPLE_NAME /path/to/SAMPLE_NAME_R1.fastq.gz /path/to/SAMPLE_NAME_R2.fastq.gz
```

## 配置

主要参数在 `config.yaml`，参见 `config.schema.yaml` 校验规则。

### 数据库路径

| 参数 | 说明 |
|------|------|
| `database.reference` | 参考基因组 FASTA |
| `database.gtf` | 基因注释 GTF |
| `database.star_index` | STAR 基因组索引目录 |
| `database.ctat_lib` | STAR-Fusion CTAT 资源库 |
| `database.arriba_*` | Arriba 参考文件（黑名单、已知融合、蛋白结构域、细胞带） |

### 工具参数

| 参数 | 说明 |
|------|------|
| `fastp.extra` | fastp 额外参数（当前：`--qualified_quality_phred 20 --length_required 36 --cut_front --cut_tail --cut_mean_quality 20 --detect_adapter_for_pe --trim_poly_g`） |
| `threads.low/medium/high` | 线程分级 |

## 输出

| 文件 | 说明 |
|------|------|
| `final/{sample}_fusion_report.tsv` | 完整融合报告（含注释和评分） |
| `final/{sample}_fusion_report.xlsx` | Excel 格式报告 |
| `final/{sample}_priority_fusions.tsv` | Tier 1/2 临床相关融合 |
| `annotation/{sample}_filtered_fusions.tsv` | 过滤后融合列表 |
| `annotation/{sample}_all_merged_fusions.tsv` | 全部合并融合 |
| `arriba/{sample}_fusions.pdf` | Arriba 融合可视化 |
| `qc/multiqc/{sample}_multiqc.html` | MultiQC 质控报告 |

### 融合过滤标准

保留满足任一条件的融合：

- 双工具（STAR-Fusion + Arriba）同时检出
- STAR-Fusion FFPM ≥ 0.1
- Arriba 可信度为 high

### 优先级评分

综合 `support_label`、`reading_frame`、`actionable`、`confidence`、`FFPM` 五个维度评分排序，Tier 1/2 为临床相关融合。
