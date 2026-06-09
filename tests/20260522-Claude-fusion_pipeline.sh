#!/usr/bin/env bash
# ============================================================
# RNA-seq 融合基因分析流程
# FASTQ → QC → STAR-Fusion + Arriba → 注释 & 整合
# ============================================================
# 依赖工具：
#   FastQC / MultiQC / Trim Galore
#   STAR (≥2.7.x)
#   STAR-Fusion (≥1.12)
#   Arriba (≥2.4)
#   samtools
#   Python3 + pandas
# ============================================================

set -euo pipefail

# ─────────────────────────────────────────────
# 0. 参数与路径配置
# ─────────────────────────────────────────────
SAMPLE="SAMPLE_NAME"            # 样本名称，按需修改
THREADS=16                       # CPU 线程数
READ1="data/${SAMPLE}_R1.fastq.gz"
READ2="data/${SAMPLE}_R2.fastq.gz"

# 参考基因组目录（需提前构建，详见 Step 2 注释）
GENOME_DIR="/ref/STAR_genome_hg38"
# STAR-Fusion CTAT 资源库（含 genome_lib_build_dir）
CTAT_LIB="/ref/GRCh38_gencode_v46_CTAT_lib_Apr012024"
# Arriba 参考文件目录
ARRIBA_DIR="/ref/arriba"
GTF="${ARRIBA_DIR}/gencode.v46.annotation.gtf"
GENOME_FA="${ARRIBA_DIR}/GRCh38.primary_assembly.genome.fa"
BLACKLIST="${ARRIBA_DIR}/blacklist_hg38_GRCh38_v2.4.0.tsv.gz"
KNOWN_FUSIONS="${ARRIBA_DIR}/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz"
PROTEIN_DOMAINS="${ARRIBA_DIR}/protein_domains_hg38_GRCh38_v2.4.0.gff3"
CYTOBANDS="${ARRIBA_DIR}/cytobands_hg38_GRCh38_v2.4.0.tsv"

OUT="results/${SAMPLE}"
mkdir -p "${OUT}"/{qc,trim,star,star_fusion,arriba,annotation,final}

echo "=============================="
echo " 样本: ${SAMPLE}"
echo " 开始时间: $(date)"
echo "=============================="

# ─────────────────────────────────────────────
# 1. 质控 — FastQC
# ─────────────────────────────────────────────
echo "[STEP 1] FastQC 原始数据质控..."

fastqc -t ${THREADS} \
    "${READ1}" "${READ2}" \
    -o "${OUT}/qc/"

# ─────────────────────────────────────────────
# 2. 接头与低质量碱基剪切 — Trim Galore
# ─────────────────────────────────────────────
echo "[STEP 2] Trim Galore 质量修剪..."

trim_galore \
    --paired \
    --cores ${THREADS} \
    --quality 20 \
    --length 36 \
    --stringency 3 \
    --fastqc \
    -o "${OUT}/trim/" \
    "${READ1}" "${READ2}"

TRIMMED_R1="${OUT}/trim/${SAMPLE}_R1_val_1.fq.gz"
TRIMMED_R2="${OUT}/trim/${SAMPLE}_R2_val_2.fq.gz"

# ─────────────────────────────────────────────
# 3. STAR 比对（双模式：Fusion + Arriba 共用）
#    - --chimSegmentMin 等参数同时满足两款工具
#    - 若尚未建好基因组索引，取消下方注释先运行
# ─────────────────────────────────────────────

# 【可选】构建 STAR 基因组索引（只需运行一次）
# STAR --runMode genomeGenerate \
#      --genomeDir ${GENOME_DIR} \
#      --genomeFastaFiles ${GENOME_FA} \
#      --sjdbGTFfile ${GTF} \
#      --runThreadN ${THREADS} \
#      --genomeSAindexNbases 14

echo "[STEP 3] STAR 比对..."

STAR \
    --runThreadN ${THREADS} \
    --genomeDir "${GENOME_DIR}" \
    --readFilesIn "${TRIMMED_R1}" "${TRIMMED_R2}" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMattributes NH HI AS NM MD \
    --outBAMsortingThreadN ${THREADS} \
    --outFileNamePrefix "${OUT}/star/${SAMPLE}_" \
    \
    # ── 融合基因检测必要参数 ──────────────────────
    --chimSegmentMin 12 \
    --chimJunctionOverhangMin 8 \
    --chimOutJunctionFormat 1 \
    --chimOutType WithinBAM SeparateSAMold \
    --chimNonchimScoreDrop 10 \
    --chimMultimapScoreRange 3 \
    --chimMultimapNmax 20 \
    --chimScoreJunctionNonGTAG -4 \
    --chimScoreSeparation 1 \
    \
    # ── Arriba 推荐参数 ───────────────────────────
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --alignIntronMax 500000 \
    --alignMatesGapMax 1000000 \
    --alignSJDBoverhangMin 3 \
    --outSJtype Standard \
    --twopassMode Basic \
    --peOverlapNbasesMin 12 \
    --peOverlapMMp 0.1 \
    \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 10 \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0 \
    --outFilterMatchNmin 0 \
    \
    --quantMode GeneCounts

samtools index -@ ${THREADS} "${OUT}/star/${SAMPLE}_Aligned.sortedByCoord.out.bam"

BAM="${OUT}/star/${SAMPLE}_Aligned.sortedByCoord.out.bam"
CHIM_JUNC="${OUT}/star/${SAMPLE}_Chimeric.out.junction"

# ─────────────────────────────────────────────
# 4. STAR-Fusion 融合基因检测
# ─────────────────────────────────────────────
echo "[STEP 4] STAR-Fusion 检测融合基因..."

STAR-Fusion \
    --genome_lib_dir "${CTAT_LIB}" \
    -J "${CHIM_JUNC}" \
    --output_dir "${OUT}/star_fusion/" \
    --CPU ${THREADS} \
    --FusionInspector validate \
    --examine_coding_effect \
    --denovo_reconstruct

# 主要输出：
#   star-fusion.fusion_predictions.tsv            → 所有预测
#   star-fusion.fusion_predictions.abridged.tsv   → 精简版
#   FusionInspector-validate/                     → 验证结果

# ─────────────────────────────────────────────
# 5. Arriba 融合基因检测
# ─────────────────────────────────────────────
echo "[STEP 5] Arriba 检测融合基因..."

arriba \
    -x "${BAM}" \
    -a "${GENOME_FA}" \
    -g "${GTF}" \
    -o "${OUT}/arriba/${SAMPLE}_fusions.tsv" \
    -O "${OUT}/arriba/${SAMPLE}_fusions_discarded.tsv" \
    -b "${BLACKLIST}" \
    -k "${KNOWN_FUSIONS}" \
    -t "${KNOWN_FUSIONS}" \
    -p "${PROTEIN_DOMAINS}" \
    -T \
    -P

# 可视化融合结果（生成 PDF）
draw_fusions \
    -g "${GTF}" \
    -a "${GENOME_FA}" \
    -f "${OUT}/arriba/${SAMPLE}_fusions.tsv" \
    -o "${OUT}/arriba/${SAMPLE}_fusions.pdf" \
    -c "${CYTOBANDS}" \
    -p "${PROTEIN_DOMAINS}"

# ─────────────────────────────────────────────
# 6. 整合 STAR-Fusion 与 Arriba 结果
#    并取并集（两者均检出或任一检出）
# ─────────────────────────────────────────────
echo "[STEP 6] 整合两种工具结果..."

python3 - <<'PYEOF'
import pandas as pd, os, re, sys

sample  = os.environ.get("SAMPLE", "SAMPLE_NAME")
out_dir = f"results/{sample}"

# ── 读取 STAR-Fusion 结果 ─────────────────────────────────────
sf_file = f"{out_dir}/star_fusion/star-fusion.fusion_predictions.abridged.tsv"
sf = pd.read_csv(sf_file, sep="\t", comment="#")
sf.columns = sf.columns.str.lstrip("#").str.strip()

# 标准化列名
sf_renamed = sf.rename(columns={
    "#FusionName": "fusion_name",
    "FusionName": "fusion_name",
    "JunctionReadCount": "junction_reads_sf",
    "SpanningFragCount": "spanning_reads_sf",
    "FFPM": "ffpm_sf",
    "LeftGene": "gene1_sf",
    "RightGene": "gene2_sf",
    "LeftBreakpoint": "breakpoint1_sf",
    "RightBreakpoint": "breakpoint2_sf",
    "LargeAnchorSupport": "large_anchor_sf",
    "annots": "annots_sf",
})
sf_renamed["tool_sf"] = True
sf_renamed["gene_pair"] = sf_renamed["fusion_name"].str.replace("--", "::")

# ── 读取 Arriba 结果 ──────────────────────────────────────────
arr_file = f"{out_dir}/arriba/{sample}_fusions.tsv"
arr = pd.read_csv(arr_file, sep="\t")

arr_renamed = arr.rename(columns={
    "#gene1": "gene1_arr",
    "gene2": "gene2_arr",
    "breakpoint1": "breakpoint1_arr",
    "breakpoint2": "breakpoint2_arr",
    "split_reads1": "split_reads1_arr",
    "split_reads2": "split_reads2_arr",
    "discordant_mates": "discordant_mates_arr",
    "confidence": "confidence_arr",
    "filters": "filters_arr",
    "type": "type_arr",
    "direction1": "direction1_arr",
    "direction2": "direction2_arr",
    "reading_frame": "reading_frame_arr",
    "peptide_sequence": "peptide_arr",
    "fusion_transcript": "transcript_arr",
    "tags": "tags_arr",
    "site1": "site1_arr",
    "site2": "site2_arr",
    "gene_id1": "gene_id1_arr",
    "gene_id2": "gene_id2_arr",
    "transcript_id1": "tx_id1_arr",
    "transcript_id2": "tx_id2_arr",
})
arr_renamed["tool_arr"] = True
arr_renamed["gene_pair"] = arr_renamed["gene1_arr"] + "::" + arr_renamed["gene2_arr"]

# ── 规范化断点（去除链信息，便于比对）────────────────────────
def normalize_bp(bp):
    return str(bp).split(":")[0] + ":" + str(bp).split(":")[1] if ":" in str(bp) else str(bp)

sf_renamed["bp1_norm"] = sf_renamed["breakpoint1_sf"].apply(normalize_bp)
sf_renamed["bp2_norm"] = sf_renamed["breakpoint2_sf"].apply(normalize_bp)
arr_renamed["bp1_norm"] = arr_renamed["breakpoint1_arr"].apply(normalize_bp)
arr_renamed["bp2_norm"] = arr_renamed["breakpoint2_arr"].apply(normalize_bp)

# ── 合并：先尝试断点级精确匹配，再回落到基因对匹配 ──────────
merged_bp = pd.merge(
    sf_renamed, arr_renamed,
    on=["bp1_norm", "bp2_norm"],
    suffixes=("_sf", "_arr"),
    how="outer",
)

# 对未匹配的行，用基因对合并
unmatched_sf  = sf_renamed[~sf_renamed.set_index(["bp1_norm","bp2_norm"]).index.isin(
    merged_bp.dropna(subset=["gene_pair_sf"]).set_index(["bp1_norm","bp2_norm"]).index)]
unmatched_arr = arr_renamed[~arr_renamed.set_index(["bp1_norm","bp2_norm"]).index.isin(
    merged_bp.dropna(subset=["gene_pair_arr"]).set_index(["bp1_norm","bp2_norm"]).index)]

merged_gene = pd.merge(
    unmatched_sf.rename(columns={"gene_pair": "gene_pair_sf"}),
    unmatched_arr.rename(columns={"gene_pair": "gene_pair_arr"}),
    left_on="gene_pair_sf", right_on="gene_pair_arr",
    suffixes=("_sf2","_arr2"),
    how="outer",
)

# 汇总
final_cols = [
    "fusion_name", "gene_pair_sf", "gene_pair_arr",
    "breakpoint1_sf", "breakpoint2_sf",
    "breakpoint1_arr", "breakpoint2_arr",
    "junction_reads_sf", "spanning_reads_sf", "ffpm_sf",
    "split_reads1_arr", "split_reads2_arr", "discordant_mates_arr",
    "confidence_arr", "reading_frame_arr", "type_arr",
    "annots_sf", "tags_arr",
    "tool_sf", "tool_arr",
]
combined_rows = []
for _, row in merged_bp.iterrows():
    combined_rows.append(row)

result = pd.DataFrame(combined_rows)

# ── 质量过滤 ─────────────────────────────────────────────────
# 双工具支持，或 FFPM ≥ 0.1，或 Arriba 可信度为 high
both_tools = result["tool_sf"].fillna(False) & result["tool_arr"].fillna(False)
high_ffpm   = result["ffpm_sf"].fillna(0).astype(float) >= 0.1
high_conf   = result["confidence_arr"].fillna("") == "high"

result["keep"] = both_tools | high_ffpm | high_conf
result["support_label"] = "single_tool_low"
result.loc[both_tools, "support_label"] = "both_tools"
result.loc[high_conf, "support_label"] = "arriba_high"
result.loc[both_tools & high_conf, "support_label"] = "both_tools_high_conf"

filtered = result[result["keep"]].copy()
print(f"  原始融合: STAR-Fusion={len(sf_renamed)}, Arriba={len(arr_renamed)}")
print(f"  整合后合并: {len(result)} 条")
print(f"  过滤后保留: {len(filtered)} 条")

# ── 保存 ────────────────────────────────────────────────────
result.to_csv(f"{out_dir}/annotation/all_merged_fusions.tsv", sep="\t", index=False)
filtered.to_csv(f"{out_dir}/annotation/filtered_fusions.tsv", sep="\t", index=False)
print(f"  已写出: {out_dir}/annotation/filtered_fusions.tsv")
PYEOF

# ─────────────────────────────────────────────
# 7. 融合基因功能注释
#    7a. 对 filtered_fusions.tsv 做已知融合注释
#    7b. 蛋白阅读框 & 结构域注释（Arriba 已内置）
#    7c. 外部数据库比对（COSMIC / ChimerDB）
# ─────────────────────────────────────────────
echo "[STEP 7] 功能注释..."

python3 - <<'PYEOF'
import pandas as pd, os, re

sample  = os.environ.get("SAMPLE", "SAMPLE_NAME")
out_dir = f"results/{sample}"

df = pd.read_csv(f"{out_dir}/annotation/filtered_fusions.tsv", sep="\t")

# ── 7a. 临床相关性标注（示例基因列表，可替换为实际数据库） ────
CLINICAL_FUSIONS = {
    # 格式: frozenset({gene1, gene2}): {"cancer_type": ..., "significance": ...}
    frozenset({"BCR","ABL1"}):    {"cancer": "CML/ALL",  "sig": "Oncogenic driver"},
    frozenset({"EML4","ALK"}):    {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"TMPRSS2","ERG"}): {"cancer": "Prostate", "sig": "Oncogenic driver"},
    frozenset({"SS18","SSX1"}):   {"cancer": "Synovial sarcoma", "sig": "Diagnostic marker"},
    frozenset({"FUS","DDIT3"}):   {"cancer": "Liposarcoma", "sig": "Diagnostic marker"},
    frozenset({"EWSR1","FLI1"}):  {"cancer": "Ewing sarcoma", "sig": "Diagnostic marker"},
    frozenset({"KMT2A","AFF1"}):  {"cancer": "ALL",      "sig": "Oncogenic driver"},
    frozenset({"NPM1","ALK"}):    {"cancer": "ALCL",     "sig": "Oncogenic driver"},
    frozenset({"PAX3","FOXO1"}):  {"cancer": "Rhabdomyosarcoma","sig": "Oncogenic driver"},
    frozenset({"RET","KIF5B"}):   {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"ROS1","CD74"}):   {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"NTRK1","TPM3"}):  {"cancer": "Various",  "sig": "Oncogenic driver"},
    frozenset({"MYB","NFIB"}):    {"cancer": "ACC",      "sig": "Diagnostic marker"},
}

TIER_LABELS = {
    "Oncogenic driver": "Tier 1 — Actionable",
    "Diagnostic marker": "Tier 2 — Diagnostic",
}

def annotate_row(row):
    # 从两个工具各自提取基因名
    genes = set()
    for col in ["gene_pair_sf", "gene_pair_arr", "fusion_name"]:
        val = str(row.get(col, ""))
        for sep in ["::", "--", "&"]:
            if sep in val:
                g1, g2 = val.split(sep, 1)
                genes.update([g1.split("^")[0].strip(), g2.split("^")[0].strip()])
    if len(genes) < 2:
        return pd.Series({"clinical_cancer": "Unknown", "clinical_sig": "VUS",
                          "tier": "Tier 4 — Unknown", "actionable": False})
    key = frozenset(list(genes)[:2])
    for k, v in CLINICAL_FUSIONS.items():
        if k.issubset(genes) or genes.issuperset(k):
            return pd.Series({
                "clinical_cancer": v["cancer"],
                "clinical_sig":    v["sig"],
                "tier":            TIER_LABELS.get(v["sig"], "Tier 3 — Uncertain"),
                "actionable":      v["sig"] == "Oncogenic driver",
            })
    return pd.Series({"clinical_cancer": "N/A", "clinical_sig": "VUS",
                      "tier": "Tier 3 — Uncertain", "actionable": False})

annot = df.apply(annotate_row, axis=1)
df = pd.concat([df, annot], axis=1)

# ── 7b. 增加阅读框评分 ───────────────────────────────────────
def frame_score(row):
    rf = str(row.get("reading_frame_arr", "")).lower()
    if "in-frame" in rf:    return 2
    if "out-of-frame" in rf: return 0
    return 1  # 未知

df["frame_score"] = df.apply(frame_score, axis=1)

# ── 7c. 综合评分 ─────────────────────────────────────────────
def priority_score(row):
    score = 0
    score += 3 if str(row.get("support_label","")).startswith("both") else 1
    score += row.get("frame_score", 0) * 2
    score += 2 if row.get("actionable", False) else 0
    score += 1 if str(row.get("confidence_arr","")).lower() == "high" else 0
    ffpm = float(row.get("ffpm_sf", 0) or 0)
    score += 2 if ffpm >= 1.0 else (1 if ffpm >= 0.1 else 0)
    return score

df["priority_score"] = df.apply(priority_score, axis=1)
df = df.sort_values("priority_score", ascending=False)

# ── 保存最终报告 ─────────────────────────────────────────────
REPORT_COLS = [
    "fusion_name", "gene_pair_sf", "gene_pair_arr",
    "breakpoint1_sf", "breakpoint2_sf",
    "breakpoint1_arr", "breakpoint2_arr",
    "junction_reads_sf", "spanning_reads_sf", "ffpm_sf",
    "split_reads1_arr", "split_reads2_arr", "discordant_mates_arr",
    "confidence_arr", "reading_frame_arr", "type_arr",
    "support_label", "annots_sf", "tags_arr",
    "clinical_cancer", "clinical_sig", "tier", "actionable",
    "frame_score", "priority_score",
]
existing_cols = [c for c in REPORT_COLS if c in df.columns]
report = df[existing_cols]

report.to_csv(f"{out_dir}/final/{os.environ.get('SAMPLE','SAMPLE')}_fusion_report.tsv",
              sep="\t", index=False)
report.to_excel(f"{out_dir}/final/{os.environ.get('SAMPLE','SAMPLE')}_fusion_report.xlsx",
                index=False)

# 仅输出 Tier 1/2
priority = report[report["tier"].str.startswith(("Tier 1","Tier 2"))]
priority.to_csv(f"{out_dir}/final/{os.environ.get('SAMPLE','SAMPLE')}_priority_fusions.tsv",
                sep="\t", index=False)

print(f"  总注释融合: {len(report)}")
print(f"  Tier 1/2 (临床相关): {len(priority)}")
print(f"  最终报告: {out_dir}/final/")
PYEOF

# ─────────────────────────────────────────────
# 8. MultiQC 汇总报告
# ─────────────────────────────────────────────
echo "[STEP 8] MultiQC 汇总..."

multiqc \
    "${OUT}/qc/" \
    "${OUT}/trim/" \
    "${OUT}/star/" \
    -n "${SAMPLE}_multiqc" \
    -o "${OUT}/qc/"

# ─────────────────────────────────────────────
# 完成
# ─────────────────────────────────────────────
echo "=============================="
echo " 流程完成: $(date)"
echo " 主要输出文件:"
echo "   QC 报告   : ${OUT}/qc/${SAMPLE}_multiqc.html"
echo "   STAR-Fusion: ${OUT}/star_fusion/star-fusion.fusion_predictions.abridged.tsv"
echo "   Arriba    : ${OUT}/arriba/${SAMPLE}_fusions.tsv"
echo "   Arriba 图 : ${OUT}/arriba/${SAMPLE}_fusions.pdf"
echo "   整合结果  : ${OUT}/annotation/filtered_fusions.tsv"
echo "   最终报告  : ${OUT}/final/${SAMPLE}_fusion_report.tsv"
echo "   优先融合  : ${OUT}/final/${SAMPLE}_priority_fusions.tsv"
echo "=============================="
