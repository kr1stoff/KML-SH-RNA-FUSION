import pandas as pd

sf = pd.read_csv(snakemake.input.sf_abridged, sep="\t", comment="#")
sf.columns = sf.columns.str.lstrip("#").str.strip()
sf_renamed = sf.rename(columns={
    "#FusionName": "fusion_name", "FusionName": "fusion_name",
    "JunctionReadCount": "junction_reads_sf", "SpanningFragCount": "spanning_reads_sf",
    "FFPM": "ffpm_sf", "LeftGene": "gene1_sf", "RightGene": "gene2_sf",
    "LeftBreakpoint": "breakpoint1_sf", "RightBreakpoint": "breakpoint2_sf",
    "LargeAnchorSupport": "large_anchor_sf", "annots": "annots_sf",
})
sf_renamed["tool_sf"] = True
sf_renamed["gene_pair"] = sf_renamed["fusion_name"].str.replace("--", "::")

arr = pd.read_csv(snakemake.input.arr_fusions, sep="\t")
arr_renamed = arr.rename(columns={
    "#gene1": "gene1_arr", "gene2": "gene2_arr",
    "breakpoint1": "breakpoint1_arr", "breakpoint2": "breakpoint2_arr",
    "split_reads1": "split_reads1_arr", "split_reads2": "split_reads2_arr",
    "discordant_mates": "discordant_mates_arr", "confidence": "confidence_arr",
    "filters": "filters_arr", "type": "type_arr",
    "direction1": "direction1_arr", "direction2": "direction2_arr",
    "reading_frame": "reading_frame_arr", "peptide_sequence": "peptide_arr",
    "fusion_transcript": "transcript_arr", "tags": "tags_arr",
    "site1": "site1_arr", "site2": "site2_arr",
    "gene_id1": "gene_id1_arr", "gene_id2": "gene_id2_arr",
    "transcript_id1": "tx_id1_arr", "transcript_id2": "tx_id2_arr",
})
arr_renamed["tool_arr"] = True
arr_renamed["gene_pair"] = arr_renamed["gene1_arr"] + "::" + arr_renamed["gene2_arr"]

def normalize_bp(bp):
    return str(bp).split(":")[0] + ":" + str(bp).split(":")[1]

sf_renamed["bp1_norm"] = sf_renamed["breakpoint1_sf"].apply(normalize_bp)
sf_renamed["bp2_norm"] = sf_renamed["breakpoint2_sf"].apply(normalize_bp)
arr_renamed["bp1_norm"] = arr_renamed["breakpoint1_arr"].apply(normalize_bp)
arr_renamed["bp2_norm"] = arr_renamed["breakpoint2_arr"].apply(normalize_bp)

merged_bp = pd.merge(sf_renamed, arr_renamed, on=["bp1_norm", "bp2_norm"], suffixes=("_sf", "_arr"), how="outer")

unmatched_sf = sf_renamed[~sf_renamed.set_index(["bp1_norm", "bp2_norm"]).index.isin(
    merged_bp.dropna(subset=["gene_pair_sf"]).set_index(["bp1_norm", "bp2_norm"]).index)]
unmatched_arr = arr_renamed[~arr_renamed.set_index(["bp1_norm", "bp2_norm"]).index.isin(
    merged_bp.dropna(subset=["gene_pair_arr"]).set_index(["bp1_norm", "bp2_norm"]).index)]

merged_gene = pd.merge(
    unmatched_sf.rename(columns={"gene_pair": "gene_pair_sf"}),
    unmatched_arr.rename(columns={"gene_pair": "gene_pair_arr"}),
    left_on="gene_pair_sf", right_on="gene_pair_arr", suffixes=("_sf2", "_arr2"), how="outer")

combined_rows = []
for _, row in merged_bp.iterrows():
    combined_rows.append(row)
for _, row in merged_gene.iterrows():
    combined_rows.append(row)

result = pd.DataFrame(combined_rows)

both_tools = result["tool_sf"].fillna(False) & result["tool_arr"].fillna(False)
high_ffpm = result["ffpm_sf"].fillna(0).astype(float) >= 0.1
high_conf = result["confidence_arr"].fillna("") == "high"

result["keep"] = both_tools | high_ffpm | high_conf
result["support_label"] = "single_tool_low"
result.loc[both_tools, "support_label"] = "both_tools"
result.loc[high_conf, "support_label"] = "arriba_high"
result.loc[both_tools & high_conf, "support_label"] = "both_tools_high_conf"

filtered = result[result["keep"]].copy()

result.to_csv(snakemake.output.all_merged, sep="\t", index=False)
filtered.to_csv(snakemake.output.filtered, sep="\t", index=False)
