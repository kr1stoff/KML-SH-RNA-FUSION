import pandas as pd

df = pd.read_csv(snakemake.input.filtered, sep="\t")

CLINICAL_FUSIONS = {
    frozenset({"BCR", "ABL1"}):    {"cancer": "CML/ALL",  "sig": "Oncogenic driver"},
    frozenset({"EML4", "ALK"}):    {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"TMPRSS2", "ERG"}): {"cancer": "Prostate", "sig": "Oncogenic driver"},
    frozenset({"SS18", "SSX1"}):   {"cancer": "Synovial sarcoma", "sig": "Diagnostic marker"},
    frozenset({"FUS", "DDIT3"}):   {"cancer": "Liposarcoma", "sig": "Diagnostic marker"},
    frozenset({"EWSR1", "FLI1"}):  {"cancer": "Ewing sarcoma", "sig": "Diagnostic marker"},
    frozenset({"KMT2A", "AFF1"}):  {"cancer": "ALL",      "sig": "Oncogenic driver"},
    frozenset({"NPM1", "ALK"}):    {"cancer": "ALCL",     "sig": "Oncogenic driver"},
    frozenset({"PAX3", "FOXO1"}):  {"cancer": "Rhabdomyosarcoma", "sig": "Oncogenic driver"},
    frozenset({"RET", "KIF5B"}):   {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"ROS1", "CD74"}):   {"cancer": "NSCLC",    "sig": "Oncogenic driver"},
    frozenset({"NTRK1", "TPM3"}):  {"cancer": "Various",  "sig": "Oncogenic driver"},
    frozenset({"MYB", "NFIB"}):    {"cancer": "ACC",      "sig": "Diagnostic marker"},
}

TIER_LABELS = {
    "Oncogenic driver": "Tier 1 \u2014 Actionable",
    "Diagnostic marker": "Tier 2 \u2014 Diagnostic",
}

def annotate_row(row):
    genes = set()
    for col in ["gene_pair_sf", "gene_pair_arr", "fusion_name"]:
        val = str(row.get(col, ""))
        for sep in ["::", "--", "&"]:
            if sep in val:
                g1, g2 = val.split(sep, 1)
                genes.update([g1.split("^")[0].strip(), g2.split("^")[0].strip()])
    if len(genes) < 2:
        return pd.Series({"clinical_cancer": "Unknown", "clinical_sig": "VUS",
                          "tier": "Tier 4 \u2014 Unknown", "actionable": False})
    for k, v in CLINICAL_FUSIONS.items():
        if k.issubset(genes) or genes.issuperset(k):
            return pd.Series({
                "clinical_cancer": v["cancer"],
                "clinical_sig": v["sig"],
                "tier": TIER_LABELS.get(v["sig"], "Tier 3 \u2014 Uncertain"),
                "actionable": v["sig"] == "Oncogenic driver",
            })
    return pd.Series({"clinical_cancer": "N/A", "clinical_sig": "VUS",
                      "tier": "Tier 3 \u2014 Uncertain", "actionable": False})

annot = df.apply(annotate_row, axis=1)
df = pd.concat([df, annot], axis=1)

def frame_score(row):
    rf = str(row.get("reading_frame_arr", "")).lower()
    if "in-frame" in rf:    return 2
    if "out-of-frame" in rf: return 0
    return 1

df["frame_score"] = df.apply(frame_score, axis=1)

def priority_score(row):
    score = 0
    score += 3 if str(row.get("support_label", "")).startswith("both") else 1
    score += row.get("frame_score", 0) * 2
    score += 2 if row.get("actionable", False) else 0
    score += 1 if str(row.get("confidence_arr", "")).lower() == "high" else 0
    ffpm = float(row.get("ffpm_sf", 0) or 0)
    score += 2 if ffpm >= 1.0 else (1 if ffpm >= 0.1 else 0)
    return score

df["priority_score"] = df.apply(priority_score, axis=1)
df = df.sort_values("priority_score", ascending=False)

REPORT_COLS = [
    "fusion_name", "gene_pair_sf", "gene_pair_arr",
    "breakpoint1_sf", "breakpoint2_sf", "breakpoint1_arr", "breakpoint2_arr",
    "junction_reads_sf", "spanning_reads_sf", "ffpm_sf",
    "split_reads1_arr", "split_reads2_arr", "discordant_mates_arr",
    "confidence_arr", "reading_frame_arr", "type_arr",
    "support_label", "annots_sf", "tags_arr",
    "clinical_cancer", "clinical_sig", "tier", "actionable",
    "frame_score", "priority_score",
]
existing_cols = [c for c in REPORT_COLS if c in df.columns]
report = df[existing_cols]

report.to_csv(snakemake.output.report_tsv, sep="\t", index=False)
report.to_excel(snakemake.output.report_xlsx, index=False)

priority = report[report["tier"].str.startswith(("Tier 1", "Tier 2"))]
priority.to_csv(snakemake.output.priority, sep="\t", index=False)
