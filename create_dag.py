#!/usr/bin/env python3
"""
Generate a publication-quality DAG figure for the Freebayes-Loop pipeline.
Requires: graphviz (Python package) + dot (system binary)
"""

import graphviz

dot = graphviz.Digraph(
    "Freebayes_Loop_DAG",
    format="png",
    graph_attr={
        "rankdir": "TB",
        "fontname": "Helvetica",
        "fontsize": "14",
        "label": "Freebayes-Loop v1.0 — Pipeline DAG\nRegion-based Parallel Variant Calling",
        "labelloc": "t",
        "labeljust": "c",
        "bgcolor": "white",
        "dpi": "200",
        "pad": "0.5",
        "nodesep": "0.6",
        "ranksep": "0.8",
    },
)

# --- Styles ---
input_style = {
    "shape": "folder",
    "style": "filled",
    "fillcolor": "#E8F5E9",
    "color": "#388E3C",
    "fontname": "Helvetica",
    "fontsize": "10",
}
process_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#E3F2FD",
    "color": "#1565C0",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
merge_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#FFF3E0",
    "color": "#E65100",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
helper_style = {
    "shape": "box",
    "style": "filled,rounded",
    "fillcolor": "#EDE7F6",
    "color": "#4527A0",
    "fontname": "Helvetica-Bold",
    "fontsize": "10",
    "penwidth": "2",
}
output_style = {
    "shape": "folder",
    "style": "filled",
    "fillcolor": "#FFEBEE",
    "color": "#C62828",
    "fontname": "Helvetica",
    "fontsize": "10",
}

edge_main = {"color": "#1565C0", "penwidth": "1.8"}
edge_merge = {"color": "#E65100", "penwidth": "1.8"}
edge_helper = {"color": "#4527A0", "penwidth": "1.8", "style": "dashed"}
edge_output = {"color": "#00695C", "penwidth": "1.8"}

# --- Input nodes ---
dot.node("input_bed", "Input BED regions\n(start_range.bed)", **input_style)
dot.node("bam_list", "BAM file list\n(bam_paths.txt)", **input_style)
dot.node("ref_fasta", "Reference FASTA\n(hpv16_ref.fa)", **input_style)

# --- Channel split ---
dot.node(
    "split",
    "splitText / splitCsv\n─────────────────────\nParse BED into\n(chr, start, end) tuples",
    shape="box",
    style="filled,rounded",
    fillcolor="#E0F2F1",
    color="#00695C",
    fontname="Helvetica-Bold",
    fontsize="10",
    penwidth="2",
)

# --- freebayes process ---
dot.node(
    "freebayes",
    "freebayes\n─────────────────────\nVariant calling per region\n(pooled, ploidy=4)\n→ variant_*.vcf.gz",
    **process_style,
)

# --- collect ---
dot.node(
    "collect",
    "collect()\n─────────────────────\nGather all region VCFs",
    shape="ellipse",
    style="filled",
    fillcolor="#F3E5F5",
    color="#7B1FA2",
    fontname="Helvetica",
    fontsize="10",
)

# --- merge_vcfs process ---
dot.node(
    "merge_vcfs",
    "merge_vcfs\n─────────────────────\nbcftools concat\n+ vcf2bed_03.py\n(deletion detection)",
    **merge_style,
)

# --- Helper scripts ---
dot.node(
    "vcf2bed",
    "vcf2bed_03.py\n─────────────────────\nExtract large deletions\nfrom merged VCF → BED",
    **helper_style,
)
dot.node(
    "make_bed",
    "make_bed.sh\n─────────────────────\nGenerate BED regions\nfrom genome size",
    **helper_style,
)
dot.node(
    "update_lib",
    "update_deletion_lib.py\n─────────────────────\nUpdate deletion library\nwith new findings",
    **helper_style,
)

# --- Outputs ---
dot.node("vcf_out", "merge.vcf.gz\n(all variants)", **output_style)
dot.node("bed_out", "merge.vcf.gz.bed\n(large deletions BED)", **output_style)
dot.node("log_out", "vcf.log\n(merge log)", **output_style)

# --- Edges ---
# Input → split
dot.edge("input_bed", "split", **edge_main)

# Helper for BED generation
dot.edge("make_bed", "input_bed", label="  generates  ", **edge_helper)

# split → freebayes
dot.edge("split", "freebayes", label="  (chr, start, end)  ", **edge_main)
dot.edge("bam_list", "freebayes", style="dotted", color="#388E3C")
dot.edge("ref_fasta", "freebayes", style="dotted", color="#388E3C")

# freebayes → collect
dot.edge("freebayes", "collect", label="  .vcf.gz + .csi  ", **edge_main)

# collect → merge_vcfs
dot.edge("collect", "merge_vcfs", label="  all VCFs  ", **edge_merge)

# vcf2bed helper
dot.edge("vcf2bed", "merge_vcfs", label="  called by  ", **edge_helper)

# merge_vcfs → outputs
dot.edge("merge_vcfs", "vcf_out", **edge_output)
dot.edge("merge_vcfs", "bed_out", **edge_output)
dot.edge("merge_vcfs", "log_out", **edge_output)

# Deletion library update (post-pipeline)
dot.edge("bed_out", "update_lib", label="  new deletions  ", **edge_helper)

# --- Legend ---
with dot.subgraph(name="cluster_legend") as legend:
    legend.attr(
        label="Legend",
        style="dashed",
        color="gray",
        fontname="Helvetica-Bold",
        fontsize="11",
    )
    legend.node("leg1", "Input / Output", shape="folder", style="filled",
                fillcolor="#E8F5E9", color="#388E3C", fontsize="9", fontname="Helvetica")
    legend.node("leg2", "Variant calling", shape="box", style="filled,rounded",
                fillcolor="#E3F2FD", color="#1565C0", fontsize="9", fontname="Helvetica")
    legend.node("leg3", "VCF merging", shape="box", style="filled,rounded",
                fillcolor="#FFF3E0", color="#E65100", fontsize="9", fontname="Helvetica")
    legend.node("leg4", "Helper scripts", shape="box", style="filled,rounded",
                fillcolor="#EDE7F6", color="#4527A0", fontsize="9", fontname="Helvetica")
    legend.edge("leg1", "leg2", style="invis")
    legend.edge("leg2", "leg3", style="invis")
    legend.edge("leg3", "leg4", style="invis")

output_path = "/tmp/freebayes_loop/dag"
dot.render(output_path, cleanup=True)
print(f"DAG saved to {output_path}.png")
