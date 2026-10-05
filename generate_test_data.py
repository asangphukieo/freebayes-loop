#!/usr/bin/env python3
"""
Generate mock test data for freebayes_loop pipeline testing.
Creates:
  - A small HPV16-like reference FASTA (500 bp)
  - Two mock BAM files with paired-end reads aligned to the reference
  - A BED file splitting the reference into regions
  - A bam_paths.txt listing the BAM files
  - An example start_range.bed for the pipeline

Requires: pysam
"""

import pysam
import os
import random

random.seed(42)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_data")
os.makedirs(OUT_DIR, exist_ok=True)

REF_NAME = "gi|333031|lcl|HPV16REF.1|"
REF_LEN = 500

# --- Generate reference FASTA ---
def random_seq(length):
    return "".join(random.choice("ACGT") for _ in range(length))

ref_seq = random_seq(REF_LEN)
ref_fasta = os.path.join(OUT_DIR, "hpv16_ref.fa")
with open(ref_fasta, "w") as f:
    f.write(f">{REF_NAME}\n")
    for i in range(0, len(ref_seq), 80):
        f.write(ref_seq[i:i+80] + "\n")

# Index with pysam
pysam.faidx(ref_fasta)

# --- Generate BAM files ---
def make_bam(bam_path, ref_name, ref_len, ref_seq, n_pairs=100):
    """Create a mock BAM with paired-end reads."""
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": ref_name, "LN": ref_len}],
    })

    read_len = 75
    insert_size = 200

    reads = []
    for i in range(n_pairs):
        pos1 = random.randint(0, ref_len - insert_size - 1)
        pos2 = pos1 + insert_size - read_len

        # Some reads have a small deletion (to test vcf2bed)
        has_deletion = (i % 20 == 0) and (pos1 + read_len + 15 < ref_len)

        # Read 1 (forward)
        a = pysam.AlignedSegment(header)
        a.query_name = f"read_{i}"
        a.query_sequence = ref_seq[pos1:pos1+read_len]
        a.flag = 99  # paired, proper pair, mate reverse, first in pair
        a.reference_id = 0
        a.reference_start = pos1
        a.mapping_quality = 60
        if has_deletion:
            # 60M 15D 15M
            a.cigartuples = [(0, 60), (2, 15), (0, 15)]
        else:
            a.cigartuples = [(0, read_len)]
        a.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        a.next_reference_id = 0
        a.next_reference_start = pos2
        a.template_length = insert_size
        reads.append(a)

        # Read 2 (reverse)
        b = pysam.AlignedSegment(header)
        b.query_name = f"read_{i}"
        b.query_sequence = ref_seq[pos2:pos2+read_len]
        b.flag = 147  # paired, proper pair, reverse, second in pair
        b.reference_id = 0
        b.reference_start = pos2
        b.mapping_quality = 60
        b.cigartuples = [(0, read_len)]
        b.query_qualities = pysam.qualitystring_to_array("I" * read_len)
        b.next_reference_id = 0
        b.next_reference_start = pos1
        b.template_length = -insert_size
        reads.append(b)

    # Sort by position
    reads.sort(key=lambda x: (x.reference_id, x.reference_start))

    with pysam.AlignmentFile(bam_path, "wb", header=header) as outf:
        for r in reads:
            outf.write(r)

    pysam.sort("-o", bam_path, bam_path)
    pysam.index(bam_path)


# Create two sample BAMs
for sample in ["sample1", "sample2"]:
    bam_file = os.path.join(OUT_DIR, f"{sample}.bam")
    make_bam(bam_file, REF_NAME, REF_LEN, ref_seq, n_pairs=100)
    print(f"Created {bam_file}")

# --- Create bam_paths.txt ---
bam_paths_file = os.path.join(OUT_DIR, "bam_paths.txt")
with open(bam_paths_file, "w") as f:
    for sample in ["sample1", "sample2"]:
        f.write(os.path.join(OUT_DIR, f"{sample}.bam") + "\n")
print(f"Created {bam_paths_file}")

# --- Create start_range.bed ---
# Split the 500 bp reference into 2 regions
bed_file = os.path.join(OUT_DIR, "start_range.bed")
with open(bed_file, "w") as f:
    f.write(f"{REF_NAME} 1 250\n")
    f.write(f"{REF_NAME} 251 {REF_LEN}\n")
print(f"Created {bed_file}")

# --- Create example output ---
example_dir = os.path.join(OUT_DIR, "example_output")
os.makedirs(example_dir, exist_ok=True)

# Example deletion BED
with open(os.path.join(example_dir, "merge.vcf.gz.bed"), "w") as f:
    f.write(f"{REF_NAME} 120 135\n")
    f.write(f"{REF_NAME} 300 315\n")
print(f"Created example output in {example_dir}/")

print("\nTest data generation complete!")
