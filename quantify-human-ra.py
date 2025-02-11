#!/usr/bin/env python3
import pysam
from collections import defaultdict
import os
import pandas as pd

delivery = "..."

output_dir = os.path.join(delivery, "output")
results_dir = os.path.join(output_dir, "results")
basic_stats_path = os.path.join(results_dir, "qc_basic_stats.tsv.gz")

basic_stats = pd.read_csv(basic_stats_path, sep='\t')
basic_stats_raw = basic_stats[basic_stats["stage"] == "raw_concat"].copy()
basic_stats_raw["group"] = basic_stats_raw["sample"].str.split("-div").str[0]
n_raw_reads_per_group = basic_stats_raw.groupby("group")["n_read_pairs"].sum()

human_alignments_sam = os.path.join(results_dir, "human_alignments.sam")

def is_primary_alignment(read):
    """Check if read is a primary alignment."""
    return not read.is_unmapped and not read.is_secondary and not read.is_supplementary


n_human_reads = defaultdict(int)
for read in pysam.AlignmentFile(human_alignments_sam):
    # Get read group tag (RG:Z:) from read
    sample = dict(read.tags).get('RG')
    group = sample.split("-div")[0]
    if is_primary_alignment(read):
        n_human_reads[group] += 1

print("Group\tRaw reads\tHuman reads\tRelative abundance")
for group in n_raw_reads_per_group.keys():

    percentage = (n_human_reads[group] / n_raw_reads_per_group[group]) * 100.
    print(f"{group}\t{n_raw_reads_per_group[group]}\t{n_human_reads[group]}\t{percentage:.2f}%")
