#!/usr/bin/env python3

import os
import subprocess
import glob
import argparse
from collections import Counter, defaultdict

def get_long_reads(fastq_file, min_length=1000):
    """Extract read IDs of sequences longer than min_length using seqtk."""
    temp_output = "temp_seqtk_output.fastq"
    subprocess.run(f"seqtk seq -L {min_length} {fastq_file} > {temp_output}", shell=True)

    read_ids = {line.strip().split()[0][1:] for i, line in enumerate(open(temp_output))
                if i % 4 == 0}  # Get header lines and extract IDs
    os.remove(temp_output)
    return read_ids

def process_assignment_file(filename, read_ids):
    """Process assignment file and return read assignments and counts."""
    read_assignments = defaultdict(list)
    taxon_counter = Counter()

    with open(filename) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[1] in read_ids:
                read_assignments[parts[1]].append(parts[2])
                taxon_counter[parts[2]] += 1

    return read_assignments, taxon_counter

def main():
    parser = argparse.ArgumentParser(description="Process long reads and analyze Kraken2 assignments")
    parser.add_argument("-c", "--cleaned-dir", required=True, help="Directory containing cleaned fastq files")
    parser.add_argument("-p", "--fastq-prefix", required=True, help="Prefix of fastq files to process")
    parser.add_argument("-k", "--kraken-dir", required=True, help="Directory containing Kraken output files")
    args = parser.parse_args()

    fastq_pattern = os.path.join(args.cleaned_dir, f"{args.fastq_prefix}*.fastq.gz")

    # Get long reads from all matching files
    all_long_reads = set()
    for fastq_file in glob.glob(fastq_pattern):
        print(f"Processing {fastq_file}...")
        all_long_reads.update(get_long_reads(fastq_file))

    print(f"\nFound {len(all_long_reads)} long reads")

    # Process assignments
    ribo_read_assignments = defaultdict(list)
    non_ribo_read_assignments = defaultdict(list)
    ribo_taxon_counter = Counter()
    non_ribo_taxon_counter = Counter()

    # Process non-ribo assignments
    non_ribo_file = os.path.join(args.kraken_dir, "PZ-250116-BoDT-NAS-R1-no-ribo.output")
    if os.path.exists(non_ribo_file):
        print(f"Processing non-ribo assignments from {non_ribo_file}...")
        non_ribo_read_assignments, non_ribo_taxon_counter = process_assignment_file(non_ribo_file, all_long_reads)
    else:
        print(f"Warning: {non_ribo_file} not found!")

    # Process ribo assignments
    ribo_file = os.path.join(args.kraken_dir, "PZ-250116-BoDT-NAS-R1-ribo.output")
    if os.path.exists(ribo_file):
        print(f"Processing ribo assignments from {ribo_file}...")
        ribo_read_assignments, ribo_taxon_counter = process_assignment_file(ribo_file, all_long_reads)
    else:
        print(f"Warning: {ribo_file} not found!")

    # Print results
    print(f"\nSummary:")
    print(f"Total long reads found: {len(all_long_reads)}")

    print(f"\nNon-ribo results:")
    print(f"Reads with assignments: {len(non_ribo_read_assignments)}")
    print(f"Unique taxa found: {len(non_ribo_taxon_counter)}")
    print("\nTop 10 non-ribo taxa by assignment count:")
    for taxon, count in non_ribo_taxon_counter.most_common(10):
        print(f"{taxon}: {count} occurrences")

    print(f"\nRibo results:")
    print(f"Reads with assignments: {len(ribo_read_assignments)}")
    print(f"Unique taxa found: {len(ribo_taxon_counter)}")
    print("\nTop 10 ribo taxa by assignment count:")
    for taxon, count in ribo_taxon_counter.most_common(10):
        print(f"{taxon}: {count} occurrences")

if __name__ == "__main__":
    main()
