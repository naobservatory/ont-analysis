#!/usr/bin/env python3

import argparse
import gzip
import json
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from collections import defaultdict


def process_fastq_file(filepath):
    """Process a single fastq file and return length distribution."""
    lengths = defaultdict(int)

    with gzip.open(filepath, "rt") as inf:
        for title, sequence, quality in FastqGeneralIterator(inf):
            seql = len(sequence)
            lengths[seql] += 1

    return dict(lengths)


def extract_barcode(filename):
    """Extract barcode base from filename."""
    # Split by dash and take everything before the div part
    parts = filename.split("-div")[0]
    return parts


def main():
    parser = argparse.ArgumentParser(
        description="Process fastq files and aggregate read lengths by barcode"
    )
    parser.add_argument("input_dir", help="Directory containing fastq.gz files")
    parser.add_argument("output_json", help="Output JSON file path")

    args = parser.parse_args()

    # Dictionary to store results for all barcodes
    barcode_lengths = {}

    # Process all fastq.gz files in the directory
    for filename in os.listdir(args.input_dir):
        if filename.endswith(".fastq.gz"):
            filepath = os.path.join(args.input_dir, filename)
            barcode = extract_barcode(filename)

            print(f"Processing {filename}...")

            # Process the file and get length distribution
            lengths = process_fastq_file(filepath)

            # Convert defaultdict to regular dict for JSON serialization
            if barcode not in barcode_lengths:
                barcode_lengths[barcode] = lengths
            else:
                # Merge length counts if we have multiple files for same barcode
                for length, count in lengths.items():
                    if length in barcode_lengths[barcode]:
                        barcode_lengths[barcode][length] += count
                    else:
                        barcode_lengths[barcode][length] = count

    # Write consolidated results to JSON file
    with open(args.output_json, "w") as outf:
        json.dump(barcode_lengths, outf, indent=2, sort_keys=True)

    print(f"Results written to {args.output_json}")


if __name__ == "__main__":
    main()
