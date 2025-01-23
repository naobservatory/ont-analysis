#! /usr/bin/env python3

import argparse
from collections import defaultdict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def identify_complex_reads(fasta_file):
    """
    Identify reads that are complex (not masked by dustmasker).
    Returns a dictionary of read IDs and their sequences that have no masked regions.
    Dustmasker marks low complexity regions with lowercase letters.
    """
    good_seqs = defaultdict(str)
    with open(fasta_file) as inf:
        for record in SeqIO.parse(inf, "fasta"):
            seq_str = str(record.seq)
            # Keep sequences that have no lowercase letters (no masked regions)
            if seq_str.isupper():
                good_seqs[record.id] = record.seq
    return good_seqs


def main():
    parser = argparse.ArgumentParser(
        description="Filter out low complexity reads from a FASTA file"
    )
    parser.add_argument("input_fasta", help="Input FASTA file with dustmasker results")
    parser.add_argument("output_fasta", help="Output FASTA file for complex reads")
    args = parser.parse_args()

    # Get complex reads
    good_seqs = identify_complex_reads(args.input_fasta)

    # Write filtered sequences to output file
    with open(args.output_fasta, "w") as fasta_handle:
        for read_id, seq in good_seqs.items():
            SeqIO.write(
                SeqRecord(seq, id=read_id, description=""), fasta_handle, "fasta"
            )

    print(f"Wrote {len(good_seqs)} complex reads to {args.output_fasta}")


if __name__ == "__main__":
    main()
