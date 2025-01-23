#! /usr/bin/env python3
import argparse
import os
import pandas as pd
import gzip
import random
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter


def has_reverse_complement_kmer(read, k):
    seq = read
    len_seq = len(seq)
    kmers = []
    reverse_complements = []
    # Generate all possible k-mers
    for i in range(len_seq - k + 1):
        kmer = seq[i : i + k]
        kmer = Seq(kmer)
        rc_kmer = kmer.reverse_complement()
        kmers.append(str(kmer))
        reverse_complements.append(str(rc_kmer))

    # Check for reverse complement k-mers
    kmer_set = set(kmers)
    rc_set = set(reverse_complements)
    kmer_rc_overlaps = kmer_set.intersection(rc_set)

    if len(kmer_rc_overlaps) == 0:
        contains_reverse_complement_kmer = False
    else:
        contains_reverse_complement_kmer = True

    # Check for duplicate k-mers
    kmer_counter = Counter(kmers)
    if any(count > 1 for count in kmer_counter.values()):
        contains_duplicate_kmers = True
    else:
        contains_duplicate_kmers = False

    return contains_reverse_complement_kmer, contains_duplicate_kmers


def randomize_and_subset_reads(sequence_file, subset):
    # if gzipped fastq, open as gzip
    if sequence_file.endswith("fastq.gz"):
        with gzip.open(sequence_file, "rt") as inf:
            records = list(SeqIO.parse(inf, "fastq"))
    # if fasta, open as fasta
    elif sequence_file.endswith(".fasta"):
        with open(sequence_file, "r") as inf:
            records = list(SeqIO.parse(inf, "fasta"))
    else:
        raise ValueError(f"Unknown file type: {sequence_file}")
    random.shuffle(records)
    records = records[:subset]
    return records


def analyse_fastq(sequence_file, k, subset):
    # In a basecalled fastq, duplex reads are not in random order.
    # Randomize the fastq to get a random subset of reads.
    records = randomize_and_subset_reads(sequence_file, subset)

    n_reads = 0
    n_rc_reads = 0
    n_duplicate_kmers = 0
    n_duplicate_rc_reads = 0

    for record in records:
        n_reads += 1
        contains_reverse_complement_kmer, contains_duplicate_kmers = (
            has_reverse_complement_kmer(str(record.seq), k)
        )

        if contains_reverse_complement_kmer:
            n_rc_reads += 1

        if contains_duplicate_kmers:
            n_duplicate_kmers += 1

        if contains_reverse_complement_kmer and contains_duplicate_kmers:
            n_duplicate_rc_reads += 1

    if n_reads == 0:
        fraction_rc_reads = "0.00%"
        fraction_duplicate_kmers = "0.00%"
        fraction_duplicate_rc_reads = "0.00%"
    else:
        fraction_rc_reads = f"{n_rc_reads / n_reads * 100:.2f}%"
        fraction_duplicate_kmers = f"{n_duplicate_kmers / n_reads * 100:.2f}%"
        fraction_duplicate_rc_reads = f"{n_duplicate_rc_reads / n_reads * 100:.2f}%"
    sequence_file = os.path.basename(sequence_file)
    return [
        sequence_file,
        n_reads,
        n_rc_reads,
        n_duplicate_kmers,
        n_duplicate_rc_reads,
        fraction_rc_reads,
        fraction_duplicate_kmers,
        fraction_duplicate_rc_reads,
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input-dir", help="Input directory")
    parser.add_argument("-d", "--delivery", help="Delivery directory for results")
    parser.add_argument("-k", "--kmer-length", help="K-mer length", type=int)
    parser.add_argument("-s", "--subset", help="Subset of reads", type=int)

    args = parser.parse_args()
    delivery = args.delivery
    input_dir = args.input_dir

    # Create output directory
    delivery_analyses_dir = "delivery_analyses"
    delivery_dir = os.path.join(delivery_analyses_dir, delivery)
    output_dir = os.path.join(delivery_dir, "rc-and-duplicates")
    os.makedirs(output_dir, exist_ok=True)

    # Analysis
    k = args.kmer_length
    subset = args.subset

    lines = []
    for sequence_file in os.listdir(input_dir):
        full_path = os.path.join(input_dir, sequence_file)
        line = analyse_fastq(full_path, k, subset)
        lines.append(line)

    df = pd.DataFrame(
        lines,
        columns=[
            "sequence_file",
            "n_reads",
            "n_rc_reads",
            "n_duplicate_kmers",
            "n_duplicate_rc_reads",
            "fraction_rc_reads",
            "fraction_duplicate_kmers",
            "fraction_duplicate_rc_reads",
        ],
    )
    df.to_csv(
        os.path.join(output_dir, "rc_and_duplicates.tsv"),
        sep="\t",
        index=False,
    )
    print(f"Saved to {os.path.join(output_dir, 'rc_and_duplicates.tsv')}")


if __name__ == "__main__":
    main()
