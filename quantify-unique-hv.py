#!/usr/bin/env python3

import argparse
import subprocess
import os
import pysam
from collections import defaultdict


def is_primary_alignment(read):
    """Check if read is a primary alignment."""
    return not read.is_unmapped and not read.is_secondary and not read.is_supplementary


def extract_reads_from_sam(sam_file, output_fasta):
    print(f"Extracting reads from {sam_file}...")

    # Use samtools fasta to extract primary alignments and convert to fasta
    cmd = [
        "samtools",
        "fasta",
        "-F",
        "0x900",  # Exclude secondary/supplementary alignments
        sam_file,
        ">",
        output_fasta,
    ]

    subprocess.run(" ".join(cmd), shell=True, check=True)


def run_flye_assembly(input_fastq, flye_dir):
    print("Running Flye assembly...")

    cmd = [
        "flye",
        "--nano-corr",  # Assumes <3% error. Filtlong (in mgs-workflow CLEAN) filters for mean 99% accuracy.
        input_fastq,
        "--out-dir",
        flye_dir,
        "--meta",  # Metagenomic mode for viral sequences
    ]

    try:
        subprocess.run(cmd, check=True)

    except subprocess.CalledProcessError as e:
        print(f"\nFlye assembly failed. No contigs constructed.\nError:{e}\n")
        return None


def align_reads_to_assembly(reads_fasta, assembly_fasta, assembly_sam):
    print("Aligning reads to assembly ...")
    if assembly_fasta is None:
        return None

    cmd = [
        "minimap2",
        "-ax",
        "lr:hq",  # Oxford Nanopore read mapping
        "-o",
        assembly_sam,
        assembly_fasta,
        reads_fasta,
    ]

    subprocess.run(cmd, check=True)


def count_aligned_bases(hv_sam, assembly_sam):
    print("Counting aligned bases...")

    ref_coverage = defaultdict(set)
    ref_read_counts = defaultdict(int)  # New counter for reads per reference
    assembly_reads = set()
    hv_lengths = []
    hv_reads = {}

    # Process assembly alignments first
    if assembly_sam is None:
        print("No assembly SAM file found, skipping assembly alignment processing...")
    else:
        print("Processing assembly SAM file...")
        assembly_sam = pysam.AlignmentFile(assembly_sam, "r")
        for read in assembly_sam:
            if is_primary_alignment(read):
                read_name = read.query_name
                assembly_reads.add(read_name)
                reference = read.reference_name
                ref_read_counts[reference] += 1  # Count this read for this reference
                start_pos = read.reference_start
                end_pos = read.reference_end

                for pos in range(start_pos, end_pos):
                    ref_coverage[reference].add(pos)
                hv_reads.append(read_name)
            hv_lengths.append(read.query_length)

    # Process original alignments for reads that didn't align to assembly
    print("Processing original alignments...")
    with pysam.AlignmentFile(hv_sam, "r") as hv_sam:
        for read in hv_sam:
            read_name = read.query_name
            if is_primary_alignment(read):
                reference = read.reference_name
                # Check if reference is one of the specified sequences
                if reference in ["MN306053.1", "KY249660.1", "AY037928.1"]:
                    print(f"{reference} full sequence: {read.seq}")
                ref_read_counts[reference] += 1  # Count this read for this reference
                start_pos = read.reference_start
                end_pos = read.reference_end
                for pos in range(start_pos, end_pos):
                    ref_coverage[reference].add(pos)
                hv_reads[read_name] = read.seq
                hv_lengths.append(read.query_length)

    return ref_coverage, hv_reads, hv_lengths, ref_read_counts


def main():
    parser = argparse.ArgumentParser(
        description="Quantify unique virus base pairs from SAM file"
    )
    parser.add_argument(
        "-i",
        "--input-sam",
        help="Input SAM file containing virus reads",
        required=True,
    )
    parser.add_argument(
        "-d", "--delivery", help="Delivery directory for results", required=True
    )
    args = parser.parse_args()
    delivery = args.delivery
    input_sam = args.input_sam

    # Create output directories
    work_dir = os.path.join("work", delivery)
    output_dir = os.path.join(work_dir, "hv-quant-output")
    flye_dir = os.path.join(output_dir, "flye_assembly")
    os.makedirs(output_dir, exist_ok=True)

    # Create delivery directory for final results
    delivery_analyses_dir = "delivery_analyses"
    delivery_dir = os.path.join(delivery_analyses_dir, delivery)
    os.makedirs(delivery_dir, exist_ok=True)

    # Execute
    reads_fasta = os.path.join(output_dir, "extracted_reads.fasta")
    assembly_fasta = os.path.join(flye_dir, "assembly.fasta")
    assembly_sam = os.path.join(output_dir, "assembly_alignments.sam")

    if os.path.exists(reads_fasta):
        print(f"{output_dir}/reads_fasta already exists, skipping extraction")
    else:
        extract_reads_from_sam(input_sam, reads_fasta)

    if os.path.exists(assembly_fasta):
        print(f"{output_dir}/assembly_fasta already exists, skipping assembly")
    else:
        assembly_fasta = run_flye_assembly(reads_fasta, flye_dir)

    if os.path.exists(assembly_sam):
        print(f"{assembly_sam} already exists, skipping alignment")
    else:
        assembly_sam = align_reads_to_assembly(
            reads_fasta, assembly_fasta, assembly_sam
        )

    ref_coverage, hv_reads, hv_lengths, ref_read_counts = count_aligned_bases(
        input_sam, assembly_sam
    )

    # Calculate metrics
    n_hv_reads = len(hv_reads)
    avg_hv_length = sum(hv_lengths) / len(hv_lengths)
    total_coverage = sum(len(coverage) for coverage in ref_coverage.values())
    reference_matches = len(ref_coverage)

    # Sort references by number of reads
    sorted_refs = sorted(ref_read_counts.items(), key=lambda x: x[1], reverse=True)

    # Print metrics
    print(f"\n###############\n")
    print(f"Delivery: {delivery}")
    print("\n")
    print(f"Number of HV reads: {n_hv_reads}")
    print(f"Average HV read length: {avg_hv_length:.2f}")
    print(f"Total coverage (bp): {total_coverage}")
    print(f"Total references: {reference_matches}")

    # Print top 10 references by read count
    print("\nTop 10 references by number of reads:")
    for genbank_id, read_count in sorted_refs[:10]:
        print(f"{genbank_id}: {read_count} reads")

    # Printing the first 5 HV reads
    print("\n")
    print("First 5 HV reads:")
    for read, seq in list(hv_reads.items())[:5]:
        print(f"{read}: {seq}")
    # Check if assembly failed
    if assembly_fasta is None:
        print(
            "\nWARNING: Assembly failed - no assembly FASTA file was produced. "
            f"Results are only for {input_sam}."
        )


if __name__ == "__main__":
    main()
