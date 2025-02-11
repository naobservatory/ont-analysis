#!/usr/bin/env python3

import os
import csv
import glob
import gzip
import argparse
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pysam

def parse_args():
    parser = argparse.ArgumentParser(description='Filter ONT output to identify human virus sequences')
    parser.add_argument('--delivery-dir', required=True,
                      help='Path to the delivery directory (e.g. NAO-ONT-20250127-Zephyr9)')
    parser.add_argument('--virus-db-index', required=True,
                      help='Path to virus database index file (gzipped TSV)')
    parser.add_argument('--virus-genome-metadata', required=True,
                      help='Path to virus genome metadata file (gzipped TSV)')
    parser.add_argument('--metadata-dir', required=True,
                      help='Directory containing sample metadata files')
    return parser.parse_args()

def get_delivery_paths(delivery_dir):
    """Get standardized paths within the delivery directory structure."""
    paths = {
        'human_alignments': os.path.join(delivery_dir, 'output', 'results', 'human_alignments.sam'),
        'virus_alignments': os.path.join(delivery_dir, 'output', 'results', 'hv_alignments.sam'),
        'fastq_dir': os.path.join(delivery_dir, 'output', 'intermediates', 'reads', 'cleaned'),
        'output_dir': os.path.join(delivery_dir, 'hvs')
    }

    # Create output directory if it doesn't exist
    os.makedirs(paths['output_dir'], exist_ok=True)

    return paths

def main():
    args = parse_args()
    paths = get_delivery_paths(args.delivery_dir)

    print("Starting virus database processing...")
    retain_taxids = set()
    with gzip.open(args.virus_db_index, "rt") as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
            if row["infection_status_human"] != "0":
                retain_taxids.add(int(row["taxid"]))

    retain_genome_ids = set()
    genome_names = {}
    genome_to_taxid = {}

    with gzip.open(args.virus_genome_metadata, "rt") as inf:
        for row in csv.DictReader(inf, delimiter='\t'):
            taxid = int(row["taxid"])
            if taxid in retain_taxids:
                retain_genome_ids.add(row["genome_id"])
                genome_names[row["genome_id"]] = row["organism_name"]
                genome_to_taxid[row["genome_id"]] = taxid


    def build_human_read_ids(human_read_ids_fname):
        print(f"Processing human alignments from {paths['human_alignments']}")
        human_read_ids = set()
        with open(paths['human_alignments']) as inf:
            for line in inf:
                if line.startswith("@"):
                    continue
                bits = line.rstrip("\n").split("\t")

                chromosome = bits[2]
                if chromosome == "*":
                    continue
                assert chromosome.startswith("chr")

                read_id = bits[0]
                human_read_ids.add(read_id)
        with open(human_read_ids_fname, "w") as outf:
            for read_id in sorted(human_read_ids):
                outf.write("%s\n" % read_id)
        print(f"Found {len(human_read_ids)} human reads")

    print(f"Processing virus alignments from {paths['virus_alignments']}")
    read_id_to_genome_id = {}

    with open(paths['virus_alignments']) as inf:
        matching_reads = 0
        for line in inf:
            if line.startswith("@"):
                continue

            bits = line.rstrip("\n").split("\t")
            genome_id = bits[2]

            if genome_id not in retain_genome_ids:
                continue

            read_id = bits[0]
            read_id_to_genome_id[read_id] = genome_id
            matching_reads += 1

        print(f"Found {matching_reads} matching virus reads")

    human_read_ids_fname = os.path.join(paths['output_dir'], "human.read_ids.txt")
    if not os.path.exists(human_read_ids_fname):
        build_human_read_ids(human_read_ids_fname)

    human_read_ids = set()
    with open(human_read_ids_fname) as inf:
        for line in inf:
            human_read_ids.add(line.strip())

    hv_read_ids_fname = os.path.join(paths['output_dir'], "hv.read_ids.txt")
    with open(hv_read_ids_fname, "w") as outf:
        for read_id in sorted(read_id_to_genome_id):
            if read_id in human_read_ids:
                continue
            outf.write("%s\n" % read_id)

    output = []

    # Process each sample's metadata
    for metadata_file in glob.glob(os.path.join(args.metadata_dir, "*.metadata.tsv")):
        with open(metadata_file) as inf:
            for row in csv.DictReader(inf, delimiter='\t'):
                sample = row["sample"]
                loc = row["fine_location"]
                date = row["date"]

                # Look for FASTQ files without the _human suffix
                sample_pattern = sample.replace("_human", "") + "-div*.fastq.gz"
                fastq_files = glob.glob(os.path.join(paths['fastq_dir'], sample_pattern))
                print(f"Found {len(fastq_files)} FASTQ files for sample {sample}")

                for fastq in fastq_files:
                    print(f"Processing FASTQ file: {fastq}")
                    fastq_subset_out = os.path.join(paths['output_dir'],
                                                  os.path.basename(fastq).replace(".fastq.gz", ".hv.fastq"))
                    if not os.path.exists(fastq_subset_out):
                        tmp_fastq_subset_out = fastq_subset_out + ".tmp"
                        with open(tmp_fastq_subset_out, "w") as outf:
                            subprocess.check_call(
                                ["seqtk", "subseq", fastq, hv_read_ids_fname],
                                stdout=outf)
                        os.rename(tmp_fastq_subset_out, fastq_subset_out)

                    with open(fastq_subset_out) as inf2:
                        for read_id, seq, qua in FastqGeneralIterator(inf2):
                            if read_id in read_id_to_genome_id:
                                genome_id = read_id_to_genome_id[read_id]
                                output.append((
                                    genome_to_taxid[genome_id],
                                    genome_names[genome_id],
                                    genome_id,
                                    sample,
                                    loc,
                                    date,
                                    read_id,
                                    seq,
                                    qua,
                                ))

    print(f"Found {len(output)} total virus reads")
    if output:
        print("Writing output file...")
        final_read_ids = set()
        with gzip.open(os.path.join(paths['output_dir'], "analysis.tsv.gz"), "wt") as outf:
            outf.write("\t".join((
                "taxid",
                "organism_name",
                "genome_id",
                "sample",
                "location",
                "date",
                "read_id",
                "sequence",
                "quality")) + "\n")
            for record in sorted(output):
                final_read_ids.add(record[6])  # Add read_id to set
                outf.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % record)
        print("Done!")
    else:
        print("No output to write")

if __name__ == "__main__":
    main()