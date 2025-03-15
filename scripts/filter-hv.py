#!/usr/bin/env python3

import pandas as pd
import os
import csv
from collections import namedtuple
import gzip
import argparse

# Define named tuple for blast matches
BlastMatch = namedtuple("BlastMatch", ["taxid", "genome_id", "bitscore"])

# Constants
VIRUS_TAXID = 10239


def parse_args():
    parser = argparse.ArgumentParser(description="Process NGS data for viral analysis")
    parser.add_argument(
        "-d",
        "--delivery",
        type=str,
        required=True,
        help="Delivery name (e.g. NAO-ONT-20250313-Zephyr12)",
    )
    parser.add_argument(
        "-i",
        "--index-date",
        type=str,
        help="Index date (e.g. 20250224 for index.20250224.taxonomy-names.dmp)",
    )
    return parser.parse_args()


def setup_directories(delivery):
    delivery_dir = os.path.join("deliveries", delivery)
    output_dir = os.path.join(delivery_dir, "output")
    results_dir = os.path.join(output_dir, "results")
    delivery_analyses_dir = os.path.join("delivery_analyses", delivery)

    os.makedirs(delivery_analyses_dir, exist_ok=True)

    return {
        "delivery_dir": delivery_dir,
        "output_dir": output_dir,
        "results_dir": results_dir,
        "delivery_analyses_dir": delivery_analyses_dir,
        "hv_fname": os.path.join(results_dir, "hv.tsv.gz"),
        "blast_fname": os.path.join(results_dir, "merged_blast_filtered.tsv.gz"),
    }


def load_taxonomy(index_date):
    parents = {}
    taxid_names = {}

    # Load parent-child relationships
    with open(os.path.join(f"index/index.{index_date}.taxonomy-nodes.dmp")) as inf:
        for line in inf:
            child_taxid, parent_taxid, rank, *_ = line.replace("\t|\n", "").split(
                "\t|\t"
            )
            parents[int(child_taxid)] = int(parent_taxid)

    # Load taxid names
    with open(os.path.join(f"index/index.{index_date}.taxonomy-names.dmp")) as inf:
        for line in inf:
            taxid, name, unique_name, name_class = line.replace("\t|\n", "").split(
                "\t|\t"
            )
            taxid = int(taxid)
            if taxid not in taxid_names or name_class == "scientific name":
                taxid_names[taxid] = name

    return parents, taxid_names


def descends_from(taxid, target, parents):
    if taxid in [0, 1] or taxid not in parents:
        return False
    if taxid == target:
        return True
    return descends_from(parents[taxid], target, parents)


def is_viral(taxid, parents, viral_taxids_cache=None):
    if viral_taxids_cache is None:
        viral_taxids_cache = {}

    if taxid not in viral_taxids_cache:
        viral_taxids_cache[taxid] = descends_from(taxid, VIRUS_TAXID, parents)

    return viral_taxids_cache[taxid]


def load_human_virus_data(index_date):
    hv_taxids = set()
    genome_to_taxid = {}

    # Load human-infecting virus taxids
    with gzip.open(
        os.path.join(f"index/index.{index_date}.total-virus-db-annotated.tsv.gz"), "rt"
    ) as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        for row in reader:
            if row["infection_status_human"] != "0":
                hv_taxids.add(int(row["taxid"]))

    # Load genome-to-taxid mapping for human viruses
    with gzip.open(
        os.path.join(f"index/index.{index_date}.virus-genome-metadata-gid.tsv.gz"), "rt"
    ) as inf:
        reader = csv.DictReader(inf, delimiter="\t")
        for row in reader:
            taxid = int(row["taxid"])
            if taxid in hv_taxids:
                genome_to_taxid[row["genome_id"]] = taxid

    return hv_taxids, genome_to_taxid


def process_blast_matches(blast_fname, hv_taxids):
    best_matches = {}
    best_hv_matches = {}
    viral_taxids_cache = {}

    with gzip.open(blast_fname, "rt") as inf:
        # Use CSV reader for consistent parsing
        reader = csv.reader(inf, delimiter="\t")
        header = next(reader)  # Skip header

        for bits in reader:
            try:
                (
                    qseqid,
                    sseqid,
                    sgi,
                    staxid,
                    qlen,
                    evalue,
                    bitscore,
                    qcovs,
                    length,
                    pident,
                    mistmatch,
                    gapopen,
                    ssstrand,
                    qstart,
                    qend,
                    sstart,
                    send,
                    bitscore_rank_dense,
                    bitscore_fraction,
                ) = bits

                read_id = qseqid
                _, _, _, genome_id, _ = sseqid.split("|")
                taxid = int(staxid)
                bitscore = float(bitscore)

                bm = BlastMatch(
                    taxid=taxid,
                    genome_id=genome_id,
                    bitscore=bitscore,
                )

                # Track best overall match
                if (
                    read_id not in best_matches
                    or best_matches[read_id].bitscore < bm.bitscore
                ):
                    best_matches[read_id] = bm

                # Track best human virus match
                if taxid in hv_taxids:
                    if (
                        read_id not in best_hv_matches
                        or best_hv_matches[read_id].bitscore < bm.bitscore
                    ):
                        best_hv_matches[read_id] = bm

            except Exception as e:
                print(f"Error processing line: {bits}")
                raise e

    return best_matches, best_hv_matches


def process_hv_file(
    hv_fname,
    best_matches,
    best_hv_matches,
    hv_taxids,
    parents,
    taxid_names,
    delivery_analyses_dir,
    delivery,
):
    out = []
    excluded_lines = []

    # Stats counters
    stats = {
        "total_reads": 0,
        "no_hv_match": 0,
        "non_viral_match": 0,
        "better_non_hv_match": 0,
        "no_judgement": 0,
        "kept_reads": 0,
    }

    no_judge_reads = []
    viral_taxids_cache = {}

    with gzip.open(hv_fname, "rt") as inf:
        reader = csv.DictReader(inf, delimiter="\t")

        # Prepare extended header
        header = reader.fieldnames.copy()
        header.extend(
            [
                "minimap2_taxid_name",
                "blast_taxid",
                "blast_taxid_name",
                "blast_genome_id",
                "blast_bitscore",
            ]
        )

        out.append(header)
        excluded_lines.append(header)

        for line in reader:
            stats["total_reads"] += 1

            # Get original fields
            bits = [line[h] for h in reader.fieldnames]
            read_id = line["query_name"]
            minimap2_taxid = int(line["minimap2_taxid_primary"])
            minimap2_taxid_name = taxid_names.get(minimap2_taxid, "Unknown")
            bits.append(minimap2_taxid_name)

            general_bm = best_matches.get(read_id)
            hv_bm = best_hv_matches.get(read_id)

            # Case 1: No human virus match
            if not hv_bm:
                stats["no_hv_match"] += 1
                if not general_bm:
                    bits.extend(
                        [
                            "No BLAST match",
                            "No BLAST match",
                            "No BLAST match",
                            "No BLAST match",
                        ]
                    )
                else:
                    taxid = general_bm.taxid
                    genome_id = general_bm.genome_id
                    bitscore = general_bm.bitscore
                    bits.extend(
                        [taxid, taxid_names.get(taxid, "Unknown"), genome_id, bitscore]
                    )
                excluded_lines.append(bits)
                continue

            # Case 2: Best match is non-viral
            if not is_viral(general_bm.taxid, parents, viral_taxids_cache):
                stats["non_viral_match"] += 1
                continue

            # Case 3: Best match is much better than best human virus match
            if general_bm.bitscore > hv_bm.bitscore * 2:
                stats["better_non_hv_match"] += 1
                continue

            # Case 4: Taxid not in human virus set
            taxid = hv_bm.taxid
            if taxid not in hv_taxids:
                stats["no_judgement"] += 1
                no_judge_reads.append([read_id, taxid, hv_bm.genome_id, hv_bm.bitscore])
                continue

            # Case 5: Valid human virus match
            genome_id = hv_bm.genome_id
            bits.extend([taxid, taxid_names[taxid], genome_id, hv_bm.bitscore])
            stats["kept_reads"] += 1
            out.append(bits)

    # Write output files
    with open(os.path.join(delivery_analyses_dir, f"{delivery}.hv.tsv"), "wt") as outf:
        for bits in out:
            outf.write("\t".join(str(x) for x in bits) + "\n")

    with open(
        os.path.join(delivery_analyses_dir, f"{delivery}-excluded.tsv"), "wt"
    ) as outf:
        for bits in excluded_lines:
            outf.write("\t".join(str(x) for x in bits) + "\n")

    return stats, no_judge_reads


def print_stats(stats, no_judge_reads):
    print(f"Total reads processed: {stats['total_reads']}")
    print(f"Reads without HV match: {stats['no_hv_match']}")
    print(f"Reads with non-viral best match: {stats['non_viral_match']}")
    print(f"Reads with much better non-HV match: {stats['better_non_hv_match']}")
    print(f"Reads kept: {stats['kept_reads']}")
    print(f"Reads with no judgement: {stats['no_judgement']}")

    # Print a sample of no-judgment reads
    reads_printed = 0
    for read in no_judge_reads:
        if reads_printed < 10:
            print(read)
            reads_printed += 1
        else:
            break


def main():
    args = parse_args()

    # Setup directories and file paths
    paths = setup_directories(args.delivery)

    # Load taxonomy and human virus data
    parents, taxid_names = load_taxonomy(args.index_date)
    hv_taxids, genome_to_taxid = load_human_virus_data(args.index_date)

    # Process blast matches
    best_matches, best_hv_matches = process_blast_matches(
        paths["blast_fname"], hv_taxids
    )

    # Process HV file
    stats, no_judge_reads = process_hv_file(
        paths["hv_fname"],
        best_matches,
        best_hv_matches,
        hv_taxids,
        parents,
        taxid_names,
        paths["delivery_analyses_dir"],
        args.delivery,
    )

    # Print statistics
    print_stats(stats, no_judge_reads)


if __name__ == "__main__":
    main()
