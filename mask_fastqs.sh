#!/bin/bash
set -e # Exit on any error

# Check if input and output directories are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <input_dir> <output_dir>"
    echo "The input directory should contain .fastq or .fastq.gz files"
    exit 1
fi

input_dir=$1
delivery=$2

# Create temporary directory for intermediate files
tmp_dir="tmp_dustmasker_processing"
mkdir -p "$tmp_dir"

# Create output directory and set output path
mkdir -p "work/$delivery"
output_dir="work/$delivery"

# Process each FASTQ file
for fastq in "$input_dir"/*.fastq "$input_dir"/*.fastq.gz; do
    if [ ! -e "$fastq" ]; then
        continue
    fi

    echo "Processing $fastq..."
    base_name=$(basename "$fastq" .fastq.gz)
    base_name=$(basename "$base_name" .fastq)

    # Decompress if gzipped
    if [[ "$fastq" == *.gz ]]; then
        echo "Decompressing gzipped file..."
        gunzip -c "$fastq" > "${tmp_dir}/${base_name}.fastq"
    else
        cp "$fastq" "${tmp_dir}/${base_name}.fastq"
    fi

    # Convert FASTQ to FASTA using seqtk
    echo "Converting FASTQ to FASTA..."
    seqtk seq -A "${tmp_dir}/${base_name}.fastq" > "${tmp_dir}/${base_name}.fasta"

    # Run dustmasker
    echo "Running dustmasker..."
    dustmasker -in "${tmp_dir}/${base_name}.fasta" \
        -out "${tmp_dir}/${base_name}_masked.fasta" \
        -outfmt fasta

    # Filter out masked sequences
    echo "Filtering masked sequences..."
    ./bin/fasta-drop-non-complex.py \
        "${tmp_dir}/${base_name}_masked.fasta" \
        "${output_dir}/${base_name}_complex.fasta"
done

# Clean up
rm -r "$tmp_dir"

echo "Processing complete. Complex reads are in ${output_dir}/"