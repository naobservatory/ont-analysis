#!/usr/bin/env python3

import json
import matplotlib.pyplot as plt
import numpy as np
import argparse


def plot_length_distributions(data, output_file):
    """Create a single plot of read length distributions across all files."""
    # Create a figure
    fig, ax = plt.subplots(figsize=(10, 4))

    # Set up linear bins
    bins = np.logspace(0, 4, 100)  # 10^0 to 10^4, 100 bins

    # Combine all lengths into single arrays
    all_lengths = []
    all_counts = []

    for file, lengths in data.items():
        lengths_array = np.array(list(lengths.keys()), dtype=int)
        counts_array = np.array(list(lengths.values()))
        all_lengths.extend(lengths_array)
        all_counts.extend(counts_array)

    # Create histogram
    ax.hist(
        all_lengths,
        bins=bins,
        weights=all_counts,
        alpha=0.75,
        edgecolor="black",
        linewidth=0.5,
    )

    ax.set_xlim(50, 10000)
    # Set log scale for x-axis
    # ax.set_xscale("log")

    # Set title and labels
    ax.set_title("Read Length Distribution Across All Files")
    ax.set_ylabel("Count")
    ax.set_xlabel("Read length")

    # Adjust layout
    plt.tight_layout()

    # Save the plot
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot read length distributions from JSON file"
    )
    parser.add_argument("input_json", help="Input JSON file with read length data")
    parser.add_argument("output_plot", help="Output plot file (e.g., plot.png)")

    args = parser.parse_args()

    # Load the JSON data
    with open(args.input_json) as f:
        data = json.load(f)

    # Create the plot
    plot_length_distributions(data, args.output_plot)
    print(f"Plot saved to {args.output_plot}")


if __name__ == "__main__":
    main()
