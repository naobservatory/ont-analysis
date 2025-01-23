#!/usr/bin/env python3

import json
import matplotlib.pyplot as plt
import numpy as np
import argparse


def plot_length_distributions(data, output_file):
    """Create a multi-panel plot of read length distributions."""
    # Get number of samples
    n_samples = len(data)

    # Create a figure with subplots
    fig, axes = plt.subplots(n_samples, 1, figsize=(10, 2 * n_samples))
    if n_samples == 1:
        axes = [axes]

    # Set up log-scale bins
    bins = np.logspace(0, 4, 100)  # 10^0 to 10^4, 100 bins

    for ax, (barcode, lengths) in zip(axes, data.items()):
        # Convert length dictionary to arrays
        lengths_array = np.array(list(lengths.keys()), dtype=int)
        counts_array = np.array(list(lengths.values()))

        # Create histogram
        ax.hist(
            lengths_array,
            bins=bins,
            weights=counts_array,
            alpha=0.75,
            edgecolor="black",
            linewidth=0.5,
        )

        # Set log scale for x-axis
        ax.set_xscale("log")

        # Set title and labels
        ax.set_title(barcode)
        ax.set_ylabel("Count")

        # Only show x-label for bottom plot
        if ax == axes[-1]:
            ax.set_xlabel("Read length")

    # Adjust layout to prevent overlap
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
