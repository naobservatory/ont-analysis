#!/usr/bin/env python3
import seaborn as sns
import json
import matplotlib.pyplot as plt
import numpy as np
import argparse
import pandas as pd

def generate_label_dict(metadata):
    metadata["date_treatment"] = metadata["date"] + " | " + metadata["notes"]
    barcode_to_date_treatment = dict(zip(metadata["sample"], metadata["date_treatment"]))
    return barcode_to_date_treatment

def plot_length_distributions(data, output_file, barcode_to_label):
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


        if 'Amicon' in barcode_to_label[barcode] or 'CP' in barcode_to_label[barcode]:
            color = sns.color_palette()[0]  # Seaborn blue
        else:
            color = sns.color_palette()[1]  # Seaborn orange

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
            color=color
        )

        # Set log scale for x-axis
        ax.set_xscale("log")

        # Set title and labels
        ax.set_title(barcode_to_label[barcode])
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
    parser.add_argument("-i", "--input_json", required=True, help="Input JSON file with read length data")
    parser.add_argument("-m", "--input_metadata", required=True, help="Input metadata file (e.g., metadata.tsv)")

    args = parser.parse_args()

    # Load the JSON data
    with open(args.input_json) as f:
        data = json.load(f)

    # Load the metadata
    metadata = pd.read_csv(args.input_metadata, sep="\t")
    barcode_to_label = generate_label_dict(metadata)

    # Drop samples with no notes


    # Create the plot
    plot_name = args.input_json.replace(".json", ".png")
    plot_length_distributions(data, plot_name, barcode_to_label)
    print(f"Plot saved to {plot_name}")


if __name__ == "__main__":
    main()
