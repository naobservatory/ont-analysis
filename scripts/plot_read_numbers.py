#!/usr/bin/env python3

import argparse
import os
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def main():
    parser = argparse.ArgumentParser(description="Plot read numbers")
    parser.add_argument("-d", "--delivery", required=True, help="Delivery name")
    args = parser.parse_args()
    delivery = args.delivery

    delivery_metadata = pd.read_csv(
        f"../mgs-restricted/dashboard/{delivery}.metadata.tsv", sep="\t"
    )

    delivery_metadata = delivery_metadata[delivery_metadata["demultiplexed"] == True]
    delivery_metadata["date_treatment"] = (
        delivery_metadata["date"] + "|" + delivery_metadata["notes"]
    )

    # Create figure and axis
    fig, ax = plt.subplots(dpi=300, figsize=(10, 4))

    # Create horizontal bar chart
    sns.barplot(
        data=delivery_metadata,
        y="date_treatment",
        x="reads",
        ax=ax,
        zorder=10,
    )

    # Customize appearance
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, axis="x", linestyle="--", alpha=0.7, zorder=-5)

    # Format x-axis to use comma separator for thousands
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ",")))

    ax.set_xlabel("Total number of reads")
    ax.set_ylabel("")

    plt.tight_layout()

    output_dir = f"delivery_analyses/{delivery}"
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(f"{output_dir}/read_numbers.png")


if __name__ == "__main__":
    main()
