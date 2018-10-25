#!/usr/bin/env python

"""
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from peppy import Project
from ngs_toolkit.chipseq import ChIPSeqAnalysis


sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


def main():
    # Start Project object and curate samples
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    prj._samples = [s for s in prj.samples if s.to_use in [1, "1", True]]
    for sample in prj.samples:
        if hasattr(sample, "protocol"):
            sample.library = sample.protocol

        if sample.library in ["ATAC-seq", "ChIP-seq", "ChIPmentation"]:
            sample.mapped = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(
                sample.paths.sample_root,
                "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(
                sample.paths.sample_root,
                "peaks", sample.name + "_peaks.narrowPeak")

    # Start Analysis object
    chip_analysis = ChIPSeqAnalysis(name=prj.project_name + "_chipseq", samples=prj.samples)
    chip_analysis.sample_variables = prj.sample_attributes
    chip_analysis.group_variables = prj.group_attributes

    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))

    chip_analysis.call_peaks_from_comparisons(comparison_table)

    chip_analysis.filter_peaks(comparison_table, filter_bed=os.path.join("data", "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))
    # filter_peaks(chip_analysis, comparison_table, filter_bed=os.path.join("data", "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))

    # Get summary of peak calls
    peak_counts = chip_analysis.summarize_peaks_from_comparisons(comparison_table, filtered=True)
    # peak_counts = summarize_peaks_from_comparisons(chip_analysis, comparison_table, filtered=True)
    peak_counts.to_csv(os.path.join(chip_analysis.results_dir, "chipseq_peaks", "peak_count_summary.csv"), index=False)


    # Plot peak numbers
    peak_counts['IP'] = list(map(lambda x: x[1], peak_counts['comparison_name'].str.split("_")))
    peak_counts['condition'] = list(map(lambda x: x[2], peak_counts['comparison_name'].str.split("_")))
    ips = ["BRD4", "MTHFD1"]

    peak_counts = peak_counts.loc[
        (((peak_counts['IP'] == "BRD4") & (peak_counts['peak_type'] == "homer_factor")) |
        ((peak_counts['IP'] == "H3K27ac") & (peak_counts['peak_type'] == "homer_factor")) |
        ((peak_counts['IP'] == "MTHFD1") & (peak_counts['peak_type'] == "homer_histone"))) &
        (peak_counts['comparison_name'].str.contains("IgG"))
    ]

    fig, axis = plt.subplots(2, len(ips), figsize=(len(ips) * 2, 2 * 2), sharey=False)
    for i, label in enumerate(['DMSO|dBET', 'WT|MTHFD1KO']):
        for j, ip in enumerate(ips):
            sns.barplot(
                data=peak_counts[(peak_counts['IP'] == ip) & (peak_counts['comparison_name'].str.contains(label))],
                x='condition', y='peak_counts', ax=axis[i, j])
            axis[i, j].set_title(ip)
            axis[i, j].set_xlabel("Condition")
            axis[i, j].set_ylabel("Peak count")
    sns.despine(fig)
    fig.savefig(os.path.join("peak_count_summary.barplot.svg"), dpi=300, bbox_inches="tight")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
