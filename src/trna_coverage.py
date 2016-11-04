
import os
import numpy as np
import pandas as pd
import pysam

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import cPickle as pickle

# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def genebody_coverage(gtf, bam_file, left=1000, right=1000, genebody_bins=100, boundary_resolution=1):
    """
    """
    bam = pysam.AlignmentFile(bam_file)
    cov = np.zeros((gtf.shape[0], left + genebody_bins + right))
    cov.fill(np.nan)

    for g, gene in enumerate(gtf.index):
        s = gtf.ix[gene].copy()
        # Get coverage at left boundary
        # handle strand
        if s['strand'] == "+":
            boundary_start = s['start'] - left
            boundary_end = s['start']
        if s['strand'] == "-":
            boundary_start = s['end']
            boundary_end = s['end'] + left
        c = np.array(bam.count_coverage(reference=s['chrom'], start=boundary_start, end=boundary_end)).sum(0)
        # interpolate to specified resolution
        size = abs(boundary_end - boundary_start)
        boundary_bins = size / float(boundary_resolution)  # for boundaries
        cov[g, :left] = np.interp(np.linspace(0, len(c), boundary_bins), range(len(c)), c)

        # Get coverage in genebody
        c = np.array(bam.count_coverage(reference=s['chrom'], start=s['start'], end=s['end'])).sum(0)
        # handle strand
        if s['strand'] == "-":
            c = c[::-1]
        # interpolate to specified bins
        cov[g, left:-right] = np.interp(np.linspace(0, len(c), genebody_bins), range(len(c)), c)

        # Get coverage at right boundary
        # handle strand
        if s['strand'] == "+":
            boundary_start = s['end']
            boundary_end = s['end'] + right
        elif s['strand'] == "-":
            boundary_start = s['start'] - right
            boundary_end = s['start']
        c = np.array(bam.count_coverage(reference=s['chrom'], start=boundary_start, end=boundary_end)).sum(0)
        # interpolate to specified resolution
        size = abs(boundary_end - boundary_start)
        boundary_bins = size / float(boundary_resolution)  # for boundaries
        cov[g, -right:] = np.interp(np.linspace(0, len(c), boundary_bins), range(len(c)), c)

    return cov


gtf_file = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.tRNAs.gtf.gz"
gtf = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, compression="gzip")
gtf[8] = pd.Series(gtf[8].str.split(";")).apply(lambda x: x[0].split(" ")[1])
gtf = gtf[[0, 3, 4, 6, 8]]
gtf.columns = ['chrom', 'start', 'end', 'strand', 'id']

samples = [
    "HAP1_ChIPmentation_H3K27ac_WT",
    "HAP1_ChIPmentation_BRD4_WT",
    "HAP1_ChIPmentation_MTHFD1_C3_WT",
    "HAP1_ChIPmentation_IgG_WT"
]

cov = dict()
for i, sample in enumerate(samples):
    # Get sample
    print(sample)
    bam_file = "/home/arendeiro/projects/mthfd1/data/{}/mapped/{}.trimmed.bowtie2.filtered.bam".format(sample, sample)

    # Get coverage
    cov[sample] = genebody_coverage(gtf, bam_file, left=3000, right=3000)
pickle.dump(cov, os.path.join("results", "tRNA_coverage.pickle"))

# Plot
colors = ["red", "brown", "green", "grey"]
fig, axis = plt.subplots(2, 2, sharex=True, sharey=False)
axis = axis.flatten()
for i, sample in enumerate(samples):
    # Normalized transform
    cov2 = pd.DataFrame(cov[sample] + 1).apply(lambda x: (x) / float(x.sum()), axis=1)

    # Plot
    axis[i].plot(range(cov2.shape[1]), cov2.mean(0).rolling(10).mean(), label=sample, linewidth=1.5, alpha=0.8)
    axis[i].fill_between(
        range(cov2.shape[1]),
        np.nanpercentile(cov2, 25, axis=0),
        np.nanpercentile(cov2, 75, axis=0),
        alpha=0.25,
        color=colors[i],
        linewidth=0.5,
        label=sample)
    axis[i].set_xlim((0 + 1000, cov2.shape[1] - 1000))
    axis[i].set_title(sample)

    # add annotations
    axis[i].axvline(3000, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
    axis[i].axvline(3000 + 100, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
sns.despine(fig)
fig.savefig(os.path.join("results", "tRNA_coverage.svg"), bbox_inches="tight")


fig, axis = plt.subplots(1)
for i, sample in enumerate(samples):
    # Normalized transform
    cov2 = pd.DataFrame(cov[sample] + 1).apply(lambda x: (x) / float(x.sum()), axis=1)

    # Plot
    axis.plot(range(cov2.shape[1]), cov2.mean(0).rolling(10).mean(), label=sample, linewidth=0.8, alpha=0.8)

    # add annotations
    axis.axvline(3000, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
    axis.axvline(3000 + 100, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
axis.set_xlim((0 + 2000, cov2.shape[1] - 2000))
axis.legend()
sns.despine(fig)
fig.savefig(os.path.join("results", "tRNA_coverage.stacked.svg"), bbox_inches="tight")


fig, axis = plt.subplots(1)
for i, sample in enumerate(samples):
    # Normalized transform
    cov2 = pd.DataFrame(cov[sample] + 1).apply(lambda x: (x) / float(x.sum()), axis=1)

    # Plot
    axis.plot(range(cov2.shape[1]), cov2.mean(0).rolling(50).mean(), label=sample, linewidth=0.8, alpha=0.8)

    # add annotations
    axis.axvline(3000, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
    axis.axvline(3000 + 100, linewidth=0.5, color="black", linestyle="--", alpha=0.5)
axis.set_xlim((0 + 2000, cov2.shape[1] - 2000))
axis.legend()
sns.despine(fig)
fig.savefig(os.path.join("results", "tRNA_coverage.stacked.smooth50.svg"), bbox_inches="tight")


# Heatmaps
fig, axis = plt.subplots(1, 4, figsize=(16, 24))
axis = axis.flatten()
for i, sample in enumerate(samples):
    print(sample)

    # Subset
    cov2 = cov[sample][:, range(0 + 2000, cov[sample].shape[1] - 2000)]

    # Normalized transform
    cov2 = pd.DataFrame(cov2 + 1).apply(lambda x: (x) / float(x.sum()), axis=1)

    # Log transform
    cov2 = np.log2(1 + (cov2 * 1e4))

    # Sort by mean
    if i == 0:
        index = cov2.mean(axis=1).sort_values(ascending=True).index
    cov2 = cov2.ix[index]

    # plot
    sns.heatmap(
        cov2,
        xticklabels=False,
        yticklabels=False,
        cmap="Blues",
        vmin=2,
        vmax=7,
        ax=axis[i]
    )
    axis[i].set_title(sample)
fig.savefig(os.path.join("results", "tRNA_coverage.heatmap.png"), bbox_inches="tight", dpi=100)
fig.savefig(os.path.join("results", "tRNA_coverage.heatmap.300dpi.png"), bbox_inches="tight", dpi=300)
fig.savefig(os.path.join("results", "tRNA_coverage.heatmap.svg"), bbox_inches="tight")
