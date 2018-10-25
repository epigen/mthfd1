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

sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


# Peak overlap
df = pd.read_csv(os.path.join("results", "pairwise_jaccard.txt"), sep=' ', index_col=0)
np.fill_diagonal(df.values, np.nan)

order = ["H3K27ac_WT_", "BRD4_WT_", "BRD4_MTHFD1KO_", "BRD4_DMSO_", "BRD4_dBET6_", "MTHFD1_WT_", "MTHFD1_MTHFD1KO_", "MTHFD1_DMSO_", "MTHFD1_dBET6_"]
df = df.loc[order, order]
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
for i, label in enumerate(['DMSO|dBET', 'WT|MTHFD1KO']):
    f = df.index.str.contains(label + '|H3K27ac')
    sns.heatmap(df.loc[f, f], ax=axis[i], square=True, cbar_kws={"label": "Jaccard index"})
plt.tight_layout()
fig.savefig(os.path.join("peak_overlap.jaccard.heatmap.svg"), dpi=300, bbox_inches="tight")


df = pd.read_csv(os.path.join("results", "pairwise_intersect.txt"), sep=' ', index_col=0).astype(float)
np.fill_diagonal(df.values, np.nan)
df = df.loc[order, order]
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
for i, label in enumerate(['DMSO|dBET', 'WT|MTHFD1KO']):
    f = df.index.str.contains(label + '|H3K27ac')
    sns.heatmap(np.log10(1 + df.loc[f, f]), ax=axis[i], square=True, cbar_kws={"label": "Peak count intersect (log10)"})
plt.tight_layout()
fig.savefig(os.path.join("results", "peak_overlap.intersect.heatmap.svg"), dpi=300, bbox_inches="tight")

# as percent of diagonal
df = pd.read_csv(os.path.join("results", "pairwise_intersect.txt"), sep=' ', index_col=0)
df = df.loc[order, order]
df2 = (df / np.diagonal(df) * 100)
df2 = pd.DataFrame(np.triu(df2).T + np.triu(df2), index=df2.index, columns=df2.columns)
df2.to_csv(os.path.join("results", "pairwise_intersect.csv"))
np.fill_diagonal(df2.values, np.nan)
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
for i, label in enumerate(['DMSO|dBET', 'WT|MTHFD1KO']):
    f = df2.index.str.contains(label + '|H3K27ac')
    sns.heatmap(df2.loc[f, f], ax=axis[i], square=True, annot=True, cbar_kws={"label": "Peak overlap (%)"})
plt.tight_layout()
fig.savefig(os.path.join("results", "peak_overlap.intersect.percentage.heatmap.svg"), dpi=300, bbox_inches="tight")

