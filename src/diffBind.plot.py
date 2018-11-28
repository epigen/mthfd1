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
from scipy.stats import mannwhitneyu


sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


output_dir = os.path.join("results", "differential_analysis_ChIP-seq")

# Read in DiffBind results
comparisons = {
    "BRD4_dBET": "DBA_BRD4_dBET_vs_BRD4_DMSO.csv",
    "MTHFD1_dBET": "DBA_MTHFD1_dBET_vs_MTHFD1_DMSO.csv",
    "BRD4_MTHFD1KO": "DBA_BRD4_MTHFD1KO_vs_BRD4_WT.csv",
    "MTHFD1_MTHFD1KO": "DBA_MTHFD1_MTHFD1KO_vs_MTHFD1_WT.csv"
}
res = pd.DataFrame()
for comp, f in comparisons.items():
    df = pd.read_csv(os.path.join(output_dir, f))
    df['comparison_name'] = comp
    res = res.append(df)

    # Plot diagnostic plots
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.axvline(0, linestyle="--", alpha=0.5, color="black")
    axis.scatter(df['Fold'], -np.log10(df['p-value']), s=0.5, alpha=0.1, rasterized=True, color="black")
    s = df[(df['FDR'] < 0.1) & (df['Fold'].abs() >= 1)].index
    axis.scatter(df.loc[s, 'Fold'], -np.log10(df.loc[s, 'p-value']), s=2, alpha=0.2, rasterized=True, color="red")
    l = df['Fold'].abs().max() * 1.05
    axis.set_xlim((-l, l))
    p = -np.log10(df.loc[:, 'p-value'])
    if (p > 4).sum() < 3:
        axis.set_ylim((0, 4))
    axis.set_xlabel("log2(Fold-Change)")
    axis.set_ylabel("-log10(p-value)")
    fig.savefig(os.path.join(output_dir, "diffBind.volcano.{}.scatter.svg".format(comp)), dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(4, 4))
    m = df.iloc[:, [4, 5]].max(1).max()
    axis.plot((0, m), (0, m), linestyle="--", alpha=0.5, color="black")
    axis.scatter(df.iloc[:, 4], df.iloc[:, 5], s=0.5, alpha=0.1, rasterized=True, color="black")
    s = df[(df['FDR'] < 0.1) & (df['Fold'].abs() >= 1)].index
    axis.scatter(df.iloc[s, 4], df.iloc[s, 5], s=2, alpha=0.2, rasterized=True, color="red")
    axis.set_xlabel(df.iloc[:, 4].name)
    axis.set_ylabel(df.iloc[:, 5].name)
    fig.savefig(os.path.join(output_dir, "diffBind.scatterplot.{}.scatter.svg".format(comp)), dpi=300, bbox_inches="tight")

res = res[~res['Chr'].str.contains("_")]
res = res.set_index(['Chr', 'Start', 'End'])

res.to_csv(os.path.join(output_dir, "diffBind.results.csv"))


# fold changes
res2 = res.pivot_table(index=['Chr', 'Start', 'End'], columns='comparison_name', values="Fold")
res3 = res2.loc[res[(res['FDR'] <= 0.1) & (res['Fold'].abs() >= 1)].index]

grid = sns.clustermap(
    res3, cmap="RdBu_r", center=0, metric="correlation",
    rasterized=True, yticklabels=False, cbar_kws={"label": "log2(Fold-change)"})
grid.savefig(os.path.join(
        output_dir,
        "diffBind.log_fold_changes.signficant.heatmap.svg"),
    dpi=300, bbox_inches="tight")

# signed p-values
res4 = res.pivot_table(
    index=['Chr', 'Start', 'End'], columns='comparison_name', values="p-value")
res5 = -np.log10(res4.loc[res[(res['FDR'] <= 0.1) & (res['Fold'].abs() >= 1)].index]) * (res3 >= 0).astype(int).replace(0, -1)

grid2 = sns.clustermap(
    res5, cmap="RdBu_r", center=0, metric="correlation",
    rasterized=True, yticklabels=False, cbar_kws={"label": "Signed log(p-value)"},
    col_linkage=grid.dendrogram_col.linkage, row_linkage=grid.dendrogram_row.linkage)
grid2.savefig(os.path.join(
        output_dir,
        "diffBind.pvalue.signficant.heatmap.svg"),
    dpi=300, bbox_inches="tight")

# concentration values
res6 = (
    res
    .reset_index()
    .melt(id_vars=['Chr', 'Start', 'End'] + ['Conc', "FDR", "p-value", "Fold", "comparison_name"])
    .pivot_table(index=['Chr', 'Start', 'End'], columns=["variable"], values="value"))
res7 = res6.loc[res[(res['FDR'] <= 0.1) & (res['Fold'].abs() >= 1)].index]

grid3 = sns.clustermap(
    res7, metric="correlation", figsize=(4, 4),
    rasterized=True, yticklabels=False, cbar_kws={"label": "Concentration"},
    # row_linkage=grid.dendrogram_row.linkage)
    )
grid3.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.heatmap.svg"),
    dpi=300, bbox_inches="tight")

order = [
    "Conc_BRD4_WT", "Conc_BRD4_MTHFD1KO", "Conc_BRD4_DMSO", "Conc_BRD4_dBET",
    "Conc_MTHFD1_WT", "Conc_MTHFD1_MTHFD1KO", "Conc_MTHFD1_DMSO", "Conc_MTHFD1_dBET"]
p = res7.loc[:, order]
# p = pd.DataFrame(scipy.stats.zscore(scipy.stats.zscore(p, 0), 1), index=p.index, columns=p.columns)
grid3 = sns.clustermap(
    p, metric="correlation", figsize=(4, 4),
    rasterized=True, yticklabels=False, cbar_kws={"label": "Concentration"},
    col_cluster=False,
    )
grid3.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.heatmap.sorted.svg"),
    dpi=300, bbox_inches="tight")
grid3 = sns.clustermap(
    p.loc[:, p.columns.str.contains("DMSO|dBET")], metric="correlation", figsize=(4, 4),
    rasterized=True, yticklabels=False, cbar_kws={"label": "Concentration"},
    col_cluster=False, robust=True)
grid3.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.heatmap.sorted.dmso_dbet.svg"),
    dpi=300, bbox_inches="tight")
grid3 = sns.clustermap(
    p.loc[:, ~p.columns.str.contains("DMSO|dBET")], metric="correlation", figsize=(4, 4),
    rasterized=True, yticklabels=False, cbar_kws={"label": "Concentration"},
    col_cluster=False, robust=True)
grid3.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.heatmap.sorted.wt_ko.svg"),
    dpi=300, bbox_inches="tight")


# res7 = res6.loc[res[res['comparison_name'].str.contains("BRD4_dBET")].sort_values('p-value').head(1000).index, order]
# p = res7.loc[:, order]
grid3 = sns.clustermap(
    p.loc[:, p.columns.str.contains("DMSO|dBET")], metric="correlation", figsize=(4, 4),
    rasterized=True, yticklabels=False, cbar_kws={"label": "Concentration"},
    col_cluster=False, robust=True)
grid3.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.heatmap.sorted.dmso_dbet-specific_regions.svg"),
    dpi=300, bbox_inches="tight")


fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 4))
sns.barplot(
    data=p.loc[:, ~p.columns.str.contains("DMSO|dBET")].melt(),
    y="variable", x="value", orient="horiz",
    ax=axis[0])
sns.barplot(
    data=p.loc[:, p.columns.str.contains("DMSO|dBET")].melt(),
    y="variable", x="value", orient="horiz",
    ax=axis[1])
fig.savefig(os.path.join(
        output_dir,
        "diffBind.concentration.signficant.barplot.svg"),
    dpi=300, bbox_inches="tight")

# calculate p-value between distribution of regions
r = list()
for a in res7.columns:
    for b in res7.columns:
        r.append([a, b] + list(mannwhitneyu(res7[a], res7[b], alternative='two-sided')))

r = pd.DataFrame(r, columns=['a', 'b', 'stat', 'p']).pivot_table(index='a', columns='b', values='p')
r.to_csv(os.path.join(output_dir, "diffBind.signficance.csv"))
