#!/usr/bin/env python

"""
"""

if __name__ == '__main__':
    import matplotlib
    matplotlib.use("Agg")

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from ngs_toolkit.general import parse_great_enrichment
from scipy.stats import zscore

sns.set_style("white")
plt.rcParams['svg.fonttype'] = 'none'


output_dir = os.path.join("results", "unsupervised_analysis_ChIP-seq")


# Get GREAT enrichments
for ip in ['BRD4', "MTHFD1"]:
    e = parse_great_enrichment(os.path.join(output_dir, "{}.greatExportAll.tsv".format(ip)))
    e.to_csv(os.path.join(output_dir, "{}.greatExportAll.parsed.csv".format(ip)), index=False)

    # Plot
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.set_xlabel("Fold Enrichment")
    axis.set_ylabel("-log10(FDR)")
    for ontology, color in zip(e['Ontology'].unique(), sns.color_palette("colorblind")):
        axis.scatter(
            e.loc[e['Ontology'] == ontology, 'GeneFoldEnrich'],
            -np.log10(e.loc[e['Ontology'] == ontology, 'HyperFdrQ']), color=color, label=ontology)
    axis.legend(loc='best')
    e['extr'] = zscore(-np.log10(e["HyperFdrQ"])) * zscore(e["GeneFoldEnrich"])
    done = list()
    for t in e.sort_values("HyperFdrQ").head(20).index:
        if e.loc[t, 'Desc'] not in done:
            axis.text(
                x=e.loc[t, 'GeneFoldEnrich'],
                y=-np.log10(e.loc[t, 'HyperFdrQ']), s=e.loc[t, 'Desc'])
            done.append(e.loc[t, 'Desc'])
    for t in e.sort_values("extr").tail(20).index:
        if e.loc[t, 'Desc'] not in done:
            axis.text(
                x=e.loc[t, 'GeneFoldEnrich'],
                y=-np.log10(e.loc[t, 'HyperFdrQ']), s=e.loc[t, 'Desc'])
            done.append(e.loc[t, 'Desc'])
    for t in e.sort_values("GeneFoldEnrich").tail(20).index:
        if e.loc[t, 'Desc'] not in done:
            axis.text(
                x=e.loc[t, 'GeneFoldEnrich'],
                y=-np.log10(e.loc[t, 'HyperFdrQ']), s=e.loc[t, 'Desc'])
            done.append(e.loc[t, 'Desc'])
    sns.despine(fig)
    fig.savefig(os.path.join(
            output_dir,
            "{}_WT.great.scatter.svg".format(ip)),
        dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(4, 4))
    e["HyperFdrQ log"] = -np.log10(e["HyperFdrQ"])
    axis.set_ylabel("HyperFdrQ log")
    axis.set_xlabel("GO terms")
    sns.barplot(
        data=e.sort_values('HyperFdrQ log', ascending=False).head(50), dodge=False,
        x='HyperFdrQ log', y='Desc', hue="Ontology", orient="horiz")
    sns.despine(fig)
    fig.savefig(os.path.join(output_dir, "{}_WT.great.barplot.svg".format(ip)), dpi=300, bbox_inches="tight")
