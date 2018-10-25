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


# Get GREAT enrichments
for ip in ['BRD4', "MTHFD1"]:
    e = pd.read_csv("{}_WT.great.csv".format(ip))
    fig, axis = plt.subplots(1, figsize=(4, 4))
    axis.set_xlabel("Fold Enrichment")
    axis.set_ylabel("-log10(FDR)")
    for ontology, color in zip(e['Ontology'].unique(), sns.color_palette("colorblind")):
        axis.scatter(
            e.loc[e['Ontology'] == ontology, 'Hyper Fold Enrichment'],
            -np.log10(e.loc[e['Ontology'] == ontology, 'Hyper FDR Q-Val']), color=color, label=ontology)
    axis.legend()
    done = list()
    for t in e.sort_values("Hyper FDR Q-Val").head(10).index:
        if e.loc[t, 'Term Name'] not in done:
            axis.text(
                x=e.loc[t, 'Hyper Fold Enrichment'],
                y=-np.log10(e.loc[t, 'Hyper FDR Q-Val']), s=e.loc[t, 'Term Name'])
            done.append(e.loc[t, 'Term Name'])
    for t in e.sort_values("Hyper Fold Enrichment").tail(10).index:
        if e.loc[t, 'Term Name'] not in done:
            axis.text(
                x=e.loc[t, 'Hyper Fold Enrichment'],
                y=-np.log10(e.loc[t, 'Hyper FDR Q-Val']), s=e.loc[t, 'Term Name'])
            done.append(e.loc[t, 'Term Name'])
    sns.despine(fig)
    fig.savefig(os.path.join("{}_WT.great.scatter.svg".format(ip)), dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, figsize=(4, 4))
    e["Hyper FDR Q-Val log"] = -np.log10(e["Hyper FDR Q-Val"])
    axis.set_ylabel("Hyper FDR Q-Val log")
    axis.set_xlabel("GO terms")
    sns.barplot(data=e.sort_values('Hyper FDR Q-Val log', ascending=False), x='Hyper FDR Q-Val log', y='Term Name', hue="Ontology", orient="horiz")
    sns.despine(fig)
    fig.savefig(os.path.join("{}_WT.great.barplot.svg".format(ip)), dpi=300, bbox_inches="tight")
