#!/usr/bin/env python

"""
This is the main script of the mthfd1 project.
"""

import os
import sys
from looper.models import Project
import pybedtools
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import cm
import multiprocessing
import parmap
import pysam
import pandas as pd
import numpy as np
import cPickle as pickle
from collections import Counter

# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


def pickle_me(function):
    """
    Decorator for some methods of Analysis class.
    """
    def wrapper(obj, *args):
        function(obj, *args)
        pickle.dump(obj, open(obj.pickle_file, 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    return wrapper


class Analysis(object):
    """
    Class to hold functions and data from analysis.
    """

    def __init__(
            self,
            data_dir=os.path.join(".", "data"),
            results_dir=os.path.join(".", "results"),
            pickle_file=os.path.join(".", "data", "analysis.pickle"),
            samples=None,
            prj=None,
            from_pickle=False,
            **kwargs):
        # parse kwargs with default
        self.data_dir = data_dir
        self.results_dir = results_dir
        self.samples = samples
        self.pickle_file = pickle_file

        for directory in [self.data_dir, self.results_dir]:
            if not os.path.exists(directory):
                os.makedirs(directory)

        # parse remaining kwargs
        self.__dict__.update(kwargs)

        # reload itself if required
        if from_pickle:
            self.__dict__.update(self.from_pickle().__dict__)

    @pickle_me
    def to_pickle(self):
        pass

    def from_pickle(self):
        return pickle.load(open(self.pickle_file, 'rb'))

    @pickle_me
    def get_consensus_sites(self, samples):
        """Get consensus (union) sites across samples"""

        for i, sample in enumerate(samples):
            print(sample.name)
            # Get peaks
            peaks = pybedtools.BedTool(sample.peaks)
            # Merge overlaping peaks within a sample
            peaks = peaks.merge()
            if i == 0:
                sites = peaks
            else:
                # Concatenate all peaks
                sites = sites.cat(peaks)

        # Merge overlaping peaks across samples
        sites = sites.merge()

        # Filter
        # remove blacklist regions
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))
        # remove chrM peaks and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(os.path.join(self.results_dir, "mthfd1_peaks.bed"))

        # Read up again
        self.sites = pybedtools.BedTool(os.path.join(self.results_dir, "mthfd1_peaks.bed"))

    @pickle_me
    def get_consensus_sites_from_comparisons(
            self, comparisons, samples,
            min_width=[(["H3K27ac"], 1000), (["BRD4", "MTHFD1"], 500)]):
        """Get consensus (union) of sites across samples"""
        self.sites = dict()

        # For each comparison
        for j, comparison in enumerate(comparisons['comparison'].drop_duplicates()):
            # Raise on comparisons with only one side (incomplete still)
            if len(set(comparisons[(comparisons["comparison"] == comparison)]["comparison_side"].tolist())) != 2:
                raise ValueError("Comparison '%s' only has one side" % comparison)

            # Get peaks
            peaks = pybedtools.BedTool(comparisons['peaks'].drop_duplicates().ix[comparisons['comparison'] == comparison].tolist()[0])
            comparisons.loc[(comparisons["comparison"] == comparison), "number_peaks"] = len(peaks)
            print(j, comparison, len(peaks))

            # Merge overlaping peaks within a sample
            peaks = peaks.merge()

            # Extend peaks if under a certain size (depending on factor)
            for factors, width in min_width:
                if any([f in comparison for f in factors]):
                    print("Selected minimal width for comparison %s is %i bp." % (comparison, width))
                    break
            peaks = extend_sites_to_minimum(peaks, min_width=width)

            if j == 0:
                sites = peaks
            else:
                # Concatenate all peaks
                sites = sites.cat(peaks)

            # Merge overlaping peaks across samples
            sites = sites.merge()
            print(j, comparison, len(sites))

        # Filter out blacklisted regions
        output_name = os.path.join(self.data_dir, "mthfd1_peaks.bed")
        blacklist = pybedtools.BedTool(os.path.join(self.data_dir, "external", "wgEncodeDacMapabilityConsensusExcludable.bed"))
        # remove blacklist regions and save
        sites.intersect(v=True, b=blacklist).filter(lambda x: x.chrom != 'chrM').saveas(output_name)
        sites.saveas(output_name)

        # Read up again
        self.sites = pybedtools.BedTool(output_name)

    def calculate_peak_support(self, comparisons, samples, min_width=[(["H3K27ac"], 1000), (["BRD4", "MTHFD1"], 500)]):
        # calculate support (number of samples overlaping each merged peak)
        for i, comparison in enumerate(comparisons['comparison'].drop_duplicates()):
            print(comparison)
            peaks = pybedtools.BedTool(comparisons['peaks'].drop_duplicates().ix[comparisons['comparison'] == comparison].tolist()[0])

            # Extend peaks if under a certain size (depending on factor)
            for factors, width in min_width:
                if any([f in comparison for f in factors]):
                    print("Selected minimal width for comparison %s is %i bp." % (comparison, width))
                    break
            peaks = extend_sites_to_minimum(peaks, min_width=width)

            if i == 0:
                support = self.sites.intersect(peaks, wa=True, c=True)
            else:
                support = support.intersect(peaks, wa=True, c=True)
        support = support.to_dataframe()
        # support = support.reset_index()
        support.columns = ["chrom", "start", "end"] + comparisons['comparison'].drop_duplicates().tolist()
        support.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.binary_overlap_support.csv"), index=False)

        # get % of total consensus regions found per sample
        m = pd.melt(support, ["chrom", "start", "end"], var_name="sample_name")
        # groupby
        n = m.groupby("sample_name").apply(lambda x: len(x[x["value"] == 1]))

        # divide sum (of unique overlaps) by total to get support value between 0 and 1
        support["support"] = support[range(len(samples))].apply(lambda x: sum([i if i <= 1 else 1 for i in x]) / float(len(self.samples)), axis=1)
        # save
        support.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.support.csv"), index=False)

        self.support = support

    @pickle_me
    def measure_coverage(self, samples):
        # Count reads with pysam
        # make strings with intervals
        sites_str = [str(i.chrom) + ":" + str(i.start) + "-" + str(i.stop) for i in self.sites]
        # count, create dataframe
        self.coverage = pd.DataFrame(
            map(
                lambda x:
                    pd.Series(x),
                    parmap.map(
                        count_reads_in_intervals,
                        [sample.filtered for sample in samples],
                        sites_str,
                        parallel=True
                    )
            ),
            index=[sample.name for sample in samples]
        ).T

        # Add interval description to df
        ints = map(
            lambda x: (
                x.split(":")[0],
                x.split(":")[1].split("-")[0],
                x.split(":")[1].split("-")[1]
            ),
            self.coverage.index
        )
        self.coverage["chrom"] = [x[0] for x in ints]
        self.coverage["start"] = [int(x[1]) for x in ints]
        self.coverage["end"] = [int(x[2]) for x in ints]

        # save to disk
        self.coverage.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.raw_coverage.tsv"), sep="\t", index=True)

    @pickle_me
    def normalize_coverage_quantiles(self, samples):
        # Normalize by quantiles
        to_norm = self.coverage.loc[:, [s.name for s in samples]]
        self.coverage_qnorm = pd.DataFrame(
            normalize_quantiles_r(np.array(to_norm)),
            index=to_norm.index,
            columns=to_norm.columns
        )
        # Log2 transform
        self.coverage_qnorm = np.log2(1 + self.coverage_qnorm)

        self.coverage_qnorm = self.coverage_qnorm.join(self.coverage[['chrom', 'start', 'end']])
        self.coverage_qnorm.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.coverage_qnorm.log2.tsv"), sep="\t", index=True)

    def get_peak_gene_annotation(self):
        """
        Annotates peaks with closest gene.
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with hg19 TSS positions
        hg19_ensembl_tss = pybedtools.BedTool(os.path.join(self.data_dir, "external", "ensembl_tss.bed"))
        # get closest TSS of each cll peak
        closest = self.sites.closest(hg19_ensembl_tss, d=True).to_dataframe()[['chrom', 'start', 'end', 'thickStart', 'blockCount']]
        closest.columns = ['chrom', 'start', 'end', 'ensembl_transcript_id', 'distance']

        # add gene name and ensemble_gene_id
        ensembl_gtn = pd.read_table(os.path.join(self.data_dir, "external", "ensemblToGeneName.txt"), header=None)
        ensembl_gtn.columns = ['ensembl_transcript_id', 'gene_name']
        ensembl_gtp = pd.read_table(os.path.join(self.data_dir, "external", "ensGtp.txt"), header=None)[[0, 1]]
        ensembl_gtp.columns = ['ensembl_gene_id', 'ensembl_transcript_id']
        ensembl = pd.merge(ensembl_gtn, ensembl_gtp)

        gene_annotation = pd.merge(closest, ensembl, how="left")

        # aggregate annotation per peak, concatenate various genes (comma-separated)
        self.gene_annotation = gene_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.gene_annotation.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.gene_annotation.csv"), index=False)

        # save distances to all TSSs (for plotting)
        self.closest_tss_distances = closest['distance'].tolist()
        pickle.dump(self.closest_tss_distances, open(os.path.join(self.results_dir, "mthfd1_peaks.closest_tss_distances.pickle"), 'wb'))

    def get_peak_genomic_location(self):
        """
        Annotates peaks with its type of genomic location.
        Needs files downloaded by prepare_external_files.py
        """
        regions = [
            "ensembl_genes.bed", "ensembl_tss2kb.bed",
            "ensembl_utr5.bed", "ensembl_exons.bed", "ensembl_introns.bed", "ensembl_utr3.bed"]

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        for i, region in enumerate(regions):
            region_name = region.replace(".bed", "").replace("ensembl_", "")
            r = pybedtools.BedTool(os.path.join(self.data_dir, "external", region))
            if region_name == "genes":
                region_name = "intergenic"
                df = self.sites.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
                dfb = background.intersect(r, wa=True, f=0.2, v=True).to_dataframe()
            else:
                df = self.sites.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
                dfb = background.intersect(r, wa=True, u=True, f=0.2).to_dataframe()
            df['genomic_region'] = region_name
            dfb['genomic_region'] = region_name
            if i == 0:
                region_annotation = df
                region_annotation_b = dfb
            else:
                region_annotation = pd.concat([region_annotation, df])
                region_annotation_b = pd.concat([region_annotation_b, dfb])

        # sort
        region_annotation.sort(['chrom', 'start', 'end'], inplace=True)
        region_annotation_b.sort(['chrom', 'start', 'end'], inplace=True)
        # remove duplicates (there shouldn't be anyway)
        region_annotation = region_annotation.reset_index(drop=True).drop_duplicates()
        region_annotation_b = region_annotation_b.reset_index(drop=True).drop_duplicates()
        # join various regions per peak
        self.region_annotation = region_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()
        self.region_annotation_b = region_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(set([str(i) for i in x]))).reset_index()

        # save to disk
        self.region_annotation.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.region_annotation.csv"), index=False)
        self.region_annotation_b.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.region_annotation_background.csv"), index=False)

    def get_peak_chromatin_state(self):
        """
        Annotates peaks with chromatin states.
        (For now states are from CD19+ cells).
        Needs files downloaded by prepare_external_files.py
        """
        # create bedtool with CD19 chromatin states
        states_cd19 = pybedtools.BedTool(os.path.join(self.data_dir, "external", "HAP1_12_segments.annotated.bed"))

        # create background
        # shuffle regions in genome to create background (keep them in the same chromossome)
        background = self.sites.shuffle(genome='hg19', chrom=True)

        # intersect with cll peaks, to create annotation, get original peaks
        chrom_state_annotation = self.sites.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]
        chrom_state_annotation_b = background.intersect(states_cd19, wa=True, wb=True, f=0.2).to_dataframe()[['chrom', 'start', 'end', 'thickStart']]

        # aggregate annotation per peak, concatenate various annotations (comma-separated)
        self.chrom_state_annotation = chrom_state_annotation.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation.columns = ['chrom', 'start', 'end', 'chromatin_state']

        self.chrom_state_annotation_b = chrom_state_annotation_b.groupby(['chrom', 'start', 'end']).aggregate(lambda x: ",".join(x)).reset_index()
        self.chrom_state_annotation_b.columns = ['chrom', 'start', 'end', 'chromatin_state']

        # save to disk
        self.chrom_state_annotation.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.chromatin_state.csv"), index=False)
        self.chrom_state_annotation_b.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.chromatin_state_background.csv"), index=False)

    @pickle_me
    def annotate(self, samples):
        # add closest gene
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm,
            self.gene_annotation, on=['chrom', 'start', 'end'], how="left")
        # add genomic location
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.region_annotation[['chrom', 'start', 'end', 'genomic_region']], on=['chrom', 'start', 'end'], how="left")
        # add chromatin state
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.chrom_state_annotation[['chrom', 'start', 'end', 'chromatin_state']], on=['chrom', 'start', 'end'], how="left")

        # add support
        self.coverage_qnorm_annotated = pd.merge(
            self.coverage_qnorm_annotated,
            self.support[['chrom', 'start', 'end', 'support']], on=['chrom', 'start', 'end'], how="left")

        # calculate mean coverage
        self.coverage_qnorm_annotated['mean'] = self.coverage_qnorm_annotated[[sample.name for sample in samples]].mean(axis=1)
        # calculate coverage variance
        self.coverage_qnorm_annotated['variance'] = self.coverage_qnorm_annotated[[sample.name for sample in samples]].var(axis=1)
        # calculate std deviation (sqrt(variance))
        self.coverage_qnorm_annotated['std_deviation'] = np.sqrt(self.coverage_qnorm_annotated['variance'])
        # calculate dispersion (variance / mean)
        self.coverage_qnorm_annotated['dispersion'] = self.coverage_qnorm_annotated['variance'] / self.coverage_qnorm_annotated['mean']
        # calculate qv2 (std / mean) ** 2
        self.coverage_qnorm_annotated['qv2'] = (self.coverage_qnorm_annotated['std_deviation'] / self.coverage_qnorm_annotated['mean']) ** 2

        # calculate "amplitude" (max - min)
        self.coverage_qnorm_annotated['amplitude'] = (
            self.coverage_qnorm_annotated[[sample.name for sample in samples]].max(axis=1) -
            self.coverage_qnorm_annotated[[sample.name for sample in samples]].min(axis=1)
        )

        # Pair indexes
        assert self.coverage.shape[0] == self.coverage_qnorm_annotated.shape[0]
        self.coverage_qnorm_annotated.index = self.coverage.index

        # Save
        self.coverage_qnorm_annotated.to_csv(os.path.join(self.results_dir, "mthfd1_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index=True)

    def plot_peak_characteristics(self):
        # Loop at summary statistics:
        # interval lengths
        fig, axis = plt.subplots()
        sns.distplot([interval.length for interval in self.sites], bins=300, kde=False, ax=axis)
        axis.set_xlim(0, 2000)  # cut it at 2kb
        axis.set_xlabel("peak width (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.lengths.svg"), bbox_inches="tight")
        plt.close("all")

        # plot support
        fig, axis = plt.subplots()
        sns.distplot(self.support["support"], bins=40, ax=axis, kde=False)
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.support.svg"), bbox_inches="tight")
        plt.close("all")

        # Plot distance to nearest TSS
        fig, axis = plt.subplots()
        sns.distplot(self.closest_tss_distances, bins=200, ax=axis, kde=False)
        axis.set_xlim(0, 100000)  # cut it at 100kb
        axis.set_xlabel("distance to nearest TSS (bp)")
        axis.set_ylabel("frequency")
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.tss_distance.svg"), bbox_inches="tight")
        plt.close("all")

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = [item for sublist in self.region_annotation['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]
        all_region_annotation_b = [item for sublist in self.region_annotation_b['genomic_region'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_region_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_region_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        # plot individually
        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((data[1] / background[1]).astype(float))]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("genomic region")
        axis[2].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.genomic_regions.svg"), bbox_inches="tight")
        plt.close("all")

        # # Plot chromatin states
        # get long list of chromatin states (for plotting)
        all_chrom_state_annotation = [item for sublist in self.chrom_state_annotation['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]
        all_chrom_state_annotation_b = [item for sublist in self.chrom_state_annotation_b['chromatin_state'].apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_chrom_state_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort([1], ascending=False)
        # also for background
        background = Counter(all_chrom_state_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((data[1] / background[1]).astype(float))]).T, ax=axis[2])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("chromatin state")
        axis[2].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.chromatin_states.svg"), bbox_inches="tight")

        # distribution of count attributes
        data = self.coverage_qnorm_annotated.copy()

        fig, axis = plt.subplots(1)
        sns.distplot(data["mean"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.mean.distplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["qv2"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.qv2.distplot.svg"), bbox_inches="tight")

        fig, axis = plt.subplots(1)
        sns.distplot(data["dispersion"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.dispersion.distplot.svg"), bbox_inches="tight")

        # this is loaded now
        df = pd.read_csv(os.path.join(self.data_dir, "mthfd1_peaks.support.csv"))
        fig, axis = plt.subplots(1)
        sns.distplot(df["support"], rug=False, ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join(self.results_dir, "mthfd1.support.distplot.svg"), bbox_inches="tight")

        plt.close("all")

    def plot_coverage(self):
        data = self.coverage_qnorm_annotated.copy()
        # (rewrite to avoid putting them there in the first place)
        variables = ['gene_name', 'genomic_region', 'chromatin_state']

        for variable in variables:
            d = data[variable].str.split(',').apply(pd.Series).stack()  # separate comma-delimited fields
            d.index = d.index.droplevel(1)  # returned a multiindex Series, so get rid of second index level (first is from original row)
            data = data.drop([variable], axis=1)  # drop original column so there are no conflicts
            d.name = variable
            data = data.join(d)  # joins on index

        variables = [
            'chrom', 'start', 'end',
            'ensembl_transcript_id', 'distance', 'ensembl_gene_id', 'support',
            'mean', 'variance', 'std_deviation', 'dispersion', 'qv2',
            'amplitude', 'gene_name', 'genomic_region', 'chromatin_state']
        # Plot
        data_melted = pd.melt(
            data,
            id_vars=variables, var_name="sample", value_name="norm_counts")

        # transform dispersion
        data_melted['dispersion'] = np.log2(1 + data_melted['dispersion'])

        # Together in same violin plot
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        # dispersion
        fig, axis = plt.subplots(1)
        sns.violinplot("genomic_region", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.per_genomic_region.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "norm_counts", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "dispersion", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        fig, axis = plt.subplots(1)
        sns.violinplot("chromatin_state", "qv2", data=data_melted, ax=axis)
        fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.chromatin_state.violinplot.png"), bbox_inches="tight", dpi=300)

        # separated by variable in one grid
        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.mean.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="genomic_region", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.support.per_genomic_region.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "mean", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.mean.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "dispersion", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.dispersion.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "qv2", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.qv2.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)

        g = sns.FacetGrid(data_melted, col="chromatin_state", col_wrap=3)
        g.map(sns.distplot, "support", hist=False, rug=False)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts.support.chromatin_state.distplot.png"), bbox_inches="tight", dpi=300)
        plt.close("all")

    def plot_variance(self, samples):

        g = sns.jointplot('mean', "dispersion", data=self.coverage_qnorm_annotated, kind="kde")
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.dispersion.png"), bbox_inches="tight", dpi=300)

        g = sns.jointplot('mean', "qv2", data=self.coverage_qnorm_annotated)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.qv2_vs_mean.png"), bbox_inches="tight", dpi=300)

        g = sns.jointplot('support', "qv2", data=self.coverage_qnorm_annotated)
        g.fig.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.support_vs_qv2.png"), bbox_inches="tight", dpi=300)

        # Filter out regions which the maximum across all samples is below a treshold
        filtered = self.coverage_qnorm_annotated[self.coverage_qnorm_annotated[[sample.name for sample in samples]].max(axis=1) > 3]

        sns.jointplot('mean', "dispersion", data=filtered)
        plt.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.dispersion.filtered.png"), bbox_inches="tight", dpi=300)
        plt.close('all')
        sns.jointplot('mean', "qv2", data=filtered)
        plt.savefig(os.path.join(self.results_dir, "norm_counts_per_sample.support_vs_qv2.filtered.png"), bbox_inches="tight", dpi=300)


def add_args(parser):
    """
    Options for project and pipelines.
    """
    # Behaviour
    parser.add_argument("-g", "--generate", dest="stats", action="store_true",
                        help="Should we generate data and plots? Default=False")

    return parser


def count_reads_in_intervals(bam, intervals):
    """
    Counts total number of reads in a iterable holding strings
    representing genomic intervals of the type chrom:start-end.
    """
    counts = dict()

    bam = pysam.AlignmentFile(bam, mode='rb')

    chroms = ["chr" + str(x) for x in range(1, 23)] + ["chrX", "chrY"]

    for interval in intervals:
        if interval.split(":")[0] not in chroms:
            continue
        counts[interval] = bam.count(region=interval.split("|")[0])
    bam.close()

    return counts


def extend_sites_to_minimum(sites, genome="hg19", min_width=1000):
    """Extend sites smaller than `min_width`"""

    # get center
    tmpsites = sites.to_dataframe()
    tmpsites["width"] = tmpsites.apply(lambda x: int(x['end']) - int(x['start']), axis=1)
    tmpsites["new_start"] = tmpsites['start'] + (tmpsites['width'] / 2.).astype(np.int64)
    tmpsites["new_end"] = tmpsites["new_start"] + 1

    # intermediate save of smaller
    tmpsites[tmpsites["width"] < min_width][['chrom', 'new_start', 'new_end']].to_csv("tmp.bed", sep="\t", header=None, index=False)

    # extend
    extended_sites = pybedtools.BedTool("tmp.bed").slop(b=min_width / 2, genome=genome).to_dataframe()
    extended_sites = extended_sites.append(tmpsites[tmpsites["width"] >= min_width][['chrom', 'start', 'end']])
    extended_sites.to_csv("tmp.bed", sep="\t", header=None, index=False)

    # overwrite assignment and bed file
    return pybedtools.BedTool("tmp.bed").sort()


def normalize_quantiles_r(array):
    # install package
    # R
    # source('http://bioconductor.org/biocLite.R')
    # biocLite('preprocessCore')

    import rpy2.robjects as robjects
    import rpy2.robjects.numpy2ri
    rpy2.robjects.numpy2ri.activate()

    robjects.r('require("preprocessCore")')
    normq = robjects.r('normalize.quantiles')
    return np.array(normq(array))


def pca_r(x, colors, output_pdf):
    import rpy2.robjects as robj
    import pandas.rpy.common as com

    # save csvs for pca
    pd.DataFrame(x).T.to_csv('pca_file.csv', index=False)
    pd.Series(colors).to_csv('colors.csv', index=False)

    result = com.convert_robj(robj.r("""
    df = read.csv('pca_file.csv')

    colors = read.csv('colors.csv', header=FALSE)

    df.pca <- prcomp(df,
                     center = TRUE,
                     scale. = TRUE)
    return(df.pca)
    """))
    x = result['x']
    variance = result['sdev']

    # plot PC1 vs PC2
    fig, axis = plt.subplots(nrows=1, ncols=2)
    fig.set_figheight(10)
    fig.set_figwidth(25)

    # 1vs2 components
    for i in range(1, x.shape[0] + 1):
        axis[0].scatter(
            x.loc[i, 'PC1'], x.loc[i, 'PC2'],
            color=colors[i - 1],
            s=50
        )
    axis[0].set_xlabel("PC1 - {0}% variance".format(variance[0]))
    axis[0].set_ylabel("PC2 - {0}% variance".format(variance[1]))

    # plot PC1 vs PC3
    for i in range(1, x.shape[0] + 1):
        axis[1].scatter(
            x.loc[i, 'PC1'], x.loc[i, 'PC3'],
            color=colors[i - 1],
            s=50
        )
    axis[1].set_xlabel("PC1 - {0}% variance".format(variance[0]))
    axis[1].set_ylabel("PC3 - {0}% variance".format(variance[2]))

    fig.savefig(output_pdf, bbox_inches='tight')


def lola(bed_files, universe_file, output_folder):
    """
    Performs location overlap analysis (LOLA) on bedfiles with regions sets.
    """
    import rpy2.robjects as robj

    run = robj.r("""
        function(bedFiles, universeFile, outputFolder) {
            library("LOLA")

            userUniverse  <- LOLA::readBed(universeFile)

            dbPath1 = "/data/groups/lab_bock/shared/resources/regions/LOLACore/hg19/"
            dbPath2 = "/data/groups/lab_bock/shared/resources/regions/customRegionDB/hg19/"
            regionDB = loadRegionDB(c(dbPath1, dbPath2))

            if (typeof(bedFiles) == "character") {
                userSet <- LOLA::readBed(bedFiles)
                lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                lolaResults[order(support, decreasing=TRUE), ]
                writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
            } else if (typeof(bedFiles) == "double") {
                for (bedFile in bedFiles) {
                    userSet <- LOLA::readBed(bedFile)
                    lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=12)
                    lolaResults[order(support, decreasing=TRUE), ]
                    writeCombinedEnrichment(lolaResults, outFolder=outputFolder)
                }
            }
        }
    """)

    # convert the pandas dataframe to an R dataframe
    run(bed_files, universe_file, output_folder)

    # for F in `find . -iname *_regions.bed`
    # do
    #     if  [ ! -f `dirname $F`/allEnrichments.txt ]; then
    #         echo $F
    #         sbatch -J LOLA_${F} -o ${F/_regions.bed/_lola.log} ~/run_LOLA.sh $F ~/projects/mthfd1/results/mthfd1_peaks.bed hg19 `dirname $F`
    #     fi
    # done


def bed_to_fasta(bed_file, fasta_file):
    # write name column
    bed = pd.read_csv(bed_file, sep='\t', header=None)
    bed['name'] = bed[0] + ":" + bed[1].astype(str) + "-" + bed[2].astype(str)
    bed[1] = bed[1].astype(int)
    bed[2] = bed[2].astype(int)
    bed.to_csv(bed_file + ".tmp.bed", sep='\t', header=None, index=False)

    # do enrichment
    cmd = "twoBitToFa ~/resources/genomes/hg19/hg19.2bit -bed={0} {1}".format(bed_file + ".tmp.bed", fasta_file)

    os.system(cmd)
    # os.system("rm %s" % bed_file + ".tmp.bed")


def meme_ame(input_fasta, output_dir, background_fasta=None):
    # shuffle input in no background is provided
    if background_fasta is None:
        shuffled = input_fasta + ".shuffled"
        cmd = """
        fasta-dinucleotide-shuffle -c 1 -f {0} > {1}
        """.format(input_fasta, shuffled)
        os.system(cmd)

    cmd = """
    ame --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 0.05 \\
    --control {0} -o {1} {2} ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme
    """.format(background_fasta if background_fasta is not None else shuffled, output_dir, input_fasta)
    os.system(cmd)

    os.system("rm %s" % shuffled)

    # for F in `find . -iname *fa`
    # do
    #     if  [ ! -f `dirname $F`/ame.txt ]; then
    #         echo $F
    #         sbatch -J MEME-AME_${F} -o ${F/fa/ame.log} ~/run_AME.sh $F human
    #     fi
    # done


def parse_ame(ame_dir):

    with open(os.path.join(ame_dir, "ame.txt"), 'r') as handle:
        lines = handle.readlines()

    output = list()
    for line in lines:
        # skip header lines
        if line[0] not in [str(i) for i in range(10)]:
            continue

        # get motif string and the first half of it (simple name)
        motif = line.strip().split(" ")[5].split("_")[0]
        # get corrected p-value
        q_value = float(line.strip().split(" ")[-2])
        # append
        output.append((motif, q_value))

    return pd.Series(dict(output))


def enrichr(dataframe, gene_set_libraries=None, kind="genes"):
    """
    Use Enrichr on a list of genes (currently only genes supported through the API).
    """
    import json
    import requests

    ENRICHR_ADD = 'http://amp.pharm.mssm.edu/Enrichr/addList'
    ENRICHR_RETRIEVE = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
    query_string = '?userListId=%s&backgroundType=%s'

    if gene_set_libraries is None:
        gene_set_libraries = [
            'GO_Biological_Process_2015',
            "ChEA_2015",
            "KEGG_2016",
            "WikiPathways_2016",
            "Reactome_2016",
            "BioCarta_2016",
            "NCI-Nature_2016"
        ]

    results = pd.DataFrame()
    for gene_set_library in gene_set_libraries:
        print("Using enricher on %s gene set library." % gene_set_library)

        if kind == "genes":
            # Build payload with bed file
            attr = "\n".join(dataframe["gene_name"].dropna().tolist())
        elif kind == "regions":
            # Build payload with bed file
            attr = "\n".join(dataframe[['chrom', 'start', 'end']].apply(lambda x: "\t".join([str(i) for i in x]), axis=1).tolist())

        payload = {
            'list': (None, attr),
            'description': (None, gene_set_library)
        }
        # Request adding gene set
        response = requests.post(ENRICHR_ADD, files=payload)
        if not response.ok:
            raise Exception('Error analyzing gene list')

        # Track gene set ID
        user_list_id = json.loads(response.text)['userListId']

        # Request enriched sets in gene set
        response = requests.get(
            ENRICHR_RETRIEVE + query_string % (user_list_id, gene_set_library)
        )
        if not response.ok:
            raise Exception('Error fetching enrichment results')

        # Get enriched sets in gene set
        res = json.loads(response.text)
        # If there's no enrichemnt, continue
        if len(res) < 0:
            continue

        # Put in dataframe
        res = pd.DataFrame([pd.Series(s) for s in res[gene_set_library]])
        if res.shape[0] == 0:
            continue
        res.columns = ["rank", "description", "p_value", "z_score", "combined_score", "genes", "adjusted_p_value"]

        # Remember gene set library used
        res["gene_set_library"] = gene_set_library

        # Append to master dataframe
        results = results.append(res, ignore_index=True)

    return results

    # for F in `find /home/arendeiro/projects/mthfd1/results -name "*symbols.txt"`
    # do
    #     if  [ ! -f ${F/symbols.txt/enrichr.csv} ]; then
    #         echo $F
    #         sbatch -J ENRICHR_${F} -o ${F/symbols.txt/enrichr.log} ~/run_Enrichr.sh $F
    #     fi
    # done


def characterize_regions_structure(df, prefix, output_dir, universe_df=None):
    # use all sites as universe
    if universe_df is None:
        universe_df = pd.read_csv(os.path.join("results", "mthfd1_peaks.coverage_qnorm.log2.annotated.tsv"), sep="\t", index_col=0)

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # compare genomic regions and chromatin_states
    enrichments = pd.DataFrame()
    for i, var in enumerate(['genomic_region', 'chromatin_state']):
        # prepare:
        # separate comma-delimited fields:
        df_count = Counter(df[var].str.split(',').apply(pd.Series).stack().tolist())
        df_universe_count = Counter(universe_df[var].str.split(',').apply(pd.Series).stack().tolist())

        # divide by total:
        df_count = {k: v / float(len(df)) for k, v in df_count.items()}
        df_universe_count = {k: v / float(len(universe_df)) for k, v in df_universe_count.items()}

        # join data, sort by subset data
        both = pd.DataFrame([df_count, df_universe_count], index=['subset', 'all']).T
        both = both.sort("subset")
        both['region'] = both.index
        data = pd.melt(both, var_name="set", id_vars=['region']).replace(np.nan, 0)

        # sort for same order
        data.sort('region', inplace=True)

        # g = sns.FacetGrid(col="region", data=data, col_wrap=3, sharey=True)
        # g.map(sns.barplot, "set", "value")
        # plt.savefig(os.path.join(output_dir, "%s_regions.%s.svg" % (prefix, var)), bbox_inches="tight")

        fc = pd.DataFrame(np.log2(both['subset'] / both['all']), columns=['value'])
        fc['variable'] = var

        # append
        enrichments = enrichments.append(fc)

    # save
    enrichments.to_csv(os.path.join(output_dir, "%s_regions.region_enrichment.csv" % prefix), index=True)


def characterize_regions_function(df, output_dir, prefix, data_dir="data", universe_file=None):
    # use all sites as universe
    if universe_file is None:
        universe_file = os.path.join(data_dir, "mthfd1_peaks.bed")

    # make output dirs
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # save to bed
    bed_file = os.path.join(output_dir, "%s_regions.bed" % prefix)
    df[['chrom', 'start', 'end']].to_csv(bed_file, sep="\t", header=False, index=False)
    # save as tsv
    tsv_file = os.path.join(output_dir, "%s_regions.tsv" % prefix)
    df[['chrom', 'start', 'end']].reset_index().to_csv(tsv_file, sep="\t", header=False, index=False)

    # export ensembl gene names
    df['gene_name'].str.split(",").apply(pd.Series, 1).stack().drop_duplicates().to_csv(os.path.join(output_dir, "%s_genes.symbols.txt" % prefix), index=False)
    # export ensembl gene names
    df['ensembl_gene_id'].str.split(",").apply(pd.Series, 1).stack().drop_duplicates().to_csv(os.path.join(output_dir, "%s_genes.ensembl.txt" % prefix), index=False)

    # Motifs
    # de novo motif finding - enrichment
    fasta_file = os.path.join(output_dir, "%s_regions.fa" % prefix)
    bed_to_fasta(bed_file, fasta_file)

    meme_ame(fasta_file, output_dir)

    # Lola
    try:
        lola(bed_file, universe_file, output_dir)
    except:
        print("LOLA analysis for %s failed!" % prefix)

    # Enrichr
    results = enrichr(df[['chrom', 'start', 'end', "gene_name"]])

    # Save
    results.to_csv(os.path.join(output_dir, "%s_regions.enrichr.csv" % prefix), index=False, encoding='utf-8')


def plot_region_types_per_sample(analysis):
    # Plot region types per sample
    g_fig, g_axis = plt.subplots(3, 1)
    c_fig, c_axis = plt.subplots(3, 1)
    for i, sample in enumerate([s for s in analysis.samples if s.ip != "IgG"]):
        # Get regions specific to this sample
        # based on support
        index = analysis.support[analysis.support[sample.name] > 0].index

        # Plot genomic regions
        # these are just long lists with genomic regions
        all_region_annotation = [item for sublist in analysis.region_annotation.ix[index]['genomic_region'].dropna().apply(lambda x: x.split(",")) for item in sublist]
        all_region_annotation_b = [item for sublist in analysis.region_annotation_b['genomic_region'].dropna().apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_region_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort_values([0])
        # also for background
        background = Counter(all_region_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        # plot individually
        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((((data[1] / data[1].sum()) * 1000) / ((background[1] / background[1].sum()) * 1000)).astype(float))]).T, ax=axis[2])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((((data[1] / data[1].sum()) * 1000) / ((background[1] / background[1].sum()) * 1000)).astype(float))]).T, ax=g_axis[i])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("genomic region")
        axis[2].set_xlabel("genomic region")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(analysis.results_dir, "mthfd1.genomic_regions.{}.svg".format(sample.name)), bbox_inches="tight")
        plt.close("all")

        # # Plot chromatin states
        # get long list of chromatin states (for plotting)
        all_chrom_state_annotation = [item for sublist in analysis.chrom_state_annotation.ix[index]['chromatin_state'].dropna().apply(lambda x: x.split(",")) for item in sublist]
        all_chrom_state_annotation_b = [item for sublist in analysis.chrom_state_annotation_b['chromatin_state'].dropna().apply(lambda x: x.split(",")) for item in sublist]

        # count region frequency
        count = Counter(all_chrom_state_annotation)
        data = pd.DataFrame([count.keys(), count.values()]).T
        data = data.sort_values([0])
        # also for background
        background = Counter(all_chrom_state_annotation_b)
        background = pd.DataFrame([background.keys(), background.values()]).T
        background = background.ix[data.index]  # same sort order as in the real data

        fig, axis = plt.subplots(3, sharex=True, sharey=False)
        sns.barplot(x=0, y=1, data=data, ax=axis[0])
        sns.barplot(x=0, y=1, data=background, ax=axis[1])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((((data[1] / data[1].sum()) * 1000) / ((background[1] / background[1].sum()) * 1000)).astype(float))]).T, ax=axis[2])
        sns.barplot(x=0, y=1, data=pd.DataFrame([data[0], np.log2((((data[1] / data[1].sum()) * 1000) / ((background[1] / background[1].sum()) * 1000)).astype(float))]).T, ax=c_axis[i])
        axis[0].set_title("ATAC-seq peaks")
        axis[1].set_title("genome background")
        axis[2].set_title("peaks over background")
        axis[1].set_xlabel("chromatin state")
        axis[2].set_xlabel("chromatin state")
        axis[0].set_ylabel("frequency")
        axis[1].set_ylabel("frequency")
        axis[2].set_ylabel("fold-change")
        fig.autofmt_xdate()
        fig.tight_layout()
        sns.despine(fig)
        fig.savefig(os.path.join(analysis.results_dir, "mthfd1.chromatin_states.{}.svg".format(sample.name)), bbox_inches="tight")

        g_axis[i].set_ylabel(sample.name)
        c_axis[i].set_ylabel(sample.name)

    g_fig.savefig(os.path.join(analysis.results_dir, "mthfd1.genomic_regions.all.svg"), bbox_inches="tight")
    c_fig.savefig(os.path.join(analysis.results_dir, "mthfd1.chromatin_states.all.svg"), bbox_inches="tight")


def main():
    # Start project
    prj = Project("metadata/project_config.yaml")
    prj.add_sample_sheet()

    # temporary:
    for sample in prj.samples:
        sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
        sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
        sample.coverage = os.path.join(sample.paths.sample_root, "coverage", sample.name + ".cov")

    # Start analysis object
    to_include = [
        "HAP1_ChIPmentation_MTHFD1_C3_WT",
        "HAP1_ChIPmentation_BRD4_WT",
        "HAP1_ChIPmentation_H3K27ac_WT",
        "HAP1_ChIPmentation_IgG_WT",
    ]
    prj.samples = [sample for sample in prj.samples if sample.name in to_include]

    analysis = Analysis(prj=prj, samples=prj.samples)
    analysis.prj = prj

    # GET CONSENSUS PEAK SET, ANNOTATE IT, PLOT FEATURES
    # GET CHROMATIN OPENNESS MEASUREMENTS
    # PLOT STUFF
    comparisons = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
    comparisons = comparisons[comparisons['comparison'].isin(to_include[:-1])]
    comparisons['peaks'] = "data/peaks/" + comparisons['comparison'] + "_allctrls/" + comparisons['comparison'] + "_allctrls_peaks.narrowPeak"
    # Get consensus peak set from all samples
    analysis.get_consensus_sites_from_comparisons(comparisons, analysis.samples)
    # Calculate peak support
    analysis.calculate_peak_support(analysis, comparisons, analysis.samples, min_width=[(["H3K27ac"], 500), (["BRD4", "MTHFD1"], 500)])
    # Annotate peaks with closest gene
    analysis.get_peak_gene_annotation()
    # Annotate peaks with genomic regions
    analysis.get_peak_genomic_location()
    # Annotate peaks with chromatin states (HAP1)
    analysis.get_peak_chromatin_state()

    # WORK WITH "OPENNESS"
    # Get coverage values for each peak in each sample
    analysis.measure_coverage(analysis.samples)
    # normalize coverage values
    analysis.normalize_coverage_quantiles(analysis.samples)
    # Annotate peaks with closest gene, chromatin state,
    # genomic location, mean and variance measurements across samples
    analysis.annotate(analysis.samples)

    # Plots
    # plot general peak set features
    analysis.plot_peak_characteristics()
    # Plot rpkm features across peaks/samples
    analysis.plot_coverage()
    analysis.plot_variance(analysis.samples)

    # Characterize all regions as a whole

    # Overlap
    from scipy.stats import fisher_exact
    h3k4me3 = analysis.support[analysis.support['HAP1_ChIPmentation_H3K27ac_WT'] > 0]
    brd4 = analysis.support[analysis.support['HAP1_ChIPmentation_BRD4_WT'] > 0]
    mthfd1 = analysis.support[analysis.support['HAP1_ChIPmentation_MTHFD1_C3_WT'] > 0]
    table = [
        [(brd4['HAP1_ChIPmentation_MTHFD1_C3_WT'] > 0).sum(), (brd4['HAP1_ChIPmentation_MTHFD1_C3_WT'] == 0).sum()],  # BRD4: overlap, non-overlaping
        [(mthfd1['HAP1_ChIPmentation_BRD4_WT'] > 0).sum(), (mthfd1['HAP1_ChIPmentation_BRD4_WT'] == 0).sum()]   # MTHFD1: overlap, non-overlaping
    ]
    print(fisher_exact(table))

    # proximity
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    brd4 = pybedtools.BedTool(os.path.join(analysis.results_dir, "HAP1_ChIPmentation_BRD4_WT", "HAP1_ChIPmentation_BRD4_WT_regions.bed"))
    mthfd1 = pybedtools.BedTool(os.path.join(analysis.results_dir, "HAP1_ChIPmentation_MTHFD1_C3_WT", "HAP1_ChIPmentation_MTHFD1_C3_WT_regions.bed"))

    closest = mthfd1.closest(brd4, d=True)

    s = pd.Series(Counter(closest.to_dataframe()['thickStart']))

    dists = dict()
    for i in range(1000):
        dists[i] = s[s.index < i * 1000].sum() / float(s.sum())

    fig, axis = plt.subplots(1)
    axis.plot(dists.keys(), dists.values())
    axis.set_xlabel("Distance to nearest BRD4 peak (kb)")
    axis.set_ylabel("Fraction of total MTHFD1 peaks (n = {})".format(len(mthfd1)))
    ax = zoomed_inset_axes(axis, zoom=5, loc=4, axes_kwargs={"aspect": 100, "xlim": (0, 50), "ylim": (0, 0.5)})
    ax.plot(dists.keys(), dists.values())
    sns.despine(fig)
    fig.savefig(os.path.join(analysis.results_dir, "mthfd1_brd4.distances.svg"), bbox_inches="tight")

    # visuzalize gene-relative coverage
    with open("config.txt", "w") as handle:
        for sample in analysis.samples[:2] + [analysis.samples[-1]]:
            handle.write("\t".join([sample.filtered, "-1", sample.name + "\n"]))

    cmd = "ngs.plot.r -G hg19 -R genebody -L 3000 -C config.txt -O genebody_coverage -GO km"
    cmd = "ngs.plot.r -G hg19 -R genebody -C config.txt -O average_profile_coverage -D ensembl -FL 3000"

    # Only H3K27ac/BRD4/MTHFD1-bound (or proximity) peaks
    with open("config.txt", "w") as handle:
        for sample1 in prj.samples:
            for sample2 in prj.samples[:2] + [prj.samples[-1]]:
                handle.write(
                    "\t".join([
                        sample1.filtered,
                        os.path.join("results", sample2.name, sample2.name + "_genes.ensembl.txt"),
                        "{}_signal_on_{}_peaks".format(sample1.name, sample2.name) + "\n"]))

    cmd = "ngs.plot.r -G hg19 -R genebody -C config.txt -O average_profile_coverage_peakgenes -D ensembl -FL 3000"

    # Plot region types per sample
    plot_region_types_per_sample(analysis)

    # Get functional
    analysis.support.index = analysis.support["chrom"] + ":" + analysis.support['start'].astype(str) + "-" + analysis.support['end'].astype(str)
    analysis.coverage_qnorm_annotated['start'] = analysis.coverage_qnorm_annotated['start'].astype(int)
    analysis.coverage_qnorm_annotated['end'] = analysis.coverage_qnorm_annotated['end'].astype(int)

    for i, sample in enumerate([sample for sample in analysis.samples if sample.ip != "IgG"]):
        print(sample.name)
        index = analysis.support[analysis.support[sample.name] > 0].index

        characterize_regions_structure(
            df=analysis.coverage_qnorm_annotated.ix[index].dropna(),
            prefix=sample.name,
            output_dir=os.path.join(analysis.results_dir, sample.name))

        characterize_regions_function(
            df=analysis.coverage_qnorm_annotated.ix[index].dropna(),
            output_dir=os.path.join(analysis.results_dir, sample.name),
            prefix=sample.name)

    # De novo motif finding
    cmd = "meme-chip "
    cmd += "-meme-p 12 "
    cmd += "-oc {} ".format(os.path.join(analysis.results_dir, "HAP1_ChIPmentation_MTHFD1_C3_WT"))
    cmd += "-db ~/resources/motifs/motif_databases/HUMAN/HOCOMOCOv10.meme"
    cmd += "-index-name meme_chip.html "
    cmd += " {}".format(os.path.join(analysis.results_dir, "HAP1_ChIPmentation_MTHFD1_C3_WT", "HAP1_ChIPmentation_MTHFD1_C3_WT_regions.fa"))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
