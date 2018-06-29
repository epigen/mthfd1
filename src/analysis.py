#!/usr/bin/env python

"""
This is the main script of the mthfd1 project.
"""


if __name__ == '__main__':
    import matplotlib
    matplotlib.use('Agg')

import os
import sys
from looper.models import Project
import pybedtools
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
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


def atacseq_analysis():
    from ngs_toolkit.general import (collect_differential_enrichment,
                                differential_analysis,
                                differential_enrichment, differential_overlap,
                                plot_differential,
                                plot_differential_enrichment)

    # ATAC-seq
    # Start project and analysis objects
    for cell_line in ["HAP1", "A549"]:
        prj = Project("metadata/project_config.yaml")
        prj.samples = [sample for sample in prj.samples if sample.library == "ATAC-seq" and sample.cell_line == cell_line]
        for sample in prj.samples:
            sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
            sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
            sample.peaks = os.path.join(sample.paths.sample_root, "peaks", sample.name + "_peaks.narrowPeak")
        atac_analysis = ATACSeqAnalysis(name="mthfd1_atac_{}".format(cell_line), prj=prj, samples=prj.samples)

        atac_analysis.chromatin_state = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.chromatin_state.csv"))
        atac_analysis.chromatin_state_background = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.chromatin_state_background.csv"))
        atac_analysis.closest_tss_distances = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.closest_tss_distances.pickle"))
        atac_analysis.gene_annotation = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.gene_annotation.csv"))
        atac_analysis.region_annotation = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.region_annotation.csv"))
        atac_analysis.region_annotation_background = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.region_annotation_background.csv"))
        atac_analysis.support = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.support.csv"))
        atac_analysis.raw_coverage = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.raw_coverage.csv"))
        atac_analysis.coverage_qnorm = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.coverage_qnorm.csv"), index_col=0)
        atac_analysis.coverage_qnorm_annotated = pd.read_csv(os.path.join("results", "mthfd1_atac_HAP1_peaks.coverage_qnorm.annotated.csv"), index_col=0)
        atac_analysis.accessibility = pd.read_csv(os.path.join("results", "mthfd1_atac.accessibility.annotated_metadata.csv"), index_col=0, header=range(5))

        # Project's attributes
        sample_attributes = ["sample_name", "cell_line", "library", "condition", "replicate"]
        attributes_to_plot = ["cell_line", "condition", "replicate"]


        # Get consensus peak set from all samples
        atac_analysis.get_consensus_sites(atac_analysis.samples)
        atac_analysis.calculate_peak_support(atac_analysis.samples)
        atac_analysis.get_peak_gene_annotation()
        atac_analysis.get_peak_genomic_location()
        atac_analysis.get_peak_chromatin_state(os.path.join(atac_analysis.data_dir, "external", "HAP1_12_segments.annotated.bed"))

        # Get coverage values for each peak in each sample of ATAC-seq
        atac_analysis.measure_coverage(atac_analysis.samples)

        # Normalize ointly (quantile normalization + GC correction)
        atac_analysis.normalize(method="quantile")

        # Annotate peaks
        atac_analysis.annotate(quant_matrix="coverage_qnorm")

        # Annotate samples
        atac_analysis.accessibility = atac_analysis.annotate_with_sample_metadata(
            quant_matrix="coverage_annotated",
            attributes=sample_attributes)

        # Unsupervised analysis
        atac_analysis.unsupervised(
            quant_matrix="accessibility", samples=None,
            attributes_to_plot=attributes_to_plot, plot_prefix="accessibility_{}".format(cell_line))

        # Supervised analysis
        comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))
        comparison_table = comparison_table[comparison_table["toggle"] == 1]
        comparison_table = comparison_table[comparison_table["data_type"] == "ATAC-seq"]

        c = comparison_table[comparison_table["comparison_name"].str.contains(cell_line)]
        _ = differential_analysis(atac_analysis, c, data_type="ATAC-seq", overwrite=False, output_prefix="differential_analysis_{}".format(cell_line))

        # Visualize regions and samples found in the differential comparisons
        atacseq_results = pd.read_csv(
            os.path.join(
                "results",
                "differential_analysis_ATAC-seq",
                "differential_analysis_{}".format(cell_line) + ".deseq_result.all_comparisons.csv"), index_col=0)

        plot_differential(
            atac_analysis,
            atacseq_results,
            c,
            samples=None,
            data_type="ATAC-seq",
            alpha=0.05,
            corrected_p_value=True,
            fold_change=None,
            output_dir="results/differential_analysis_{data_type}",
            output_prefix="differential_analysis_{}".format(cell_line))

    atacseq_results = pd.DataFrame()
    for comp in comparison_table["comparison_name"].drop_duplicates():
        atacseq_results = atacseq_results.append(pd.read_csv(
            os.path.join("results", "differential_analysis_ATAC-seq", "differential_analysis.deseq_result.{}.csv".format(comp)),
            index_col=0).reset_index(), ignore_index=True)
    atacseq_results = atacseq_results.set_index("index")
    atacseq_results.to_csv(os.path.join("results", "differential_analysis_ATAC-seq", "differential_analysis" + ".deseq_result.all_comparisons.csv"), index=True)


def get_gene_expression_data():
    from pypiper import NGSTk
    import textwrap
    tk = NGSTk()

    gene_exp_diff_file = "/scratch/lab_bsf/projects/BSA_0197_Compound_shRNA_/hg38/rnaseq_cuffdiff_global/gene_exp.diff"
    header = ["test_id", "gene_id", "gene", "locus", "sample_1", "sample_2", "status", "value_1", "value_2", "log2(fold_change)", "test_stat", "p_value", "q_value", "significant"]
    output_dir = os.path.join("results", "gene_expression_cuffdiff")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    cell_lines = ['HAP1', 'A549', 'K562']
    comparisons = [
        ('DMSO_1h', 'JQ1_1h'),
        ('DMSO_6h', 'JQ1_6h'),
        ('DMSO_1h', 'JQMT_1h'),
        ('DMSO_6h', 'JQMT_6h'),
        ('DMSO_1h', 'MTX_1h'),
        ('DMSO_6h', 'MTX_6h'),
        ('DMSO_1h', 'dBet6_1h'),
        ('DMSO_6h', 'dBet6_6h'),
        ('shBRD4', 'shPLKO'),
        ('shMTHBRD', 'shPLKO'),
        ('shMTHFD1', 'shPLKO'),
    ]
    for cell_line in cell_lines:
        for condition, control in comparisons:
            a = "{}_{}".format(cell_line, control)
            b = "{}_{}".format(cell_line, condition)
            output_file = os.path.join(output_dir, "diff_expression.{}_vs_{}.csv".format(a, b))

            cmd = tk.slurm_header(output_file, output_file + ".log", queue='shortq', n_tasks=1, time='4:00:00', cpus_per_task=1, mem_per_cpu=4000)
            cmd += """grep -P "{}\\t{}" {} > {}\n""".format(a, b, gene_exp_diff_file, output_file)
            output_file2 = os.path.join(output_dir, "diff_expression.{}_vs_{}.csv".format(b, a))
            cmd += """\t\tgrep -P "{}\\t{}" {} > {}\n""".format(b, a, gene_exp_diff_file, output_file2)
            cmd += tk.slurm_footer()

            with open(output_file + ".sh", "w") as handle:
                handle.writelines(textwrap.dedent(cmd))
            tk.slurm_submit_job(output_file + ".sh")

    # Colllect, label
    df = pd.DataFrame()
    for cell_line in cell_lines:
        for condition, control in comparisons:
            a = "{}_{}".format(cell_line, control)
            b = "{}_{}".format(cell_line, condition)
            output_file = os.path.join(output_dir, "diff_expression.{}_vs_{}.csv".format(b, a))

            df2 = pd.read_csv(output_file, header=None, sep="\t", names=header)
            df = df.append(df2[['gene', 'sample_1', 'sample_2', 'value_1', 'value_2', 'q_value', 'log2(fold_change)']], ignore_index=True)
    df["comparison_name"] = (
        pd.Series(map(lambda x: x[0], df['sample_1'].str.split("_"))) + "_" +
        pd.Series(map(lambda x: x[1], df['sample_1'].str.split("_"))) + "_" +
        pd.Series(map(lambda x: x[1], df['sample_2'].str.split("_"))))
    df.to_csv(os.path.join(output_dir, "gene_exp.diff.relevant.csv"), index=False)

    df['mean'] = np.log2(1 + df[['value_1', 'value_2']].mean(1))
    df[df['q_value'] < 0.05].to_csv(os.path.join(output_dir, "gene_exp.diff.relevant.differential.csv"), index=False)

    return df


def main():
    from ngs_toolkit.atacseq import ATACSeqAnalysis
    from ngs_toolkit.chipseq import ChIPSeqAnalysis
    from ngs_toolkit.chipseq import homer_peaks_to_bed
    from ngs_toolkit.general import bed_to_fasta, homer_motifs, meme_ame

    # Start project and analysis objects
    prj = Project(os.path.join("metadata", "project_config.yaml"))
    for sample in prj.samples:
        sample.mapped = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.bam")
        sample.filtered = os.path.join(sample.paths.sample_root, "mapped", sample.name + ".trimmed.bowtie2.filtered.bam")
    analysis = ChIPSeqAnalysis(name="mthfd1", prj=prj, samples=prj.samples)

    # read in comparison table
    comparison_table = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))

    analysis.call_peaks_from_comparisons(
        comparison_table=comparison_table[comparison_table['comparison_type'] == 'peaks'], overwrite=False)

    # analysis.get_consensus_sites(
    #     comparison_table=comparison_table[comparison_table['comparison_type'] == 'peaks'],
    #     region_type="peaks", blacklist_bed="wgEncodeDacMapabilityConsensusExcludable.bed")

    # analysis.calculate_peak_support(
    #     comparison_table=comparison_table[comparison_table['comparison_type'] == 'peaks'])

    # # # Select only peaks called in all H3K27ac comparisons
    # # sup = analysis.support.loc[:, analysis.support.columns.str.contains("H3K27ac")].sum(axis=1)
    # # analysis.support.loc[sup == analysis.support.columns.str.contains("H3K27ac").sum(), :].index

    diffbind_annotation = pd.read_csv("metadata/diffBind_design.csv")
    for _, row in diffbind_annotation.iterrows():
        signal_sample = [s for s in analysis.samples if s.name == row["SampleID"]][0]
        control_sample = [s for s in analysis.samples if s.name == row["ControlID"]][0]
        # homer_call_chipseq_peak_job(
        #     signal_samples=[signal_sample], control_samples=[control_sample],
        #     output_dir=os.path.join(signal_sample.paths.sample_root, "peaks"), name=row["SampleID"])

        file = os.path.join(signal_sample.paths.sample_root, "peaks", row["SampleID"], row["SampleID"] + "_homer_peaks.")
        homer_peaks_to_bed(file + "narrowPeak", file + "bed")

    # # Get coverage
    # analysis.measure_coverage()

    # cov = analysis.coverage.loc[:, analysis.coverage.columns.str.contains("ChIP-seq")]
    # cov_tpm = np.log2(0.1 + (cov / cov.sum(axis=0)) * 1e6)

    # from ngs_toolkit.general import normalize_quantiles_r
    # coverage_qnorm = pd.DataFrame(
    #     normalize_quantiles_r(cov_tpm.values),
    #     index=cov_tpm.index,
    #     columns=cov_tpm.columns
    # )

    # # BRD4 vs MTHFD1
    # samples <- read.csv("metadata/diffBind_design.csv")
    # samples["Peaks"] <- samples["ReplicatePeaks"]
    # comparison <- "MTHFD1_vs_BRD4"
    # diffbind_factors <- dba(sampleSheet=samples[samples["Condition"] == "DMSO", ])
    # diffbind_factors <- dba.count(diffbind_factors, bCorPlot=FALSE)
    # diffbind_factors <- dba.contrast(diffbind_factors, categories=DBA_FACTOR, block=DBA_REPLICATE, minMembers=2)
    # diffbind_factors <- dba.analyze(diffbind_factors, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_factors.DB <- dba.report(diffbind_factors, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_factors.DB), paste("results/diffbind_analysis.replicate_peaks", comparison, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # samples <- read.csv("metadata/diffBind_design.csv")
    # samples["Peaks"] <- samples["H3K27acPeaks"]
    # comparison <- "MTHFD1_vs_BRD4"
    # diffbind_factors <- dba(sampleSheet=samples[samples["Condition"] == "DMSO", ])
    # diffbind_factors <- dba.count(diffbind_factors, bCorPlot=FALSE)
    # diffbind_factors <- dba.contrast(diffbind_factors, categories=DBA_FACTOR, block=DBA_REPLICATE, minMembers=2)
    # diffbind_factors <- dba.analyze(diffbind_factors, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_factors.DB <- dba.report(diffbind_factors, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_factors.DB), paste("results/diffbind_analysis.H3K27ac_peaks", comparison, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # # dBET treatment
    # library("DiffBind")
    # samples <- read.csv("metadata/diffBind_design.csv")
    # samples["Peaks"] <- samples["ReplicatePeaks"]
    # ip <- "BRD4"
    # diffbind_brd4 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    # diffbind_brd4 <- dba.count(diffbind_brd4, bCorPlot=FALSE)
    # diffbind_brd4 <- dba.contrast(diffbind_brd4, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    # diffbind_brd4 <- dba.analyze(diffbind_brd4, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_brd4.DB <- dba.report(diffbind_brd4, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_brd4.DB), paste("results/diffbind_analysis.replicate_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # ip <- "MTHFD1"
    # diffbind_mthfd1 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    # diffbind_mthfd1 <- dba.count(diffbind_mthfd1, bCorPlot=FALSE)
    # diffbind_mthfd1 <- dba.contrast(diffbind_mthfd1, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    # diffbind_mthfd1 <- dba.analyze(diffbind_mthfd1, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_mthfd1.DB <- dba.report(diffbind_mthfd1, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_mthfd1.DB), paste("results/diffbind_analysis.replicate_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # samples <- read.csv("metadata/diffBind_design.csv")
    # samples["Peaks"] <- samples["ComparisonPeaks"]
    # ip <- "BRD4"
    # diffbind_brd4 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    # diffbind_brd4 <- dba.count(diffbind_brd4, bCorPlot=FALSE)
    # diffbind_brd4 <- dba.contrast(diffbind_brd4, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    # diffbind_brd4 <- dba.analyze(diffbind_brd4, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_brd4.DB <- dba.report(diffbind_brd4, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_brd4.DB), paste("results/diffbind_analysis.comparison_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # ip <- "MTHFD1"
    # diffbind_mthfd1 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    # diffbind_mthfd1 <- dba.count(diffbind_mthfd1, bCorPlot=FALSE)
    # diffbind_mthfd1 <- dba.contrast(diffbind_mthfd1, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    # diffbind_mthfd1 <- dba.analyze(diffbind_mthfd1, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_mthfd1.DB <- dba.report(diffbind_mthfd1, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_mthfd1.DB), paste("results/diffbind_analysis.comparison_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)

    # dBET treatment
    samples <- read.csv("metadata/diffBind_design.csv")
    samples["Peaks"] <- samples["H3K27acPeaks"]
    ip <- "BRD4"
    diffbind_brd4 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    diffbind_brd4 <- dba.count(diffbind_brd4, bCorPlot=FALSE)
    diffbind_brd4 <- dba.contrast(diffbind_brd4, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    diffbind_brd4 <- dba.analyze(diffbind_brd4, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    diffbind_brd4.DB <- dba.report(diffbind_brd4, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    write.table(as.data.frame(diffbind_brd4.DB), paste("results/diffbind_analysis.H3K27ac_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)

    ip <- "MTHFD1"
    diffbind_mthfd1 <- dba(sampleSheet=samples[samples["Factor"] == ip, ])
    diffbind_mthfd1 <- dba.count(diffbind_mthfd1, bCorPlot=FALSE)
    diffbind_mthfd1 <- dba.contrast(diffbind_mthfd1, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    diffbind_mthfd1 <- dba.analyze(diffbind_mthfd1, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    diffbind_mthfd1.DB <- dba.report(diffbind_mthfd1, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    write.table(as.data.frame(diffbind_mthfd1.DB), paste("results/diffbind_analysis.H3K27ac_peaks", ip, "differential.csv", sep="."), sep=",", row.names=FALSE)


    # Curate peaks and get top N
    n = 500

    for factor in ["BRD4", "MTHFD1"]:
        output_bed = os.path.join("results", "diffbind_analysis.H3K27ac_peaks.{}.differential.top{}.bed".format(factor, 100))
        diff = pd.read_csv(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.{}.differential.csv".format(factor)))
        diff['name'] = diff['seqnames'] + ":" + diff['start'].astype(str) + "-" + diff['end'].astype(str)
        diff['score'] = -np.log10(diff['p.value'])
        diff.sort_values('score', ascending=False)[['seqnames', 'start', 'end', 'name', 'score']].head(n).to_csv(output_bed, sep="\t", header=None, index=False)
        cmd = "bedtools intersect -v -a {} -b {} > {}".format(
            output_bed, os.path.join("data", "external", "wgEncodeDacMapabilityConsensusExcludable.bed"), output_bed + 'no_blacklist.bed')
        os.system(cmd)

    # Plot pvalue of BRD4 changing regions vs MTHFD1
    brd4 = pd.read_csv(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.{}.differential.csv".format("BRD4")))
    mthfd1 = pd.read_csv(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.{}.differential.csv".format("MTHFD1")))
    brd4['-log_pvalue'] = -np.log10(brd4['p.value'])
    mthfd1['-log_pvalue'] = -np.log10(mthfd1['p.value'])

    fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3))
    axis[0, 1].scatter(brd4['-log_pvalue'], mthfd1['-log_pvalue'], alpha=0.5, s=2, rasterized=True)
    lim_max = max(brd4['-log_pvalue'].max(), mthfd1['-log_pvalue'].max())
    axis[0, 1].plot((0, lim_max), (0, lim_max), linestyle="--", color="black", alpha=0.5)
    axis[0, 1].set_xlim((0, lim_max))
    axis[0, 1].set_ylim((0, lim_max))
    axis[0, 1].set_xlabel("BRD4")
    axis[0, 1].set_ylabel("MTHFD1")

    sns.distplot(brd4['-log_pvalue'], ax=axis[0, 0], kde=False)
    lim_max = max(brd4['-log_pvalue'].max(), mthfd1['-log_pvalue'].max())
    axis[0, 0].set_xlim((0, lim_max))
    sns.distplot(mthfd1['-log_pvalue'], ax=axis[1, 1], kde=False)
    lim_max = max(mthfd1['-log_pvalue'].max(), mthfd1['-log_pvalue'].max())
    axis[1, 1].set_xlim((0, lim_max))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.BRD4_vs_MTHFD1.-log_pvalue.svg"), bbox_inches="tight", dpi=300)

    # Plot fold-change of BRD4 changing regions vs MTHFD1
    fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3))
    axis[0, 1].scatter(brd4['Fold'], mthfd1['Fold'], alpha=0.1, s=2, rasterized=True)
    lim_min = min(brd4['Fold'].min(), mthfd1['Fold'].min())
    lim_max = max(brd4['Fold'].max(), mthfd1['Fold'].max())
    axis[0, 1].plot((lim_min, lim_max), (lim_min, lim_max), linestyle="--", color="black", alpha=0.5)
    axis[0, 1].set_xlim((lim_min, lim_max))
    axis[0, 1].set_ylim((lim_min, lim_max))
    axis[0, 1].set_xlabel("BRD4")
    axis[0, 1].set_ylabel("MTHFD1")

    sns.distplot(brd4['Fold'], ax=axis[0, 0], kde=False)
    lim_max = max(brd4['Fold'].max(), mthfd1['Fold'].max())
    axis[0, 0].set_xlim((lim_min, lim_max))
    sns.distplot(mthfd1['Fold'], ax=axis[1, 1], kde=False)
    lim_max = max(mthfd1['Fold'].max(), mthfd1['Fold'].max())
    axis[1, 1].set_xlim((lim_min, lim_max))
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.BRD4_vs_MTHFD1.fold_change.svg"), bbox_inches="tight", dpi=300)
    
    # Combine data, plot heatmap of changes in both sets of changing regions
    brd4.index = brd4['seqnames'] + ":" + brd4['start'].astype(str) + "-" + brd4['end'].astype(str)
    mthfd1.index = mthfd1['seqnames'] + ":" + mthfd1['start'].astype(str) + "-" + mthfd1['end'].astype(str)
    occ = brd4.join(mthfd1, lsuffix="_brd4", rsuffix="_mthfd1")
    g = sns.clustermap(
        occ.loc[
            occ.sort_values("-log_pvalue_mthfd1").tail(500).index.tolist() + 
            occ.sort_values("-log_pvalue_brd4").tail(500).index.tolist(),
            ["Conc_DMSO_brd4", "Conc_dBET_brd4", "Conc_DMSO_mthfd1", "Conc_dBET_mthfd1"]
        ], cmap="inferno", vmin=0, robust=True, metric="correlation", yticklabels=False, figsize=(4, 4), rasterized=True
    )
    g.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500_both_factors.diffbind_values.clustermap.svg"), bbox_inches="tight", dpi=300)

    # Enrichment in the two sets of changing regions
    occ_enr = pd.DataFrame(columns=['ip', 'condition', 'region_set', 'enrichment', 'index'])
    for ip in ["BRD4", "MTHFD1"]:
        for condition in ['DMSO', "dBET"]:
            for region_set in ["MTHFD1", "BRD4"]:
                e = occ.loc[occ.sort_values("-log_pvalue_{}".format(region_set.lower())).tail(500).index.tolist(), "Conc_{}_{}".format(condition, ip.lower())].to_frame(name="enrichment").reset_index()
                e['ip'] = ip
                e['condition'] = condition
                e['region_set'] = region_set
                occ_enr = occ_enr.append(e)

    g = sns.FacetGrid(data=occ_enr, col="region_set", row="ip", hue="condition")
    g.map(sns.distplot, "enrichment")
    g.add_legend()
    g.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.distplot.svg"), bbox_inches="tight", dpi=300)
    g = sns.FacetGrid(data=occ_enr, col="region_set", row="ip", hue="condition")
    g.map(sns.violinplot, "enrichment", orientation="vertical")
    g.add_legend()
    g.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.violinplot.svg"), bbox_inches="tight", dpi=300)
    g = sns.FacetGrid(data=occ_enr, row=["region_set", "ip"], hue="condition")
    g.map(sns.violinplot, "enrichment", orientation="vertical")
    g.add_legend()
    g.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.violinplot.svg"), bbox_inches="tight", dpi=300)

    fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4))
    sns.violinplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "BRD4"], split=True, ax=axis[0])
    axis[0].set_title("BRD4 regions")
    sns.violinplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "MTHFD1"], split=True, ax=axis[1])
    axis[1].set_title("MTHFD1 regions")
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.violinplot.separate.svg"), bbox_inches="tight", dpi=300)
    fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4))
    sns.boxplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "BRD4"], ax=axis[0])
    axis[0].set_title("BRD4 regions")
    sns.boxplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "MTHFD1"], ax=axis[1])
    axis[1].set_title("MTHFD1 regions")
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.boxplot.separate.svg"), bbox_inches="tight", dpi=300)
    fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4))
    sns.barplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "BRD4"], ax=axis[0])
    axis[0].set_title("BRD4 regions")
    sns.barplot(x="ip", y="enrichment", hue="condition", data=occ_enr[occ_enr['region_set'] == "MTHFD1"], ax=axis[1])
    axis[1].set_title("MTHFD1 regions")
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.diffbind_values.barplot.separate.svg"), bbox_inches="tight", dpi=300)

    # LOLA enrichments
    library("LOLA")
    genome = "hg19"
    dbPath1 = paste0("/data/groups/lab_bock/shared/resources/regions/LOLACore/", genome)
    dbPath2 = paste0("/data/groups/lab_bock/shared/resources/regions/customRegionDB/", genome)
    bed_file <- "/home/arendeiro/projects/mthfd1/results/all_diff.bed"
    universe_file <- "/home/arendeiro/projects/mthfd1/results/chipseq_peaks/H3K27ac_WT_IgG/H3K27ac_WT_IgG_homer_peaks.bed"
    output_dir <- "/home/arendeiro/projects/mthfd1/results/diffbind_analysis.H3K27ac_peaks.top500.lola"

    regionDB = loadRegionDB(c(dbPath1, dbPath2))
    userSet <- LOLA::readBed(bed_file)
    userUniverse  <- LOLA::readBed(universe_file)
    lolaResults = runLOLA(list(userSet), userUniverse, regionDB, cores=8)
    writeCombinedEnrichment(lolaResults, outFolder=output_dir)

    # Motif analysis
    c = homer_motifs(
        "/home/arendeiro/projects/mthfd1/results/chipseq_peaks/H3K27ac_WT_IgG/H3K27ac_WT_IgG_homer_peaks.bed",
        "/home/arendeiro/projects/mthfd1/results/diffbind_analysis.H3K27ac_peaks.top500.lola", genome="hg19")
    !findMotifsGenome.pl /home/arendeiro/projects/mthfd1/results/all_diff.bed hg19r /home/arendeiro/projects/mthfd1/results/diffbind_analysis.H3K27ac_peaks.top500.lola -size 1000 -h -p 2 -len 8,10,12,14 -noknown

    !cut -f 1,2,3 /home/arendeiro/projects/mthfd1/results/all_diff.bed \
    > /home/arendeiro/projects/mthfd1/results/all_diff.3col.bed
    bed_to_fasta(
        "/home/arendeiro/projects/mthfd1/results/all_diff.3col.bed",
        "/home/arendeiro/projects/mthfd1/results/all_diff.fa")
    cmd = meme_ame(
        input_fasta="/home/arendeiro/projects/mthfd1/results//all_diff.fa",
        output_dir="/home/arendeiro/projects/mthfd1/results/diffbind_analysis.H3K27ac_peaks.top500.lola")

    # Enrichments
    # Upload BED file to Enrichr, download Reactome_2016 table
    df = pd.read_csv(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.condition.differential.enrichr.Reactome_2016.tsv"), sep="\t")
    df["Term"] = df["Term"].str.replace("_Homo .*", "")

    fig, axis = plt.subplots(1, figsize=(4, 4))
    top_n = 20
    sns.barplot(data=df.sort_values("Combined Score", ascending=False).head(top_n), y="Term", x="Combined Score", orient="horizontal", ax=axis)
    sns.despine(fig)
    fig.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.enrichr.Reactome.barplot.svg"), bbox_inches="tight", dpi=300)

    term_gene_assoc = df.set_index("Term")["Genes"].str.split(";").apply(pd.Series).stack().reset_index(level=0)
    term_gene_assoc = term_gene_assoc.rename(columns={0: "Gene"})
    term_gene_assoc['indent'] = 1
    piv = pd.pivot_table(
        term_gene_assoc,
        index="Term",
        columns="Gene",
        values="indent", fill_value=0)

    sig_piv = piv.loc[df.sort_values(["Combined Score"], ascending=False).head(top_n)['Term'], :]

    g = sns.clustermap(sig_piv.loc[:, sig_piv.sum() != 0], square=True, figsize=(4, 4), row_cluster=False)
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), rotation=90)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), rotation=0)
    g.savefig(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.top500.enrichr.Reactome.gene-pathway.heatmap.svg"), bbox_inches="tight", dpi=300)


    # # BRD4 vs MTHFD1
    # samples <- read.csv("metadata/diffBind_design.complete.csv")
    # samples["Peaks"] <- samples["H3K27acPeaks"]
    # diffbind_factors <- dba(sampleSheet=samples)
    # diffbind_factors <- dba.count(diffbind_factors)
    # diffbind_factors <- dba.contrast(diffbind_factors, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers=2)
    # diffbind_factors <- dba.analyze(diffbind_factors, bFullLibrarySize=TRUE, bCorPlot=FALSE, bParallel=TRUE)  # , method=DBA_ALL_METHODS
    # diffbind_factors.DB <- dba.report(diffbind_factors, bNormalized=TRUE, bCalled=TRUE, th=1.0, bUsePval=TRUE, fold=0)  # , method=DBA_ALL_METHODS
    # write.table(as.data.frame(diffbind_factors.DB), paste("results/diffbind_analysis.H3K27ac_peaks", "condition", "differential.csv", sep="."), sep=",", row.names=FALSE)


    # # QC plots
    # diff = pd.read_csv(os.path.join("results", "diffbind_analysis.H3K27ac_peaks.condition.differential.csv"))
    # diff['-log_pvalue'] = -np.log10(diff["p.value"])
    # fig, axis = plt.subplots(1, 3, figsize=(3 * 4, 1 * 4))
    # # Scatter
    # lim_max = max(diff["Conc_DMSO"].max(), diff["Conc_dBET"].max())
    # axis[0].scatter(diff["Conc_DMSO"], diff["Conc_dBET"], alpha=0.2, s=2, rasterized=True)
    # axis[0].set_xlim((0, lim_max))
    # axis[0].set_ylim((0, lim_max))
    # # MA
    # axis[1].scatter(diff["Conc"], diff["Fold"], alpha=0.2, s=2, rasterized=True)
    # # Volcano
    # lim_max = max(diff["Fold"].abs())
    # axis[2].scatter(diff["Fold"], diff["-log_pvalue"], alpha=0.2, s=2, rasterized=True)
    # sig = diff[diff['FDR'] < 0.05].index
    # axis[2].scatter(diff.loc[sig, "Fold"], diff.loc[sig, "-log_pvalue"], alpha=0.5, s=2, color="red", rasterized=True)
    # axis[2].set_xlim((-lim_max, lim_max))


    # Correlate occupancy and changes in gene expression

    ## get gene expression changes
    expr = get_gene_expression_data()

    ## get coverage in H3K27ac peaks
    analysis.set_consensus_sites(os.path.join("results", "chipseq_peaks", "H3K27ac_WT_IgG", "H3K27ac_WT_IgG_homer_peaks.bed"))
    analysis.coverage()
    coverage_rpm = analysis.normalize(method="total")
    analysis.get_peak_gene_annotation()
    coverage_rpm_gene = get_gene_level_accessibility(analysis, matrix="coverage_rpm").drop(['start', 'end'], axis=1)

    g = analysis.gene_annotation['gene_name'].str.split(",").apply(pd.Series).stack()
    g.index = g.index.droplevel(1)
    g.name = "gene_name"

    occ2 = occ.join(g).drop("gene_name", axis=1)
    occ2.index = occ.join(g).reset_index().set_index(['index', 'gene_name']).index
    occ2.columns = occ.columns
    occ3 = occ2.groupby(level="gene_name").mean()
    occ3 = occ3[occ3.columns[occ3.columns.str.contains("Conc_dBET|Conc_DMSO")]]

    ## For each comparison get occupancy in
    fig, axis = plt.subplots(3, 3, figsize=(4 * 3, 4 * 3), sharey=True, sharex=True)
    axis = axis.flatten()
    threshold = 0.05
    results = pd.DataFrame()
    for i, comparison in enumerate(expr['comparison_name'].unique()):
        if "HAP1" not in comparison:
            continue
        print(comparison)
        diff = expr[(expr['comparison_name'] == comparison) & (expr['q_value'] < threshold) & (expr['mean'] > 1)]

        up = diff[diff['log2(fold_change)'] > 0]
        down = diff[diff['log2(fold_change)'] < 0]

        if (comparison.endswith("_DMSO") or comparison.endswith("shPLKO")):
            print("reversing")
            down_ = up
            up_ = down
            down = down_
            up = up_

        up_genes = up['gene'].str.split(",").apply(pd.Series).stack().drop_duplicates()
        down_genes = down['gene'].str.split(",").apply(pd.Series).stack().drop_duplicates()

        # Get random set of genes with same size
        random_genes = (
            expr[(expr['comparison_name'] == comparison) & (expr['q_value'] > threshold) & (expr['mean'] > 1)]
            ['gene'].str.split(",").apply(pd.Series).stack().drop_duplicates()
        )
        random_up_genes = random_genes.sample(up_genes.shape[0] + down_genes.shape[0])
        random_down_genes = random_genes.sample(up_genes.shape[0] + down_genes.shape[0])

        # Use RPM values
        u = pd.melt(coverage_rpm_gene.ix[up_genes].dropna())
        d = pd.melt(coverage_rpm_gene.ix[down_genes].dropna())
        u['direction'] = 'up'
        d['direction'] = 'down'
        c = u.append(d)

        fig, axis = plt.subplots(1, 1, figsize=(6, 12))
        sns.boxplot(data=c, x="value", y="variable", hue="direction", orient="horizontal", ax=axis)
        sns.despine(fig)
        fig.savefig(os.path.join("results", "ChIP-RNA.diff.{}.{}.boxplot.svg".format(threshold, comparison)), dpi=300, bbox_inches="tight")

        # Use concentration values estimated by DiffBind
        u = pd.melt(occ3.ix[up_genes].dropna())
        d = pd.melt(occ3.ix[down_genes].dropna())
        u['direction'] = 'up'
        d['direction'] = 'down'
        r_u = pd.melt(occ3.ix[random_up_genes].dropna())
        r_d = pd.melt(occ3.ix[random_down_genes].dropna())
        r_u['direction'] = 'random_up'
        r_d['direction'] = 'random_down'
        c = u.append(d)
        c = c.append(r_u)
        c = c.append(r_d)
        sns.boxplot(data=c, x="value", y="variable", hue="direction", orient="horizontal", ax=axis[i])
        axis[i].set_title(comparison + "\nup: {}; down: {}".format(up_genes.shape[0], down_genes.shape[0]))

        for group in d['variable'].unique():
            from scipy.stats import ks_2samp
            a = d.loc[d['variable'] == group, 'value']
            b = u.loc[u['variable'] == group, 'value']
            t = ks_2samp(a, b)
            results = results.append(pd.Series([comparison, group, np.log2(a.mean() / b.mean()), t[0], t[1]]), ignore_index=True)

    sns.despine(fig)
    fig.savefig(os.path.join("results", "ChIP-RNA.diff.{}.all.comparisons.chipseq_diffbind_estimates.boxplot.svg".format(threshold)), dpi=300, bbox_inches="tight")

    results.columns = ['comparison', 'group', 'logfoldchange', 'stat', 'p-value']
    results.to_csv(os.path.join("results", "ChIP-RNA.diff.{}.all.comparisons.chipseq_diffbind_estimates.test_results.csv".format(threshold)))

    # # Correlation of fold-changes
    # comparison = "HAP1_DMSO_dBet6"
    # e = expr.loc[(expr['comparison_name'] == comparison), 'gene'].str.split(",").apply(pd.Series).stack()
    # e.index = e.index.droplevel(1)
    # e.name = "gene_name"
    # e = expr.loc[(expr['comparison_name'] == comparison), :].join(e)
    # e = e.set_index("gene_name")

    # o = occ2.groupby(level="gene_name").mean().drop_duplicates()

    # e = e.loc[o.index, 'log2(fold_change)'].dropna()
    # o = o.loc[e.index, 'Fold_mthfd1'].dropna()

    # fig, axis = plt.subplots(1)
    # axis.scatter(o, e, s=3, alpha=0.3, rasterized=True)
    # sns.despine(fig)
    # fig.savefig(os.path.join("results", "ChIP-RNA.diff.mthfd1.fold_changes.scatter.svg"), dpi=300, bbox_inches="tight")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)


def ci(x, ci=.95, N=1000):
    from scipy.stats import norm
    i = norm.interval(ci, loc=x.mean(), scale=x.std() / np.sqrt(N))
    return i[1] - i[0]

