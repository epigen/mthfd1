
import os
import pandas as pd


from peppy import Project, Sample
from ngs_toolkit.general import unsupervised_analysis
from ngs_toolkit.rnaseq import RNASeqAnalysis
from ngs_toolkit.rnaseq import knockout_plot
from ngs_toolkit.graphics import add_colorbar_to_axis
from ngs_toolkit.general import fix_batch_effect_limma, subtract_principal_component
from ngs_toolkit.general import differential_enrichment, collect_differential_enrichment, plot_differential_enrichment

import numpy as np
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt


from statsmodels.stats.multitest import multipletests
import scipy


# Set settings
pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("Paired"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


# Start project
prj = Project(os.path.join("metadata", "project_config.yaml"))
analysis = RNASeqAnalysis(name="mthfd1.rnaseq", prj=prj, samples=[s for s in prj.samples if s.library == "RNA-seq"])


# Get RNA-seq quantification from Cufflinks estimates
comps = pd.read_table(os.path.join("metadata", "rnaseq_comparisons.tsv")).squeeze().tolist()
diff = pd.DataFrame()
diff_files = [
    "/scratch/lab_bsf/projects/BSA_0197_Compound_shRNA_/hg38/rnaseq_cuffdiff_global/gene_exp.diff",
    "/scratch/lab_bsf/projects/BSA_0226_HAP1_lines/hg38/rnaseq_cuffdiff_global/gene_exp.diff"]
for f in diff_files:
    print(f)
    d = pd.read_csv(f, sep='\t')
    d['comparison_name'] = d['sample_2'] + "_vs_" + d['sample_1']
    d['mean_expression'] = d.loc[:, ['value_2', 'value_1']].mean(axis=1)
    d = d.loc[d['comparison_name'].isin(comps), :]
    diff = diff.append(d)
diff = diff.set_index("gene")
# Invert direction for certain comparisons
diff.loc[diff['sample_2'].str.contains("DMSO|WT|shPLKO"), 'log2FoldChange'] = -diff.loc[
    diff['sample_2'].str.contains("DMSO|WT|shPLKO"), 'log2FoldChange']
# rename columns to match DESeq2
diff = diff.rename(columns={
    "mean_expression": "baseMean",
    "log2(fold_change)": "log2FoldChange",
    "p_value": "pvalue",
    "q_value": "padj"})
analysis.differential_results = diff
# save selected comparisons
diff.to_csv(os.path.join(
    analysis.results_dir, "differential_analysis_RNA-seq",
    analysis.name + ".cuffdiff_analysis.selected_comparisons.csv"))

diff = analysis.differential_results
one = diff.pivot_table(index="gene", columns="sample_1", values="value_1")
two = diff.pivot_table(index="gene", columns="sample_2", values="value_2")
three = one.join(two).sort_index(axis=1)
exp = np.log2(1 + three)
exp = exp[(~exp.index.str.startswith("-")) & (~exp.index.str.contains(","))]

exp_norm = exp # normalize_quantiles_p(exp)
analysis.coverage_annotated = exp_norm - exp_norm.min().min()

analysis.coverage_annotated.to_csv(os.path.join(
    analysis.results_dir,
    analysis.name + ".cufflinks.quantification.fpkm.qnorm.csv"))

# variables = ['sample_name', 'cell_line', 'condition', 'timepoint', 'knockout_clone', 'replicate', 'experiment_name']
# analysis.annotate_with_sample_metadata(attributes=variables, quant_matrix="coverage_annotated")

analysis.expression = analysis.coverage_annotated.loc[~(analysis.coverage_annotated.sum(axis=1) == 0)].dropna().drop_duplicates()


analysis.expression = analysis.expression.rename(
    columns={"HAP1_WT": "HAP1_WT-WT", "C8": "HAP1_MTHFD1KO-C8", "D3": "HAP1_MTHFD1KO-D3"})

# Add just minimal sample info (name)
samples = list()
for i, col in enumerate(analysis.expression.columns):
    samples.append(
        Sample(
            pd.Series([col, col.split("_")[0], "_".join(col.split("_")[1:]), "b1" if "-" in col else "b2"],
                      index=['sample_name', 'cell_line', 'perturbation', 'batch'])))

analysis._samples = samples
analysis.expression.columns = pd.MultiIndex.from_arrays(
    [
        [s.name for s in analysis._samples],
        [s.cell_line for s in analysis._samples],
        [s.perturbation for s in analysis._samples],
        [s.batch for s in analysis._samples]],
    names=['sample_name', "cell_line", "perturbation", "batch"])
analysis.to_pickle()


# Unsupervised analysis
unsupervised_analysis(
    analysis,
    samples=analysis._samples,
    data_type="RNA-seq",
    quant_matrix="expression", attributes_to_plot=['cell_line', 'perturbation', 'batch'])

for cell_line in ["HAP1", "A549", "K562"]:
    unsupervised_analysis(
        analysis,
        samples=[s for s in analysis._samples if s.cell_line == cell_line],
        data_type="RNA-seq",
        quant_matrix="expression", attributes_to_plot=['perturbation', 'batch'],
        plot_prefix=cell_line)

    if cell_line == "HAP1":
        unsupervised_analysis(
            analysis,
            samples=[s for s in analysis._samples if (s.cell_line == cell_line) and ("-" not in s.name)],
            data_type="RNA-seq",
            quant_matrix="expression", attributes_to_plot=['perturbation', 'batch'],
            plot_prefix="HAP1_noKO")  # , plot_group_centroids=False, legends=True)

# Try to use limma to remove batch effect
fix = fix_batch_effect_limma(
    analysis.expression.loc[:, analysis.expression.columns.get_level_values("cell_line") == "HAP1"],
    batch_variable="batch", covariates=['cell_line'])
analysis.expression_limmafixed = fix
unsupervised_analysis(
    analysis,
    samples=[s for s in analysis._samples if s.name in analysis.expression_limmafixed.columns],
    data_type="RNA-seq",
    quant_matrix="expression_limmafixed", attributes_to_plot=['perturbation', 'batch'],
    plot_prefix="HAP1_limmafixed", test_pc_association=False)

# Try to remove first PC
to_fix = analysis.expression.loc[:, analysis.expression.columns.get_level_values("cell_line") == "HAP1"]
pcafix = subtract_principal_component(
    to_fix.T * 1000, pc=1, norm=False, plot=True,
    plot_name="PCA_based_batch_correction.svg", pcs_to_plot=6)
analysis.expression_pcafixed = pcafix.T
unsupervised_analysis(
    analysis,
    samples=[s for s in analysis._samples if s.name in analysis.expression_pcafixed.columns],
    data_type="RNA-seq",
    quant_matrix="expression_pcafixed", attributes_to_plot=['perturbation', 'batch'],
    plot_prefix="HAP1_pcafixed", test_pc_association=False)

# Knockout plot
chrom_prots = pd.read_csv(os.path.join("metadata", "chromatin_proteins.csv"), squeeze=True)
knockout_plot(
    analysis, knockout_genes=np.random.choice(chrom_prots.tolist(), 50).tolist() + ['MTHFD1'],
    # ["ACTB", "TUBB", "E2F1", "E2F4", "RB1", "CCND1", "CCND2", "MTHFD1", "MTHFD2"]
    expression_matrix=analysis.expression)


# Supervised analysis
diff = pd.read_csv(os.path.join(
    analysis.results_dir, "differential_analysis_RNA-seq",
    analysis.name + ".cuffdiff_analysis.selected_comparisons.csv"), index_col=0)

# Data as is from Cuffdiff
# # diagnostic plots
for comp in comps:
    diff3 = diff.loc[diff['comparison_name'] == comp, :].reset_index()
    s = diff3.loc[
        (diff3['padj'] < 0.05) &
        (diff3['baseMean'] >= 1)].index

    if ("C8" in comp) or ("D3" in comp):
        s = diff3.loc[
            (diff3['padj'] < 0.005) &
            (diff3['log2FoldChange'].abs() > 1) &
            (diff3['baseMean'] >= 3)].index
    print(comp, len(s))

    fig, axis = plt.subplots(1, 3, figsize=(3 * 3.5, 1 * 3))  # , subplot_kw={"aspect": 'equal'})
    fig.suptitle(comp)
    # Scatter
    axis[0].scatter(
        np.log2(1 + diff3['value_2']),
        np.log2(1 + diff3['value_1']),
        alpha=0.1, s=1, rasterized=True, color="grey")
    collection = axis[0].scatter(
        np.log2(1 + diff3.loc[s, 'value_2']),
        np.log2(1 + diff3.loc[s, 'value_1']),
        alpha=0.2, s=2, rasterized=True, c=-np.log10(diff3.loc[s, 'pvalue']), vmin=0)
    axis[0].plot(axis[0].get_ylim(), axis[0].get_ylim(), linestyle="--", color="black", alpha=0.5)
    axis[0].set_xlabel(diff3['sample_2'].unique()[0])
    axis[0].set_ylabel(diff3['sample_1'].unique()[0])
    add_colorbar_to_axis(collection, label="-log10(p-value)")
    # MA plot
    axis[1].scatter(
        np.log2(1 + diff3['baseMean']),
        diff3['log2FoldChange'],
        alpha=0.1, s=1, rasterized=True, color="grey")
    collection = axis[1].scatter(
        np.log2(1 + diff3.loc[s, 'baseMean']),
        diff3.loc[s, 'log2FoldChange'],
        alpha=0.2, s=2, rasterized=True, c=-np.log10(diff3.loc[s, 'pvalue']), vmin=0)
    axis[1].axhline(0, linestyle="--", color="black", alpha=0.5)
    axis[1].set_xlabel("Mean expression")
    axis[1].set_ylabel("log2FoldChange")
    m = diff3['log2FoldChange']
    m = m[(m != np.inf) & (m != -np.inf)].abs().max()
    axis[1].set_ylim((-m, m))
    add_colorbar_to_axis(collection, label="-log10(p-value)")

    # Volcano
    axis[2].scatter(
        diff3['log2FoldChange'],
        -np.log10(diff3['pvalue']),
        alpha=0.1, s=1, rasterized=True, color="grey")
    collection = axis[2].scatter(
        diff3.loc[s, 'log2FoldChange'],
        -np.log10(diff3.loc[s, 'pvalue']),
        alpha=0.2, s=2, rasterized=True, c=np.log2(diff3.loc[s, 'baseMean']), vmin=0)
    axis[2].set_xlabel("log2FoldChange")
    axis[2].set_ylabel("-log10(p value)")
    axis[2].axvline(0, linestyle="--", color="black", alpha=0.5)
    axis[2].set_xlim((-m, m))
    add_colorbar_to_axis(collection, label="Mean Expression")

    fig.savefig(os.path.join(analysis.results_dir, "differential_analysis_RNA-seq", "{}.{}.comparison.{}.svg".format(
        'mthfd1', "fpkm", comp)), bbox_inches="tight", dpi=300, tight_layout=True)


# # Function of differential genes
diff = analysis.differential_results
# Data as is from Cuffdiff
diff_sig = pd.DataFrame()
comps = pd.read_table(os.path.join("metadata", "rnaseq_comparisons.tsv")).squeeze().tolist()
for comp in comps:
    diff3 = diff.loc[diff['comparison_name'] == comp, :]
    s = diff3.loc[
        (diff3['padj'] < 0.05) &
        (diff3['baseMean'] >= 1)]

    if ("C8" in comp) or ("D3" in comp):
        s = diff3.loc[
            (diff3['padj'] < 0.005) &
            (diff3['log2FoldChange'].abs() > 1) &
            (diff3['baseMean'] >= 3)]
    print(comp, s.shape[0])
    s = s.loc[(~s.index.str.startswith("-")) & (~s.index.str.contains(","))]
    diff_sig = diff_sig.append(s)

diff_sig.index.name = "gene_name"
one = diff.pivot_table(index="gene", columns="sample_1", values="value_1")
two = diff.pivot_table(index="gene", columns="sample_2", values="value_2")
three = one.join(two).sort_index(axis=1)
analysis.expression = np.log2(1 + three)

enrichment_dir = os.path.join("results", "differential_analysis_RNA-seq", "enrichment")

differential_enrichment(
    analysis, diff_sig, data_type="RNA-seq",
    output_dir=enrichment_dir,
    genome="hg38", directional=True)


collect_differential_enrichment(
    diff, directional=True, data_type="RNA-seq",
    output_dir=enrichment_dir)

enrichment_table = pd.read_csv(os.path.join(enrichment_dir, "differential_analysis.enrichr.csv"))

plot_differential_enrichment(
    enrichment_table,
    "enrichr",
    data_type="RNA-seq",
    output_dir=enrichment_dir,
    output_prefix="differential_analysis",
    direction_dependent=True,
    top_n=5)


# # Start relating with MTHFD1 and BRD4 binding

# # # First, observe values of expression at steady state
exp = analysis.expression.loc[:, analysis.expression.columns.get_level_values("cell_line") == cell_line]
# # # # separate expressed from non-expressed genes
fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 2.5))
for ax in axis:
    sns.distplot(
        exp.mean(axis=1),
        ax=ax, kde=False, color=sns.color_palette("colorblind")[0], norm_hist=True)
    ax.axvline(0.5, linestyle="--", color="black")
    ax.axvline(2, linestyle="--", color="black")
    ax.set_xlabel("log2(mean expression)")
    ax.set_ylabel("Transcript density")
axis[1].set_yscale("log")
fig.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.expressed_vs_nonexpressed.from_diff.distplot.svg"), dpi=300, bbox_inches="tight")

g = exp.mean(axis=1)
not_expressed = g[g <= 0.5].index
lowly_expressed = g[(g < 2) & (g >= 0.5)].index
highly_expressed = g[g >= 2].index

# # # use GREAT peak - region assignments # get bound genes
bound_BRD = pd.read_csv(os.path.join(
        analysis.results_dir, "differential_analysis_ChIP-seq", "BRD4.greatExportAll.parsed.csv"))
bound_BRD = bound_BRD['Genes'].str.split(",").apply(pd.Series).stack().drop_duplicates().reset_index(drop=True)
bound_MTHFD = pd.read_csv(os.path.join(
        analysis.results_dir, "differential_analysis_ChIP-seq", "MTHFD1.greatExportAll.parsed.csv"))
bound_MTHFD = bound_MTHFD['Genes'].str.split(",").apply(pd.Series).stack().drop_duplicates().reset_index(drop=True)
bound_both = bound_BRD[bound_BRD.isin(bound_MTHFD.tolist())]

res = pd.DataFrame()
for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
    print(comp, label)
    d = exp.loc[ip].dropna().reset_index()
    print(comp, d.shape[0])
    d2 = d.melt(id_vars="gene", var_name="comparison_name", value_name="expression")
    d2['bound'] = label
    res = res.append(d2)

    print(comp, label)
    d3 = exp.loc[np.random.choice(exp.index, d2.shape[0])].dropna().reset_index()
    print(comp, d.shape[0])
    d3 = d3.melt(id_vars="gene", var_name="comparison_name", value_name="expression")
    d3['bound'] = "random of size " + label
    res = res.append(d3)

all_ips = bound_BRD.append(bound_MTHFD).append(bound_both)
for expr, label in [(not_expressed, "not expressed"), (lowly_expressed, "lowly expressed"), (highly_expressed, "highly expressed")]:

    # get genes which are not bound, but which are in different levels of expression
    d = (
        exp.loc[
            ~exp.index.isin(all_ips.tolist()) &
            exp.index.isin(expr.tolist())]
        .dropna()
        .reset_index()
        .melt(id_vars="gene", var_name="comparison_name", value_name="expression")
    )
    d['bound'] = "Neither " + label
    res = res.append(d)

grid = sns.catplot(
    data=res,
    y="comparison_name", x="expression", hue="bound",
    orient="horiz", kind="box", showfliers=False, height=18, aspect=0.3)
grid.axes.flatten()[0].set_xlabel("Expression log2(RPKM)")

grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.steady_state_expression.from_diff.on_bound_genes.boxplot.svg"),
    dpi=300, bbox_inches="tight")


# Match data to MTHFD1 binding


# # Are bound genes loosing expression upon KO, sh or drug?
total_genes = len(diff[(~diff.index.str.startswith("-")) & (~diff.index.str.contains(","))].index.unique())
res = pd.DataFrame()
for comp in diff['comparison_name'].unique():
    for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
        print(comp, label)
        d = diff[diff['comparison_name'] == comp]

        # if ("C8" in comp) or ("D3" in comp):
        #     d_sig = d.loc[
        #         (d['padj'] < 0.005) &
        #         (d['log2FoldChange'].abs() > 1) &
        #         (d['baseMean'] >= 3)]
        # else:
        d_sig = d[d['padj'] < 0.05]

        # cleanup
        d_sig = d_sig[(~d_sig.index.str.startswith("-")) & (~d_sig.index.str.contains(","))]

        # Test overlap
        a = d_sig.index.isin(ip).sum()  # intersection of differential and bound
        b = d_sig.shape[0] - a  # differential but not bound
        c = (~ip.isin(d_sig.index.tolist())).sum()  # bound not differential
        d = total_genes - sum([a, b, c])  # rest of genome
        s, p = scipy.stats.fisher_exact([[a, b], [c, d]], alternative="greater")
        res = res.append(
            pd.Series([label, a, b, c, d, s, p],
                      index=['ip', 'a_intersection', 'b_diffexp_only', 'c_bound_only', 'd_neither', 'log_odds', 'fisher_p_value'],
                      name=comp))

res.index.name = 'comparison_name'

res['fisher_q_value'] = multipletests(res['fisher_p_value'], method='fdr_bh')[1]
res.to_csv(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.fisher.csv"))


res['log_fisher_p_value'] = -np.log10(res['fisher_p_value'].replace(0, 1e-300))
res['log_fisher_q_value'] = -np.log10(res['fisher_q_value'].replace(0, 1e-300))
grid = sns.catplot(data=res.reset_index(), y="comparison_name", x="log_fisher_p_value", hue="ip", orient="horiz")
grid.savefig(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.fisher.log_fisher_p_value.stripplot.svg"))
grid = sns.catplot(data=res.reset_index(), y="comparison_name", x="log_fisher_q_value", hue="ip", orient="horiz")
grid.savefig(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.fisher.log_fisher_q_value.stripplot.svg"))
grid = sns.catplot(data=res.reset_index(), y="comparison_name", x="log_odds", hue="ip", orient="horiz")
grid.savefig(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.fisher.log_odds.stripplot.svg"))

fig, axis = plt.subplots(1, 1, figsize=(12, 8))
for comp in diff['comparison_name'].unique():
    axis.scatter(res.loc[comp, 'log_odds'], res.loc[comp, "log_fisher_p_value"], label=comp)
    for n, d in res.loc[comp].iterrows():
        axis.text(d['log_odds'], d["log_fisher_p_value"], s=n + " - " + d['ip'])
axis.legend(loc='lower left')
fig.savefig(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.fisher.scatter.svg"))


# Now look at the quantitative differences

# # See distribution of absolute changes

# # all genes
diff2 = diff.pivot_table(index="gene", columns="comparison_name", values="log2FoldChange")
diff2 = diff2[diff2 != np.inf]
diff2 = diff2[diff2 != -np.inf]
diff2["bound"] = False
for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
    diff2.loc[
        (diff2.index.isin(ip.tolist())), "bound"] = label

diff2.loc[:, diff2.dtypes == float] = diff2.loc[:, diff2.dtypes == float].abs()
grid = sns.catplot(
    data=diff2.melt(id_vars=['bound']),
    y='comparison_name', x='value', hue='bound',
    orient="horiz", kind="box", showfliers=False, height=16, aspect=0.5)
grid.axes.flatten()[0].set_xlabel("absolute log2(fold-change)")
grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.bound.absolute_fold_changes.boxplot.svg"), dpi=300, bbox_inches="tight")

# # separate expressed from non-expressed genes
diff2 = diff.pivot_table(index="gene", columns="comparison_name", values="log2FoldChange")
diff2 = diff2[diff2 != np.inf]
diff2 = diff2[diff2 != -np.inf]
diff2["bound"] = False
for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
    diff2.loc[
        (diff2.index.isin(ip.tolist())) &
        (diff2.index.isin(highly_expressed.tolist())), "bound"] = label + " + high expression"
    diff2.loc[
        (diff2.index.isin(ip.tolist())) &
        (~diff2.index.isin(highly_expressed.tolist())), "bound"] = label + " + low expression"
diff2.loc[
    (~diff2.index.isin(bound_BRD.tolist() + bound_MTHFD.tolist() + bound_both.tolist())) &
    (diff2.index.isin(highly_expressed.tolist())), "bound"] = "Not bound + high expression"
diff2.loc[
    (~diff2.index.isin(bound_BRD.tolist() + bound_MTHFD.tolist() + bound_both.tolist())) &
    (~diff2.index.isin(highly_expressed.tolist())), "bound"] = "Not bound + low expression"

vc = diff2['bound'].value_counts()
fig, axis = plt.subplots(1, figsize=(4, 4))
sns.barplot(vc, vc.index, ax=axis, orient="horiz")
axis.set_xscale('log')
axis.set_xlabel("Number of genes")
fig.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.bound.gene_counts.barplot.svg"), dpi=300, bbox_inches="tight")


diff2.loc[:, diff2.dtypes == float] = diff2.loc[:, diff2.dtypes == float].abs()
grid = sns.catplot(
    data=diff2.melt(id_vars=['bound']),
    y='comparison_name', x='value', hue='bound',
    orient="horiz", kind="box", showfliers=False, height=16, aspect=0.5)
grid.axes.flatten()[0].set_xlabel("absolute log2(fold-change)")
grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.bound.absolute_fold_changes.expression_dependent.boxplot.svg"), dpi=300, bbox_inches="tight")


# # See distribution of relative changes
res = pd.DataFrame()
res2 = pd.DataFrame()
for comp in diff['comparison_name'].unique():
    for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
        print(comp, label)
        d = diff[diff['comparison_name'] == comp]
        d = d[(~d.index.str.startswith("-")) & (~d.index.str.contains(","))]

        # get fold changes
        # genes bound
        a = d.loc[d.index.isin(ip), "log2FoldChange"]
        # genes not bound
        b = d.loc[~d.index.isin(ip), "log2FoldChange"]
        a = a[~((a == np.inf) | (a == -np.inf))]
        b = b[~((b == np.inf) | (b == -np.inf))]
        a_up = a[a > 0]
        a_down = a[a < 0]
        b_up = b[b > 0]
        b_down = b[b < 0]

        s, p = scipy.stats.mannwhitneyu(a, b)

        # # fig, axis = plt.subplots(1, 1, figsize=(3, 3))
        # # sns.distplot(a, kde=False, label=label, ax=axis, color="red")
        # # sns.distplot(b, kde=False, label="not " + label, ax=axis, color="blue")
        # # axis.legend()
        res = res.append(
            pd.Series([label, a.mean(), b.mean(), s, p],
                      index=['ip', 'bound', 'not_bound', 'stat', 'mannwhitneyu_p_value'],
                      name=comp))
        a = a.to_frame()
        a['ip'] = label
        a['bound'] = True
        a['comparison_name'] = comp
        a['direction'] = "all"
        a_up = a_up.to_frame()
        a_up['ip'] = label
        a_up['bound'] = True
        a_up['comparison_name'] = comp
        a_up['direction'] = "up"
        a_down = a_down.to_frame()
        a_down['ip'] = label
        a_down['bound'] = True
        a_down['comparison_name'] = comp
        a_down['direction'] = "down"
        b = b.to_frame()
        b['ip'] = label
        b['bound'] = False
        b['comparison_name'] = comp
        b['direction'] = "all"
        b_up = b_up.to_frame()
        b_up['ip'] = label
        b_up['bound'] = False
        b_up['comparison_name'] = comp
        b_up['direction'] = "up"
        b_down = b_down.to_frame()
        b_down['ip'] = label
        b_down['bound'] = False
        b_down['comparison_name'] = comp
        b_down['direction'] = "down"
        res2 = res2.append(a)
        res2 = res2.append(b)
        res2 = res2.append(a_up)
        res2 = res2.append(a_down)
        res2 = res2.append(b_up)
        res2 = res2.append(b_down)

res['mannwhitneyu_q_value'] = multipletests(res['mannwhitneyu_p_value'], method='fdr_bh')[1]
res.to_csv(os.path.join(analysis.results_dir, "differential_RNA-ChIP_overlap.mannwhitneyu.csv"))


# # plot distributions
# # # only mean
p = res[['ip', 'bound', 'not_bound']].reset_index().melt(id_vars=['index', 'ip'])
p['group'] = p['index'] + " - " + p['ip']
grid = sns.catplot(
    data=p,
    y="group", x="value", hue="variable", height=8, aspect=1)
grid.ax.set_xlabel("log2(fold-change)")
grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.bound_{}.fold_changes.stripplot.svg"
        .format(label.replace(" ", "_"))), dpi=300)

# # # whole distributions
# res2['bound'] = res2['bound'].astype(str).replace("False", "not bound").replace("True", "bound")
res2['group'] = res2['comparison_name'] + " - " + res2['ip'] + " - " + res2['direction']
grid = sns.catplot(
    data=res2,
    y="group", x="log2FoldChange", hue="bound", height=32, aspect=0.31, kind="box", showfliers=False)
grid.ax.set_xlabel("log2(fold-change)")
grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.bound_{}.fold_changes.boxenplot.svg"
        .format(label.replace(" ", "_"))), dpi=300)

# # # visualize fold-changes in bound genes
for ip, label in [(bound_BRD, "BRD4"), (bound_MTHFD, "MTHFD1"), (bound_both, "BRD4 and MTHFD1")]:
    diff2 = (
        diff.loc[
            (diff.index.isin(ip)) &
            (diff['comparison_name'] != "HAP1_MTX_6h_vs_HAP1_DMSO_6h")]
        .reset_index()
        .pivot_table(index="gene", columns="comparison_name", values="log2FoldChange"))
    diff3 = diff2[diff2 != np.inf]
    diff3 = diff3[diff3 != -np.inf]
    diff3 = diff3[diff3 != np.nan]
    diff4 = diff3.fillna(0)
    diff4 = diff4.loc[diff4.sum(1) != 0]
    grid = sns.clustermap(
        diff4.T, cmap="RdBu_r", center=0, robust=True, metric="correlation",
        cbar_kws={"label": "log2(fold-change)"}, rasterized=True)
    grid.ax_heatmap.set_xlabel("Genes bound by BRD4 and MTHFD1 (n={})".format(diff4.shape[0]))
    grid.savefig(
        os.path.join(
            analysis.results_dir,
            "differential_RNA-ChIP_overlap.bound_{}.fold_changes.clustermap.svg"
            .format(label.replace(" ", "_"))), dpi=300)

    diff5 = diff4.loc[:, diff4.columns.str.contains("HAP1")]
    diff5 = diff5.loc[diff5.sum(1) != 0]
    grid = sns.clustermap(
        diff5.T, cmap="RdBu_r", center=0, robust=True, metric="correlation",
        cbar_kws={"label": "log2(fold-change)"}, rasterized=True)
    grid.ax_heatmap.set_xlabel("Genes bound by BRD4 and MTHFD1 (n={})".format(diff4.shape[0]))
    grid.savefig(
        os.path.join(
            analysis.results_dir,
            "differential_RNA-ChIP_overlap.bound_{}.fold_changes.clustermap.only_HAP1.svg"
            .format(label.replace(" ", "_"))), dpi=300)

#


# Now look into the amount of ChIP-seq signal in the genes with deregulated expression
df = pd.read_csv(os.path.join("results", "differential_analysis_ChIP-seq", "diffBind.results.csv"))


# # get locations to annotate with genes
ids = ['Chr', 'Start', 'End', 'Conc', 'FDR', 'Fold', 'comparison_name', 'p-value']
df2 = df.melt(id_vars=ids).pivot_table(index=ids[:3], columns='variable', values='value')
locs = os.path.join("results", "differential_analysis_ChIP-seq", "diffBind.results.index.bed")
df2.index.to_frame().to_csv(locs, index=False, header=False, sep="\t")
annot = os.path.join("results", "differential_analysis_ChIP-seq", "diffBind.results.index.annotated.bed")
if not os.path.exists(annot):
    os.popen("annotatePeaks.pl {} hg19 > {}".format(locs, annot))

peaks_annotated = pd.read_csv(annot, sep="\t")
peaks_annotated = peaks_annotated.set_index(ids[:3])['Gene Name'].sort_index()

chip_per_gene = (
    df2.reset_index(drop=True).join(
        peaks_annotated.reset_index(drop=True))
    .groupby('Gene Name').mean())

# # observe signal for different sets of genes
chip_res = pd.DataFrame()
comps = pd.read_table(os.path.join("metadata", "rnaseq_comparisons.tsv")).squeeze().tolist()
for comp in comps:
    diff2 = diff.loc[diff['comparison_name'] == comp, :]
    # Up-regulated
    s = diff2.loc[
        (diff2['padj'] < 0.05)
        & (diff2['baseMean'] >= 1)
        & (diff2['log2FoldChange'] > 0)
        ].index
    if ("C8" in comp) or ("D3" in comp):
        s = diff2.loc[
            (diff2['padj'] < 0.005) &
            (diff2['log2FoldChange'].abs() > 1) &
            (diff2['baseMean'] >= 3)
            & (diff2['log2FoldChange'] > 0)].index
    print(comp, len(s))

    try:
        c = chip_per_gene.loc[s].dropna().melt()
    except KeyError:
        continue
    c['comparison_name'] = comp + ", up-regulated"
    chip_res = chip_res.append(c)

    # Down-regulated
    s = diff2.loc[
        (diff2['padj'] < 0.05)
        & (diff2['baseMean'] >= 1)
        & (diff2['log2FoldChange'] < 0)
        ].index
    if ("C8" in comp) or ("D3" in comp):
        s = diff2.loc[
            (diff2['padj'] < 0.005) &
            (diff2['log2FoldChange'].abs() > 1) &
            (diff2['baseMean'] >= 3)
            & (diff2['log2FoldChange'] < 0)].index
    print(comp, len(s))

    try:
        c = chip_per_gene.loc[s].dropna().melt()
    except KeyError:
        continue
    c['comparison_name'] = comp + ", down-regulated"
    chip_res = chip_res.append(c)

order = ['Conc_BRD4_WT', 'Conc_BRD4_MTHFD1KO', 'Conc_BRD4_DMSO', 'Conc_BRD4_dBET',
         'Conc_MTHFD1_WT', 'Conc_MTHFD1_MTHFD1KO', 'Conc_MTHFD1_DMSO', 'Conc_MTHFD1_dBET']
grid = sns.catplot(
    data=chip_res,
    y="comparison_name", x="value", hue="variable", hue_order=order,
    orient="horiz", kind="box", showfliers=False, height=22, aspect=0.3)
grid.axes.flatten()[0].set_xlabel("ChIP-seq signal")
grid.savefig(
    os.path.join(
        analysis.results_dir,
        "differential_RNA-ChIP_overlap.ChIP-seq_signal_on_diff_genes.boxplot.svg"),
    dpi=300, bbox_inches="tight")
