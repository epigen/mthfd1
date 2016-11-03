#!/usr/bin/env python

"""
This script calls peaks for the mthfd1 project.
"""

import os
import pandas as pd
from looper.models import Project
from pypiper import NGSTk
import textwrap


def macs2CallPeaks(treatmentBams, outputDir, sampleName, genome, controlBams=None, broad=False, paired=False):
    """
    Use MACS2 to call peaks.
    """
    sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}

    cmd = "macs2 callpeak -t {0}".format(treatmentBams if type(treatmentBams) is str else " ".join(treatmentBams))
    if controlBams is not None:
        cmd += " -c {0}".format(controlBams if type(controlBams) is str else " ".join(controlBams))
    if paired:
        cmd += "-f BAMPE "
    if not broad:
        cmd += " --fix-bimodal --extsize 180 --bw 200"
    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3"
    cmd += " -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)

    return cmd


def macs2CallPeaks_weak(treatmentBams, outputDir, sampleName, genome, controlBams=None, broad=False):
    """
    Use MACS2 to call peaks.
    """

    sizes = {"hg38": 2.7e9, "hg19": 2.7e9, "mm10": 1.87e9, "dr7": 1.412e9}

    if not broad:
        cmd = "macs2 callpeak -t {0}".format(treatmentBams if type(treatmentBams) is str else " ".join(treatmentBams))
        if controlBams is not None:
            cmd += " -c {0}".format(controlBams if type(controlBams) is str else " ".join(controlBams))
        cmd += " --nomodel --extsize 73 --pvalue 1e-3 -g {0} -n {1} --outdir {2}".format(sizes[genome], sampleName, outputDir)
    else:
        # Parameter setting for broad factors according to Nature Protocols (2012)
        # Vol.7 No.9 1728-1740 doi:10.1038/nprot.2012.101 Protocol (D) for H3K36me3
        cmd = "macs2 callpeak -t {0}".format(treatmentBams if type(treatmentBams) is str else " ".join(treatmentBams))
        if controlBams is not None:
            cmd += " -c {0}".format(controlBams if type(controlBams) is str else " ".join(controlBams))
        cmd += " --broad --nomodel --extsize 73 --pvalue 1e-3 -g {0} -n {1} --outdir {2}".format(
            sizes[genome], sampleName, outputDir
        )

    return cmd


def call_peaks(samples, name, controls=None):
    """
    Call H3K27ac peaks for pairs of sample and igg control.
    """
    output_dir = os.path.join("data", "peaks", name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    job_file = os.path.join(output_dir, "call_peaks.sh")
    # prepare slurm job header
    cmd = tk.slurm_header(
        name + "_MACS2", os.path.join(output_dir, "call_peaks.log"),
        cpus_per_task=8,
        queue="shortq")

    #
    cmd += macs2CallPeaks_weak(
        treatmentBams=[sample.filtered for sample in samples],
        controlBams=[sample.filtered for sample in controls] if controls is not None else None,
        outputDir=os.path.join("data", "peaks", name),
        sampleName=name,
        genome=[sample.genome for sample in samples][0],
        broad=False)

    # write job to file
    with open(job_file, 'w') as handle:
        handle.writelines(textwrap.dedent(cmd))

    tk.slurm_submit_job(job_file)


# Start project
prj = Project("metadata/project_config.yaml")
prj.add_sample_sheet()
[sample.get_read_type() for sample in prj.samples]

# Read comparison table
comparisons = pd.read_csv(os.path.join("metadata", "comparison_table.csv"))


tk = NGSTk()

# call peaks as defined in comparison table
for comparison in comparisons["comparison"].drop_duplicates():
    sub = comparisons[comparisons["comparison"] == comparison]

    treatment_samples = [s for s in prj.samples if s.name in sub[sub["comparison_side"] == 1]['sample'].tolist()]
    control_samples = [s for s in prj.samples if s.name in sub[sub["comparison_side"] == 0]['sample'].tolist()]

    call_peaks(treatment_samples, comparison + "_ctrl", control_samples)


# call peaks as defined in comparison table but with all inputs as background
for comparison in comparisons["comparison"].drop_duplicates():
    sub = comparisons[comparisons["comparison"] == comparison]

    treatment_samples = [s for s in prj.samples if s.name in sub[sub["comparison_side"] == 1]['sample'].tolist()]
    control_samples = [s for s in prj.samples if s.ip in ["IgG", "Input"]]

    call_peaks(treatment_samples, comparison + "_allctrls", control_samples)


# call peaks without control for all samples
for sample in [s for s in prj.samples if s.ip != "IgG"]:
    # Call peaks for all ChIPmentation samples
    call_peaks([sample], sample.name)
