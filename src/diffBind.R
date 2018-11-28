"""
module unload R
module load R/3.5.1
R
"""

library("DiffBind")

outputDir = "results/differential_analysis_ChIP-seq"

# Create DBA object
d <- dba(sampleSheet="metadata/diffBind_design.csv")
# Get count data
d <- dba.count(d)
# Add contrasts
d <- dba.contrast(
    d,
    dba.mask(d, DBA_CONDITION, 'MTHFD1KO') + dba.mask(d, DBA_FACTOR, 'BRD4') == 2,
    dba.mask(d, DBA_CONDITION, 'WT') + dba.mask(d, DBA_FACTOR, 'BRD4') == 2,
    "BRD4_MTHFD1KO",
    "BRD4_WT",
    )
d <- dba.contrast(
    d,
    dba.mask(d, DBA_CONDITION, 'MTHFD1KO') + dba.mask(d, DBA_FACTOR, 'MTHFD1') == 2,
    dba.mask(d, DBA_CONDITION, 'WT') + dba.mask(d, DBA_FACTOR, 'MTHFD1') == 2,
    "MTHFD1_MTHFD1KO",
    "MTHFD1_WT",
    )
d <- dba.contrast(
    d,
    dba.mask(d, DBA_CONDITION, 'dBET') + dba.mask(d, DBA_FACTOR, 'BRD4') == 2,
    dba.mask(d, DBA_CONDITION, 'DMSO') + dba.mask(d, DBA_FACTOR, 'BRD4') == 2,
    "BRD4_dBET",
    "BRD4_DMSO",
    )
d <- dba.contrast(
    d,
    dba.mask(d, DBA_CONDITION, 'dBET') + dba.mask(d, DBA_FACTOR, 'MTHFD1') == 2,
    dba.mask(d, DBA_CONDITION, 'DMSO') + dba.mask(d, DBA_FACTOR, 'MTHFD1') == 2,
    "MTHFD1_dBET",
    "MTHFD1_DMSO",
    )

# Run comparison of contrasts
d <- dba.analyze(d, method=DBA_ALL_METHODS, filter=0, bSubControl=FALSE)

# Save
for (c in seq(length(d$contrasts))){
    s <- paste0(
            dba.show(d, bContrasts=TRUE)[c, "Group1"],
            "_vs_",
            dba.show(d, bContrasts=TRUE)[c, "Group2"]
            )
    d.DB <- dba.report(
        d, c, 
        bNormalized=TRUE,
        th=1,
        file=s)
    # save in 
    file.rename(
        paste0("DBA_", s, ".csv"),
        paste0(outputDir, "/DBA_", s, ".csv"))
}
