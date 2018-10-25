"""
module unload R
module load R/3.5.1
R
"""

library("DiffBind")
d <- dba(sampleSheet="metadata/diffBind_design.csv")
d <- dba.count(d)
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
d <- dba.analyze(d, method=DBA_ALL_METHODS, filter=0, bSubControl=FALSE)

for (c in seq(4)){
    d.DB <- dba.report(
        d, c, 
        bNormalized=TRUE,
        th=1,
        file=paste0(
            dba.show(d, bContrasts=TRUE)[c, "Group1"],
            "_vs_",
            dba.show(d, bContrasts=TRUE)[c, "Group2"]
            ))
}