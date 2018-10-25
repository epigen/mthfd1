# Deeptools plots


BIGWIGS=`ls merged/HAP1_*.bigwig`
PEAKS=(
results/chipseq_peaks/HAP1_H3K27ac_WT_IgG-background/HAP1_H3K27ac_WT_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_BRD4_WT_IgG-background/HAP1_BRD4_WT_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_MTHFD1_WT_IgG-background/HAP1_MTHFD1_WT_IgG-background_homer_peaks.factor.filtered.bed
)
PEAK_NAMES=(
H3K27ac
BRD4
MTHFD1
)

declare -a SCALES
SCALES[1]='2.5 2.5 5 5 50 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5 2.5'
SCALES[2]='4 4 11 11 100 4 4 4 4 4 4 4 4 4 4 4 4'
SCALES[3]='6 6 11 11 100 7 7 7 7 7 7 7 7 7 7 7 7'


for I in `seq 3`; do
NAME=${PEAK_NAMES[$I]}
PEAK=${PEAKS[$I]}
PEAK_NUMBER=`wc -l $PEAK | sed 's/ .*//g'`

# get coverage matrix
computeMatrix reference-point \
--referencePoint center \
-p max \
-S ${BIGWIGS[@]} \
-R $PEAK \
-a 2500 -b 2500 \
-bs 50 \
-out HAP1_${NAME}_WT_IgG.mat.gz
# pretify names
gzip -d HAP1_${NAME}_WT_IgG.mat.gz
sed -i 's/.merged.sorted.bam.deeptools//g' HAP1_${NAME}_WT_IgG.mat
sed -i 's/HAP1_//g' HAP1_${NAME}_WT_IgG.mat
gzip HAP1_${NAME}_WT_IgG.mat

# plot line profile
plotProfile \
--matrixFile HAP1_${NAME}_WT_IgG.mat.gz \
--plotTitle "${NAME} peaks (n="$PEAK_NUMBER")" \
--numPlotsPerRow 4 \
--yMin 1 \
--yMax ${SCALES[$I]}
--plotType lines \
-o HAP1_${NAME}_WT_IgG.lineplot.svg
sed -i 's/genes//g' HAP1_${NAME}_WT_IgG.lineplot.svg

# plot heatmaps
plotHeatmap \
--matrixFile HAP1_${NAME}_WT_IgG.mat.gz \
--plotTitle "${NAME} peaks (n="$PEAK_NUMBER")" \
--colorMap RdBu_r \
--heatmapHeight 12 \
--yMin 1 \
--yMax ${SCALES[$I]}
--zMin 1 \
--zMax ${SCALES[$I]}
--dpi 300 \
-o HAP1_${NAME}_WT_IgG.heatmap.svg
sed -i 's/genes//g' HAP1_${NAME}_WT_IgG.heatmap.svg
sed -i 's/gene distance/distance/g' HAP1_${NAME}_WT_IgG.heatmap.svg

done
