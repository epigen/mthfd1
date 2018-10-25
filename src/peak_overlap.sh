PEAKS=(
results/chipseq_peaks/HAP1_H3K27ac_WT_IgG-background/HAP1_H3K27ac_WT_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_BRD4_WT_IgG-background/HAP1_BRD4_WT_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_MTHFD1_WT_IgG-background/HAP1_MTHFD1_WT_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_H3K27ac_MTHFD1KO_IgG-background/HAP1_H3K27ac_MTHFD1KO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_BRD4_MTHFD1KO_IgG-background/HAP1_BRD4_MTHFD1KO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_MTHFD1_MTHFD1KO_IgG-background/HAP1_MTHFD1_MTHFD1KO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_H3K27ac_DMSO_IgG-background/HAP1_H3K27ac_DMSO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_BRD4_DMSO_IgG-background/HAP1_BRD4_DMSO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_MTHFD1_DMSO_IgG-background/HAP1_MTHFD1_DMSO_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_H3K27ac_dBET6_IgG-background/HAP1_H3K27ac_dBET6_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_BRD4_dBET6_IgG-background/HAP1_BRD4_dBET6_IgG-background_homer_peaks.factor.filtered.bed
results/chipseq_peaks/HAP1_MTHFD1_dBET6_IgG-background/HAP1_MTHFD1_dBET6_IgG-background_homer_peaks.factor.filtered.bed
)

for PEAK in ${PEAKS[@]}; do
cp $PEAK 
done

# Peak overlap
PEAKS=`ls *.filtered.bed | \
sed -e 's/IgG-background_homer_peaks.//g' | \
sed -e 's/factor//g' | \
sed -e 's/histone//g' | \
sed -e 's/.filtered.bed//g' | sed -e 's/HAP1_//g'`

for file1 in ${PEAKS[@]}; do
bedtools sort -i $file1 > ${file1}.sorted.bed
mv ${file1}.sorted.bed $file1
done

echo name" "$PEAKS > results/pairwise_jaccard.txt
for file1 in ${PEAKS[@]}; do
    # make reasonable file labels
    file1_short=`echo $file1 | \
        sed -e 's/IgG-background_homer_peaks.//g' | \
        sed -e 's/factor//g' | \
        sed -e 's/histone//g' | \
        sed -e 's/.filtered.bed//g' | sed -e 's/HAP1_//g'`
    echo -n $file1_short >> results/pairwise_jaccard.txt
    for file2 in ${PEAKS[@]};
    do
        echo $file1 $file2
        # compute the jaccard stat for these two files.
        jaccard=`bedtools jaccard \
                   -a ${file1} \
                   -b ${file2} | tail -n +2 | cut -f 3`
        # report the jaccard stat for these two files
        echo -n " "$jaccard >> results/pairwise_jaccard.txt
    done
    echo >> results/pairwise_jaccard.txt
done

echo name" "$PEAKS > results/pairwise_intersect.txt
for file1 in ${PEAKS[@]}; do
    # make reasonable file labels
    file1_short=`echo $file1 | \
        sed -e 's/IgG-background_homer_peaks.//g' | \
        sed -e 's/factor//g' | \
        sed -e 's/histone//g' | \
        sed -e 's/.filtered.bed//g' | sed -e 's/HAP1_//g'`
    echo -n $file1_short >> results/pairwise_intersect.txt
    for file2 in ${PEAKS[@]};
    do
        echo $file1 $file2
        # compute the jaccard stat for these two files.
        jaccard=`bedtools jaccard \
                   -a ${file1} \
                   -b ${file2} | tail -n +2 | cut -f 4`
        # report the jaccard stat for these two files
        echo -n " "$jaccard >> results/pairwise_intersect.txt
    done
    echo >> results/pairwise_intersect.txt
done
