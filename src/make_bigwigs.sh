PROJECT_NAME=mthfd1
cd ~/projects/$PROJECT_NAME


# For each sample
readarray SAMPLES < <(find data -maxdepth 1 ! -name external ! -name data -type d)

for SAMPLE in ${SAMPLES[@]}; do
SAMPLE=${SAMPLE/data\//}
echo $SAMPLE
sbatch -J $SAMPLE.deeptools.bigwig --time "07:00:00" -p shortq -c 4 --mem 16000 \
-o /scratch/lab_bock/shared/projects/${PROJECT_NAME}/data/${SAMPLE}/coverage/${SAMPLE}.trimmed.bowtie2.deeptools.log \
--wrap "bamCoverage \
--bam /scratch/lab_bock/shared/projects/${PROJECT_NAME}/data/${SAMPLE}/mapped/${SAMPLE}.trimmed.bowtie2.bam \
-o /scratch/lab_bock/shared/projects/${PROJECT_NAME}/data/${SAMPLE}/coverage/${SAMPLE}.trimmed.bowtie2.deeptools.bigwig \
-p max \
--binSize 10 \
--normalizeUsing RPGC \
--effectiveGenomeSize 3000000000 \
--extendReads 175"
# sleep 360
done

# For merged samples
readarray SAMPLE_GROUPS < <(find merged -maxdepth 1 ! -name external ! -name merged -type f -name "*.merged.sorted.bam")

for SAMPLE in ${SAMPLE_GROUPS[@]}; do
SAMPLE=${SAMPLE/merged\//}
echo $SAMPLE
sbatch -J $SAMPLE.deeptools.bigwig --time "07:00:00" -p shortq -c 4 --mem 16000 \
-o /scratch/lab_bock/shared/projects/${PROJECT_NAME}/merged/${SAMPLE}.deeptools.bigwig.log \
--wrap "bamCoverage \
--bam /scratch/lab_bock/shared/projects/${PROJECT_NAME}/merged/${SAMPLE} \
-o /scratch/lab_bock/shared/projects/${PROJECT_NAME}/merged/${SAMPLE}.deeptools.bigwig \
-p max \
--binSize 10 \
--normalizeUsing RPGC \
--effectiveGenomeSize 3000000000 \
--extendReads 175"
# sleep 360
done
