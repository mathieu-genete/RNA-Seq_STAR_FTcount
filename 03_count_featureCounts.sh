#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --output=LOG/%j-%a_03_FeatureCounts

source fonctions.sh
source variables_RNASeq.sh

module load samtools/1.10
module load subread/2.0.1

affiche_txt "=== STEP 3 : featuresCounts ==="

fcFolder=$outputFolder"/featureCounts/"
mkdir -p $fcFolder

bam_list=($(find $outputFolder"/STAR_results/" -name "*.bam"))

bam=${bam_list[$SLURM_ARRAY_TASK_ID]} 

sname=$(basename -s "_Aligned.sortedByCoord.out.bam" $bam)

outcounts=$fcFolder/$sname-counts
affiche_txt "Launch featureCounts for $sname"

test_path $outcounts
test_counts=$?
if [ $test_counts -ne 1 ] || [ $force -eq 1 ]
then
    srun -J "samtools index" samtools index $bam
    srun -J "" featureCounts -T $SLURM_CPUS_PER_TASK -O -g $fcAttributeType -t $fcFeatureType -a $gtfRef -o $outcounts $bam > $logdir/Fcounts_out-$sname.log
else
    affiche_txt "PASS counts for $sname -- featureCounts file already exists"
fi
