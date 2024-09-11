#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --output=LOG/%j_01_Index_REF

source fonctions.sh

module load star/2.7.5a

indexedOUT=$1
fastaRef=$2
gtfRef=$3
readslength=$4
logdir=$5

test_path $indexedOUT"/SAindex"
test_SAindex=$?

if [ $test_SAindex -eq 1 ]
then
    indexPath=$indexedOUT
fi

if [ ! $indexPath ]
then
    affiche_txt "=== STEP 1 : index reference genome ==="
    
    #Create index output directory if not exists
    if [ ! -d $indexedOUT ]
    then
        affiche_txt "create output directory $indexedOUT"
        mkdir -p $indexedOUT
    fi
    
    #launch STAR command to index the genome
    affiche_txt "Launch index generation"
    
    srun -J "Index_STAR" STAR --runThreadN $SLURM_CPUS_PER_TASK --runMode genomeGenerate --genomeDir $indexedOUT --genomeFastaFiles $fastaRef --sjdbGTFfile $gtfRef --sjdbOverhang $(($readslength-1)) --genomeSAindexNbases 12 2> $logdir/STAR_index_errlog.log | tee $logdir/STAR_index_stdout.log

    indexPath=$indexedOUT
fi

affiche_txt "END Genome Index"
