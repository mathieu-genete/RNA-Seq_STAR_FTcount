#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=2G
#SBATCH --output=LOG/%j-%a_02_Align_STAR

source fonctions.sh
source variables_RNASeq.sh

module load star/2.7.5a

fastq_array=($(tail -n +2 $slist))
l=${fastq_array[$SLURM_ARRAY_TASK_ID]}

sname=$(echo $l| cut -d ',' -f1)
l=$(echo $l| cut -d ',' -f2-)
commacount=$(echo $l | grep -o "," | wc -l)
fqReads=""
fq_fwd=""
fq_rev=""
c=0
for fq in $(echo $l | tr ',' ' ')
do
    read_fw=$((c%2))
    fqpath=$rawReadsPath"/"$fq
    test_path $fqpath
    fileexist=$?
    if [ $fileexist -eq 0 ]
    then
	test_errors 10 "impossible d'acceder au fichier $fqpath" 1
    fi

    if [ $read_fw -eq 0 ]
    then
	if [ ${#fq_fwd} -eq 0 ]
	then
	    fq_fwd=$fqpath
	else
	    fq_fwd=$fq_fwd","$fqpath
	fi
    else
	if [ ${#fq_rev} -eq 0 ]
	then
	    fq_rev=$fqpath
	else
	    fq_rev=$fq_rev","$fqpath
	fi
    fi
    c=$(($c+1))
done

fqReads=$fq_fwd" "$fq_rev

affiche_txt "Launch alignment for $sname"
STARoutfolder=$outputFolder"/STAR_results/"$sname"_"$indexName
NamePrefix=$STARoutfolder"/"$sname"_"

if [ ! -d $STARoutfolder ]
then
    mkdir -p $STARoutfolder
fi

if [ $zcatoption -eq 1 ]
then
    readFilescmd="--readFilesCommand zcat "
else
    readFilescmd=""
fi

test_path $NamePrefix"Aligned.sortedByCoord.out.bam"
test_bam=$?
if [ $test_bam -ne 1 ] || [ $force -eq 1 ]
then
    srun -J "aln_STAR_$sname" STAR --genomeDir $indexPath --runThreadN $SLURM_CPUS_PER_TASK $readFilescmd--outFileNamePrefix $NamePrefix --readFilesIn $fqReads --outSAMtype BAM SortedByCoordinate --sjdbGTFfile $gtfRef > $logdir/STAR_Align_out-$sname.log
else
    affiche_txt "PASS Alignemnt for $sname -- bam already exists"
fi
