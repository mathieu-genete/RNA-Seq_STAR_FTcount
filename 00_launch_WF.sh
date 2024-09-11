#!/bin/bash

#SBATCH --partition=SLURM_PARTITION
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --output=LOG/%j-_00_Launch_WF

partition="SLURM_PARTITION"
outputFolder="results_RNASeq" #chemin vers le dossier de sortie pour vos résultats
                                                          
slist="./samples_list" # chemin vers le fichier texte contenant la liste de vos samples et des fichiers fastq
#   le fichier doit être au format suivant:
#   /chemin/vers/les/fichiers/fastq
#   nom_sample1,fichier_sample1_R1.fastq,fichier_sample1_R2.fasq
#   nom_sample2,fichier_sample2_R1.fastq,fichier_sample2_R2.fasq
#   ....
#   ....

readslength=51 # taille de vos reads (vous trouverez l'information dans le controle qualité avec fsatqc et multiqc)

fastaRef="None" # chemin vers le fichier fasta de votre référence
gtfRef="ref/gtf/Homo_sapiens.GRCh38.101.gtf" # chemin vers le fichier d'annoation (gtf) de votre référence

#Variables pour la création de l'index
indexPath="ref/hg38_star" #mettre ici le chemin vers le genome indexé si il existe... ça évite de refaire l'indexation du génome à chaque fois. Après vous pouvez aller récupérer le génome index sur le cluster (/shared/bank/homo_sapiens/version_genome/star-2.7.2b) en choisissant la bonne version
indexName="hg38_star" # nom que vous souhaitez donner à votre référence indexée

#IndexingCPUS=20 #nombre de CPU pour l'indexation du génome
#IndexMEMSTAR=2GB #taille de la mémoire à utiliser (j'ai laissé la valeur par défaut, il faudra peut-être l'augmenter)

#variables pour le mapping
AlignParallelJobs=3 #Nombre d'alignements STAR lancés en parralèle => exemple si vous mettez AlignParallelJobs=3 et AlignCPUS=10 vous allez lancer en paralèle 3 processus STAR qui concommeront 10 cpus chaqun. Soit au total 30 CPUs. Même calcul pour la mémoire
AlignCPUS=10 #Nombre de CPUs par alignement
AlignMEM=2GB #mémoire alouée par aligmenent
zcatoption=1

#variables pour le comptage
featureCountsParallelJobs=3 #même chose que pour le mapping mais pour featureCounts
featureCountsCPUS=10 #idem mapping mais pour featureCounts
fcFeatureType="gene" # option -t featureCounts (gene,exon,...)
fcAttributeType="gene_id" # option -g featureCounts (gene_id,gene_name,...)

if [ $indexPath ]
then
    indexedOUT=$indexPath
else
    indexedOUT=$outputFolder"/indexed_"$indexName
    indexPath=$indexedOUT
fi

logdir=$outputFolder"/logs"
mkdir -p $logdir

#Export variables dans un autre fichier
rawReadsPath=$(head -n1 $slist)
variable_file="variables_RNASeq.sh"
echo "slist=\"$slist\"">$variable_file
echo "outputFolder=\"$outputFolder\"">>$variable_file
echo "rawReadsPath=\"$rawReadsPath\"">>$variable_file
echo "zcatoption=$zcatoption">>$variable_file
echo "indexPath=\"$indexPath\"">>$variable_file
echo "gtfRef=\"$gtfRef\"">>$variable_file
echo "logdir=\"$logdir\"">>$variable_file
echo "indexName=\"$indexName\"">>$variable_file
echo "fcFeatureType=\"$fcFeatureType\"">>$variable_file
echo "fcAttributeType=\"$fcAttributeType\"">>$variable_file

#réalise l'index du génome
job_ID_Index=$(sbatch --wait --partition=$partition --parsable 01_index_REF.sh $indexedOUT $fastaRef $gtfRef $readslength $logdir)
echo "Job_ID Index STAR: $job_ID_Index"

#Alignement avec STAR
fastq_array=($(tail -n +2 $slist))
nbr_samples=$((${#fastq_array[@]}-1))
echo $nbr_samples "0-$nbr_samples%$AlignParallelJobs"

job_ID_Index=$(sbatch --wait --partition=$partition --array="0-$nbr_samples%$AlignParallelJobs" --parsable 02_Align_STAR.sh)
echo "Job_ID Align STAR: $job_ID_Index"

#Comptage avec FeatureCounts
bam_list=($(find $outputFolder"/STAR_results/" -name "*.bam"))
nbr_bams=$((${#bam_list[@]}-1))
echo $nbr_bams "0-$nbr_bams%$featureCountsParallelJobs"

job_ID_Index=$(sbatch --wait --partition=$partition --array="0-$nbr_bams%$featureCountsParallelJobs" --parsable 03_count_featureCounts.sh)
echo "Job_ID FeatureCounts: $job_ID_Index"
