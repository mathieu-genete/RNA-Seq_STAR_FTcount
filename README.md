# README

## Description

Ce script Bash est conçu pour automatiser le processus de traitement des données RNA-Seq en utilisant SLURM pour la gestion des tâches. Il inclut des étapes pour l'indexation du génome, l'alignement des reads avec STAR, et le comptage des features avec FeatureCounts.

## Prérequis

- SLURM
- STAR
- FeatureCounts
- Un fichier de référence FASTA
- Un fichier d'annotation GTF

## Utilisation

### Configuration

1. **Partition SLURM**: Définissez la partition SLURM à utiliser.
2. **Dossier de sortie**: Spécifiez le chemin vers le dossier de sortie pour les résultats.
3. **Liste des échantillons**: Fournissez le chemin vers un fichier texte contenant la liste de vos échantillons et des fichiers FASTQ. Le fichier doit être au format suivant :
    ```
    /chemin/vers/les/fichiers/fastq
    nom_sample1,fichier_sample1_R1.fastq,fichier_sample1_R2.fastq
    nom_sample2,fichier_sample2_R1.fastq,fichier_sample2_R2.fastq
    ...
    ```
4. **Taille des reads**: Indiquez la taille de vos reads.
5. **Fichier de référence FASTA**: Spécifiez le chemin vers le fichier FASTA de votre référence.
6. **Fichier d'annotation GTF**: Spécifiez le chemin vers le fichier d'annotation GTF de votre référence.

### Variables

- `partition`: Partition SLURM à utiliser.
- `outputFolder`: Chemin vers le dossier de sortie pour les résultats.
- `slist`: Chemin vers le fichier texte contenant la liste des échantillons et des fichiers FASTQ.
- `readslength`: Taille des reads.
- `fastaRef`: Chemin vers le fichier FASTA de votre référence.
- `gtfRef`: Chemin vers le fichier d'annotation GTF de votre référence.
- `indexPath`: Chemin vers le génome indexé (si disponible).
- `indexName`: Nom de la référence indexée.
- `AlignParallelJobs`: Nombre d'alignements STAR lancés en parallèle.
- `AlignCPUS`: Nombre de CPUs par alignement.
- `AlignMEM`: Mémoire allouée par alignement.
- `zcatoption`: Option pour zcat.
- `featureCountsParallelJobs`: Nombre de comptages FeatureCounts lancés en parallèle.
- `featureCountsCPUS`: Nombre de CPUs par comptage.
- `fcFeatureType`: Type de feature pour FeatureCounts.
- `fcAttributeType`: Type d'attribut pour FeatureCounts.

## Notes

- Assurez-vous que tous les chemins et fichiers nécessaires sont correctement définis avant d'exécuter le script.
- Les paramètres de mémoire et de CPU peuvent être ajustés en fonction des ressources disponibles et des besoins spécifiques de votre analyse.
