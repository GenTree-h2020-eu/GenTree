#!/bin/bash 
#SBATCH --job Fagus_exclude_samples
#SBATCH -o /scratch/project_2001229/sandra/Fagus/Fagus_exclude_samples_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Fagus/Fagus_exclude_samples_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0

srun gatk SelectVariants --variant Fsylvatica_GM_raw_SNPs.vcf.gz --output Fsylvatica_GM_samples_KEPT.vcf.gz --exclude-sample-name samples_to_exclude_Fsylvatica.list


