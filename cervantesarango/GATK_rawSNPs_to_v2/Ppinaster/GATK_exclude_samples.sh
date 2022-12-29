#!/bin/bash 
#SBATCH --job Ppinaster_GM_GATK_exclude_samples
#SBATCH -o /scratch/project_2001229/sandra/Ppinaster/exclude_samples_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Ppinaster/exclude_samples_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16000


module load gatk/4.1.4.0

srun gatk SelectVariants --variant Ppinaster_raw_SNPs.vcf.gz  --output Ppinaster_TRIAL.vcf.gz --exclude-sample-name Ppinaster_exclude_SORTED.list


