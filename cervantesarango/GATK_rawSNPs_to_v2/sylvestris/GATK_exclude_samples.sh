#!/bin/bash 
#SBATCH --job Psylvestris_GM_GATK_exclude_samples
#SBATCH -o /scratch/project_2001229/sandra/Psylvestris/exclude_samples_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Psylvestris/exclude_samples_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=02:30:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0

srun gatk SelectVariants --variant Psylvestris_GM_raw_snps.vcf.gz  --output Psylvestris_GM_samples_kept.vcf.gz --exclude-sample-name samples_to_exclude_Psylvestris_GM.list


