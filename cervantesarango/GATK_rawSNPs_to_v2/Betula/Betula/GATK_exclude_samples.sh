#!/bin/bash 
#SBATCH --job Bpendula_GM_GATK_exclude_samples
#SBATCH -o /scratch/project_2001229/sandra/Betula_exclude_samples_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Betula_exclude_samples_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000

module load gatk/4.1.4.0

srun gatk SelectVariants --variant Bpendula_GM_raw_snps.vcf.gz  --output Bpendula_GM_samples_kept.vcf.gz --exclude-sample-name Bpendula_exclude.list


