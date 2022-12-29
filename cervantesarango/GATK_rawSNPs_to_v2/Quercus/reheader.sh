#!/bin/bash 
#SBATCH --job Quercus_GM_reheader
#SBATCH -o /scratch/project_2001229/sandra/Quercus/Quercus_reheader.out.txt
#SBATCH -e /scratch/project_2001229/sandra/Quercus/Quercus_reheader.err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=test
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000

module load biokit

bcftools reheader Quercus_TotalVariants.vcf.gz -s reheader_samples_Quercus_TotalVariants.txt > Quercus_TotalVariants_RENAMED.vcf.gz
