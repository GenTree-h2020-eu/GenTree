#!/bin/bash 
#SBATCH --job Fagus_GM_reheader
#SBATCH -o /scratch/project_2001229/sandra/Fagus/reheader.out.txt
#SBATCH -e /scratch/project_2001229/sandra/Fagus/reheader.err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=test
#SBATCH --time=00:10:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load biokit

bcftools reheader Fagus_TotalVariants.vcf.gz -s Fsylvatica_reheader.txt > Fagus_TotalVariants_RENAMED.vcf.gz
