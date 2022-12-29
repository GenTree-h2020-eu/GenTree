#!/bin/bash 
#SBATCH --job Fagus_GM_select_snps
#SBATCH -o /scratch/project_2001229/sandra/Fagus/Fagus_select_snps_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Fagus/Fagus_select_snps_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0

srun gatk SelectVariants --variant Fagus_TotalVariants_RENAMED.vcf.gz -select-type SNP --output Fsylvatica_GM_raw_SNPs.vcf.gz


