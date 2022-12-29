#!/bin/bash 
#SBATCH --job Quercus_GM_select_snps
#SBATCH -o /scratch/project_2001229/sandra/Quercus/Quercus_select_snps_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Quercus/Quercus_select_snps_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=00:50:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0

srun gatk SelectVariants --variant Quercus_TotalVariants_RENAMED.vcf.gz -select-type SNP --output Qpetraea_GM_raw_SNPs.vcf.gz


