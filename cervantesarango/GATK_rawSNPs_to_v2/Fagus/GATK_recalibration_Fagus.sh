#!/bin/bash
#SBATCH --job Fagus_GM_GATK_recalibration
#SBATCH -o /scratch/project_2001229/sandra/Fagus/Fagus_recalibration_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Picea/Fagus_recalibration_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0


gatk SelectVariants --variant Fsylvatica_FLAGGED.vcf.gz --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC --output Fsylvatica_GM_Oulu_filtered_v2._RENAMED.vcf.gz




