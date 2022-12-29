#!/bin/bash
#SBATCH --job recalibrate_samples
#SBATCH -o /scratch/project_2001229/sandra/Psylvestris/recalibrate_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Psylvestris/recalibration_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0


gatk SelectVariants --variant Psylvestris_GM_FLAGGED_no_GRS.vcf.gz --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC --output Psylvestris_GM_Oulu_filtered_v2.vcf.gz




