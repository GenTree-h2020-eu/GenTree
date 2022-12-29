#!/bin/bash
#SBATCH --job Bpendula_GM_GATK_exclude_samples
#SBATCH -o /scratch/project_2001229/sandra/Betula/Betula_recalibration_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Betula/Betula_recalibration_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000


module load gatk/4.1.4.0


gatk SelectVariants --variant Bpendula_FLAGGED.vcf.gz --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC --output Bpendula_GM_Oulu_filtered_v2.vcf.gz




