#!/bin/bash
#SBATCH --job Populus_GM_GATK_recalibration
#SBATCH -o /scratch/project_2001229/sandra/Populus/Pnigra_recalibration_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Populus/Pnigra_recalibration_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=test
#SBATCH --time=00:15:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000



module load gatk/4.1.4.0


gatk SelectVariants --variant Pnigra_FLAGGED.vcf.gz --exclude-filtered true --exclude-non-variants true --restrict-alleles-to BIALLELIC --output Pnigra_GM_Oulu_filtered_v2.vcf.gz




