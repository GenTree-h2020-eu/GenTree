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

srun gatk VariantFiltration --variant Bpendula_GM_samples_kept.vcf.gz  -filter "MQ < 30.0" --filter-name "MQ30"  -filter "QD < 0.25" --filter-name "QD0.25" -filter "QUAL < 20.0" --filter-name "QUAL20" -filter "SOR > 3.0" --filter-name "SOR3" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -filter "AN < 507.0" --filter-name "AN507" --output Bpendula_FLAGGED.vcf.gz





