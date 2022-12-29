#!/bin/bash
#SBATCH --job pinaster_GM_hard_filters
#SBATCH -o /scratch/project_2001229/sandra/Ppinaster/hard_filters_out.txt
#SBATCH -e /scratch/project_2001229/sandra/Ppinaster/hard_filters_err.txt
#SBATCH --account=project_2001229
#SBATCH --partition=small
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4000



module load gatk/4.1.4.0

srun gatk VariantFiltration --variant Ppinaster_GM_samples_kept.vcf.gz  -filter "QD < 0.25" --filter-name "QD0.25" -filter "MQ < 30.0" --filter-name "MQ30"  -filter "QUAL < 20.0" --filter-name "QUAL20" -filter "SOR > 3.0" --filter-name "SOR3" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -filter "AN < 582.0" --filter-name "AN582" --output Ppinaster_FLAGGED.vcf.gz





