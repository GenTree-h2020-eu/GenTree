#how to run HDplot in Python (and Bash and R) for the GenTree data
#christian rellstab 2020
#example on BP samples

#we need the following py-scripts
#HDplot_python.p645.py (plot functions are not working, but we don't care)
#vcf_to_depth_p645.py

gunzip Bpendula_GM_Oulu_filtered_v2.vcf.gz

#load modules (this is cluster specific)
module load gdc python/2.7.6
source /cluster/project/gdc/shared/tools/miniconda/bin/activate statsmodels
module load gcc/4.8.2 gdc vcftools/0.1.15 perl/5.18.4 zlib/1.2.8 

#run HDplot
xvfb-run python ./HDplot_python.p645.py Bpendula_GM_Oulu_filtered_v2.vcf > hdplot.log 2>hdplot.err 


#count missing data
vcftools --vcf Bpendula_GM_Oulu_filtered_v2.vcf --out Bpendula_GM_Oulu_filtered_v2.vcf.miss --missing-site

#run R-script "HDplot_analyze.R" to create the final data file and diagnostic plots
#the file "HDplot_results_column_names.txt" delivers details on the different columns in the final HDplot result file created in the R script

