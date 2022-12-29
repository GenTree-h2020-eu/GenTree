# GenTree
Repository to share and document GenTree H2020 project bioinformatic pipelines. The repository contains code used for genomic analysis of 7 GenTree project species (Betula pendula, Fagus sylvatica, Picea abies, Pinus pinaster, Pinus sylvestris, Populus nigra and Quercus petraea).

Raw short read data:

Picea abies: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602465

Betula pendula: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602466

Quercus petraea: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602467

Pinus pinaster: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602468

Fagus sylvatica: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602470

Populus nigra: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602471

Pinus sylvestris: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602473

## Bioinformatic pipeline
1.	Mapping, https://github.com/GenTree-h2020-eu/GenTree/blob/master/Alignment_commands.txt
2.	SNP calling, https://github.com/GenTree-h2020-eu/GenTree/blob/master/SNP_calling_commands.txt
3.	Initial quality control with scikit-allele. Jupyter notebooks available at https://github.com/GenTree-h2020-eu/GenTree/tree/master/cervantesarango/JupyterNotebooks. Based on code available at: http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html and http://alimanfoo.github.io/2015/09/28/fast-pca.html
4.	Formal variant filtering via GATK
    - From raw SNPs to v2: https://github.com/GenTree-h2020-eu/GenTree/tree/master/cervantesarango/GATK_rawSNPs_to_v2

5.	Filtering of paralog regions
6. Genotype filtering
   -  from v3.1 to v5.1: https://github.com/GenTree-h2020-eu/GenTree/blob/master/GenTree_7species_GenotypeSNPfiltering_v5.1.txt
   -  from v3.1 to v5.3: https://github.com/GenTree-h2020-eu/GenTree/blob/master/GenTree_7species_GenotypeSNPfiltering_v5.3.txt
   -  from v5.3 to v5.3.1 https://github.com/GenTree-h2020-eu/GenTree/blob/master/GenTree_7species_GenotypeSNPfiltering_v5.3.1.txt

