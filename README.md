# GenTree
Repository to share and document GenTree H2020 project bioinformatic pipelines. The repository contains code used for genomic analysis of 7 GenTree project species (_Betula pendula, Fagus sylvatica, Picea abies, Pinus pinaster, Pinus sylvestris, Populus nigra_ and _Quercus petraea_).

## Raw short read data

_Picea abies_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602465

_Betula pendula_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602466

_Quercus petraea_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602467

_Pinus pinaster_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602468

_Fagus sylvatica_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602470

_Populus nigra_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602471

_Pinus sylvestris_: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA602473

## Bioinformatic pipeline
1.	Mapping, https://github.com/GenTree-h2020-eu/GenTree/blob/master/Alignment_commands.txt
2.	SNP calling, https://github.com/GenTree-h2020-eu/GenTree/blob/master/SNP_calling_commands.txt
3.	Initial quality control with scikit-allele. Jupyter notebooks available at https://github.com/GenTree-h2020-eu/GenTree/tree/master/cervantesarango/JupyterNotebooks. Based on code available at: http://alimanfoo.github.io/2016/06/10/scikit-allel-tour.html and http://alimanfoo.github.io/2015/09/28/fast-pca.html
4.	Formal variant filtering via GATK
    - From raw SNPs to v2: https://github.com/GenTree-h2020-eu/GenTree/tree/master/cervantesarango/GATK_rawSNPs_to_v2
5.	Filtering of paralog regions
    - https://github.com/GenTree-h2020-eu/GenTree/tree/master/rellstab
    - https://github.com/GenTree-h2020-eu/GenTree/tree/master/kastally/paralog_window_filtering
    - from v2 to v3.1 https://github.com/GenTree-h2020-eu/GenTree/tree/master/cervantesarango/ParalogRemoval
6. Genotype filtering
   -  from v3.1 to v5.3: https://github.com/GenTree-h2020-eu/GenTree/blob/master/GenTree_7species_GenotypeSNPfiltering_v5.3.txt
   -  from v5.3 to v5.3.1 https://github.com/GenTree-h2020-eu/GenTree/blob/master/GenTree_7species_GenotypeSNPfiltering_v5.3.1.txt

## Dataset descriptions
As different types of analyses require SNP sets curated and filtered based on their specific needs, we provide four versions of SNP sets as vcf files available at (link TBA). Note that since several species were observed to have experienced various levels and extents of hybridization, some admixed populations were included/excluded based on the purpose of a particular analysis. For example, in the analyses of population structure, we aimed to identify the potential admixed individuals and populations. In contrast, admixed individuals and populations can have a disproportionate effect on the measures of diversity and on the site frequency spectrum and were thus excluded from the corresponding analyses. 

- v.5.3 Known other species and clear hybrids removed, samples and SNPs with poor coverage or other low quality removed, organelle contigs removed (described in SNP filtering), _Populus nigra_ clones and cultivars removed.
- v.5.3.1 Master dataset derived from v.5.3, without samples with incorrect taxon assignment as indicated by genetic analysis for _F. sylvatica_ and _Q. petraea_, identical to v.5.3 for the other five species.
- v.5.3.2 Derived from v.5.3.1, excludes _Picea abies_ populations RU_PA_19, RU_PA_20 and RU_PO.
- v.6.3.1 Only includes four-fold degenerate sites, intron and intergenic sites, SNPs in high LD (1 kb windows, r > 0.5, PLINK v.1.90b4.9) were excluded, derived from v.5.3.1.

See https://github.com/GenTree-h2020-eu/GenTree/tree/master/ID_tracks_of_samples to follow the name changes and inclusion/exclusion of samples through the data processign and filtering.
