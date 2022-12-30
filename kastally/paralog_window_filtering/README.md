# Description of files

1. paralog_window_filtering.R:

This script marks SNPs in windows with high proportions of likely erroneous calls.

Usage: ./paralog_window_filtering.R -i input_file -o output_file

This script takes as input a TSV table with SNPs in rows and at least 4 columns: contig, pos, H and D, with H the heterozigosity, D the deviation to allele ratio balance and contig and pos the contig and position of the variant.

It will then output the same table with 2 or 3 additional columns:

- filter: TRUE (retained) or FALSE (excluded), based on H and D
- window_score_proportion_good_snp: window score, based on neighbour SNPs
- (optional) window_filter_TrueIsRetained: using the window score and if a threshold is given, TRUE (retained) or FALSE (excluded)

2. compute_excluded_areas.sh:

Takes as input the output of paralog_window_filtering.R and output a .bed file with the coordinates of the windows to exclude. This assumes a quality score threshold of 0.9 (ie, if the number of bad SNPs is above 10%, the window is excluded) and a window size of 250.
