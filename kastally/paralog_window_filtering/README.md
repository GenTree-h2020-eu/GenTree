# Paralog window filtering steps

Identification of genomic windows to exclude in two steps, to run for each species:

1. run `paralog_window_filtering.R` using as input the output of the HDplot analysis ([see here](https://github.com/GenTree-h2020-eu/GenTree/blob/master/rellstab))
2. run `compute_excluded_areas.sh` to produce a .bed file with the genomic coordinates of the excluded windows.

# Script description

## paralog_window_filtering.R:

This script marks SNPs in windows with high proportions of likely erroneous calls.

Usage:

`./paralog_window_filtering.R -i input_file -o output_file`

This script takes as input a TSV table with SNPs in rows and at least 4 columns: _contig_, _pos_, _H_ and _D_, with H the heterozigosity, D the deviation to allele ratio balance and contig and pos the contig and position of the variant.

It will then output the same table with 2 or 3 additional columns:

- _filter_: TRUE (retained) or FALSE (excluded), based on H and D
- _window_score_proportion_good_snp_: window score, based on neighbour SNPs
- (optional) _window_filter_TrueIsRetained_: using the window score and if a threshold is given, TRUE (retained) or FALSE (excluded)

## compute_excluded_areas.sh:

This script takes as input the output of `paralog_window_filtering.R` and outputs a .bed file with the coordinates of the windows to exclude.

By default, this script sets the _quality score threshold_ at `0.9` (ie, if the number of bad SNPs is above 10%, the window is excluded) and the _window size_ at `250`.
