#!/usr/bin/env Rscript

# Window-based paralog filtering for GenTree
#
# Author Chedly Kastally <ckastall@gmail.com>
# Copyright (C) 2020 Chedly Kastally <ckastall@gmail.com>
# Modified On 2020-02-26 19:07
# Created  2020-02-18 16:30
# Licence: GPL3

# use --help for usage.

# function ----

eval_Windows <- function(dat_single_contig,
                         # window_overlap = window_overlap_default,
                         w_size = 100) {

    # Overlapping is not implemented at the moment

    # This function look at the proportion of TRUE values in column "filter" at a
    # distance of windiow_size / 2 up and downstream for each "pos".

    # checks
    stopifnot(class(dat_single_contig) == "data.frame")

    expected_columns <- c("pos", "filter")
    stopifnot(expected_columns %in% colnames(dat_single_contig))

    result <- lapply(dat_single_contig$pos, function(x) {

               xi <- x - w_size/2
               xj <- x + w_size/2

               snps_in_window <-
                   which(dat_single_contig$pos >= xi &
                         dat_single_contig$pos <= xj)

               nSNPs_in_window <- 
                   length(snps_in_window)

               nFiltered_in_window <-
                   sum(dat_single_contig[snps_in_window, "filter"])

               nFiltered_in_window / nSNPs_in_window

             })

    unlist(result)
}


# Optparse ----

library("optparse")
option_list <-
    list(
         make_option(c("-i", "--input_file"),
         type    = "character",
         default = NULL,
         help    = "input_file [Default: %default]"),
         make_option(c("-o", "--output_file"),
         type    = "character",
         default = "window_based_filtering.tsv",
         help    = "output_file [Default: %default]"),
         make_option(c("-t", "--window_filtering_threshold"),
         type    = "numeric",
         default = NULL,
         metavar = "numeric",
         help    = "Window threshold: below that value, all positions in a given window are filtered out [Default: %default]"),
         make_option(c("-w", "--window_size"),
         type    = "integer",
         default = 100,
         metavar = "integer",
         help    = "Window size [Default: %default]")
         )

usage <- "%prog -i input_file -o output_file
This script takes as input a TSV table with at least 3 columns:
H, D and pos
 
It will then output the same table with 2 or 3 additional columns:
    - filter: TRUE (retained) or FALSE (excluded), based on H and D
    - window_score_proportion_good_snp: window score, based on neighbour SNPs
    - (optional) window_filter_TrueIsRetained: using the window score
      and if a threshold is given, TRUE (retained) or FALSE (excluded)
"


parser         <- OptionParser(usage = usage,
                   option_list = option_list)

opt            <- parse_args(parser)

input_file                 <- opt$input_file
output_file                <- opt$output_file
window_size                <- opt$window_size
window_filtering_threshold <- opt$window_filtering_threshold

# if no filtering threshold is specified, the logical window_filter won't be
# printed in the output column. It's actually a bit redundant to output the
# score and the window_filter, and since the score is more informative, I'm leaving
# the logical to be optional

# cat(sprintf("Window size is: %s\nWindow filtering threshold is: %s\n",
#             window_size, window_filtering_threshold))

stopifnot(is.integer(window_size) && window_size > 2)
if (!is.null(window_filtering_threshold)) {
    stopifnot(is.numeric(window_filtering_threshold) &&
              window_filtering_threshold > 0 &&
              window_filtering_threshold < 1)
}

required_files <- c(input_file)

# check whether the output file is in a subdir, and that the directory exists
if (!(dir.exists(dirname(output_file)))) {
    stop("Output file's directory does not exist. You have to create it first, or give a valid path.")
}

for (i in required_files) {
    if(!file.exists(i)) {
    stop(sprintf("ERROR: file %s does not exist; check the path of the input file(s).\n",
                 i))
    }
}

# library ----

library("dplyr")

# setup ----

H_upper_limit <- .6
D_upper_limit <- 20
D_lower_limit <- -20

# main ----

# Read the data
HD <- read.table(input_file, sep = "\t", header = T, stringsAsFactors = F)

# Some checks
required_colnames <- c("H", "D", "pos", "contig")
stopifnot(required_colnames %in% colnames(HD))

# Compute the proportion of reference alleles

HD <- HD %>%
    na.omit(.)

stopifnot(!any(is.na(HD)))

# By default, everything is retained, and retain = TRUE
HD$filter<-TRUE

# We tag FALSE all those with bad properties based on H and D
# TODO: check the effects of H and D parameters?
HD$filter[HD$H > H_upper_limit |
          HD$D > D_upper_limit | 
          HD$D < D_lower_limit] <- FALSE

total_number_contigs <- length(unique(HD$contig))

all_scores <-
    lapply(1:total_number_contigs, function(x) {

               cat(sprintf("Proceeding with contig number: %s on %s\n",
                           x, total_number_contigs))

               contigID <- unique(HD$contig)[x]

               # Computation is done here: return the proportion of GOOD SNPs
               # in a window centered on each pos considered
                eval_Windows(HD[HD$contig == contigID, ],
                             w_size = window_size)
       })

all_scores_2 <- unlist(all_scores)

stopifnot(nrow(HD) == length(all_scores_2))

HD$window_score_proportion_good_snp <- all_scores_2

if (!is.null(window_filtering_threshold)) {
    HD$window_filter_TrueIsRetained     <- all_scores_2 > window_filtering_threshold & HD$filter

    # TODO: if a window filterering threshold is used, output also a summary of
    # the results obtained!

}

write.table(HD, file = output_file, sep = "\t", row.names = F, quote = F)
