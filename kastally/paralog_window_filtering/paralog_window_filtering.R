#!/usr/bin/env Rscript

# Window-based paralog filtering for GenTree
#
# Author Chedly Kastally <ckastall@gmail.com>
# Copyright (C) 2020 Chedly Kastally <ckastall@gmail.com>
# Modified On 2020-04-14 18:05
# Created  2020-02-18 16:30
# Licence: GPL3

# function ----

eval_Windows <- function(dat_single_contig,
                         w_size = 100) {

    if (length(dat_single_contig$pos) == 1) {

        ifelse(dat_single_contig$filter != "excluded", 1, 0)

    } else if (length(dat_single_contig$pos) > 5000) {
        # try to split dataframe if it's too big: matrix size is pos * pos

        # Split done at positions where distances between 2 consecutive positions are > (10 * windowSize)/2
        # Pos needs to be sorted here; I do that at the begining of the script.
        consecutive_distances <-
            c(0, dat_single_contig$pos[-1] - dat_single_contig$pos[-length(dat_single_contig$pos)])

        factor_div <- as.factor(cumsum(consecutive_distances > (2 * w_size / 2)))

        stopifnot(length(factor_div) == nrow(dat_single_contig))

        sub_list <- split(dat_single_contig, factor_div)

        lapply(seq_along(sub_list),
               function(y) eval_Windows(sub_list[[y]], w_size = window_size))

    } else {

        xi <- dat_single_contig$pos - w_size/2
        xj <- dat_single_contig$pos + w_size/2

        win_matrix <- matrix(c(xi, xj), byrow = F, ncol = 2)

        nSNPs_in_windows <- apply(win_matrix, 1, function(x) {
                                      dat_single_contig$pos >= x[1] &
                                          dat_single_contig$pos <= x[2]
                  })

        total    <- colSums(nSNPs_in_windows &
                            dat_single_contig$filter != "unknown")

        retained <- colSums(nSNPs_in_windows &
                            dat_single_contig$filter == "retained")

        ifelse(total == 0, 1, retained/total)

}}

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
         make_option(c("-u", "--unknown_as_good"),
         action  = "store_true",
         type    = "logical",
         default = FALSE,
         metavar = "logical",
         help    = "Keep unknown/unscores positions by HD as good positions [Default: %default]"),
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
unknown_as_good            <- opt$unknown_as_good

# if no filtering threshold is specified, the logical window_filter won't be
# printed in the output column. It's actually a bit redundant to output the
# score and the window_filter, and since the score is more informative, I'm
# leaving the logical to be optional

# cat(sprintf("Window size is: %s\nWindow filtering threshold is: %s\n",
#             window_size, window_filtering_threshold))

stopifnot(is.integer(window_size) && window_size > 2)
if (!is.null(window_filtering_threshold)) {
    stopifnot(is.numeric(window_filtering_threshold) &&
              window_filtering_threshold > 0 &&
              window_filtering_threshold < 1)
}

if (is.null(input_file)) {
    stop("No input file provided")
}

required_files <- c(input_file)

# check whether the output file is in a subdir, and that the directory exists
if (!(dir.exists(dirname(output_file)))) {
    stop("Output file's directory does not exist. You have to create it first, or give a valid path.")
}

for (i in required_files) {
    if(!file.exists(i)) {
    stop(sprintf("File %s does not exist; check the path of the input file(s).\n",
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

# If there are NA rows. I have to derive the contig and pos
if (any(is.na(HD$contig))) {

    # NOTE: !weak spot of the script!
    get_contig_pos <- function(string) {

        contig <- gsub("(.*)_([0-9]+)$", "\\1", string)
        pos    <- gsub("(.*)_([0-9]+)$", "\\2", string)
        return(c(contig, pos))

    }

    HD[is.na(HD$contig), c("contig", "pos")] <-
        do.call(rbind,
                lapply(HD[is.na(HD$contig), "SNP.missing"],
                       get_contig_pos))

}

HD$pos <- as.numeric(HD$pos)
stopifnot(!is.na(HD$pos))

HD <- HD[order(HD$contig, HD$pos), ]

# By default, everything is retained
HD$filter <- "retained"

# We tag "excluded" all those with bad properties based on H and D 
# TODO: check the effects of H and D parameters/add this as an option probably
HD$filter[HD$H > H_upper_limit |
          HD$D > D_upper_limit | 
          HD$D < D_lower_limit] <- "excluded"

# positions with NAs are "unknown"
HD$filter[is.na(HD$H)] <- "unknown"

if (unknown_as_good) {
    # In some circumnstances, we might want to treat unknowns as good positions,
    # the reasonning being that HD can't excluded them, so they are good. One key
    # aspect of this reasonning is that it seems that HD fails to compute the
    # scores when there is a very low level of He on a locus. Note that this has
    # not been thoroughly investigated, so I do not put this as the default option.

    HD$filter[is.na(HD$H)] <- "retained"
}

HD$filter <- as.factor(HD$filter)

total_number_contigs <- length(unique(HD$contig))

HD$contig <- as.factor(HD$contig)
list_subsets <- split(HD, HD$contig)

all_scores <-
    lapply(seq_along(list_subsets), function(x) {

               cat(sprintf("Proceeding with contig number: %s on %s\n",
                           x, total_number_contigs))

               # Computation is done here: return the proportion of GOOD SNPs
               # in a window centered on each pos considered
                eval_Windows(list_subsets[[x]],
                             w_size = window_size)
       })

all_scores_2 <- unlist(all_scores)

stopifnot(nrow(HD) == length(all_scores_2))

HD$window_score_proportion_good_snp <- all_scores_2

if (!is.null(window_filtering_threshold)) {

    # A position is retained if it was not "excluded": keeping "retained" and "unknown"
    HD$window_filter_TrueIsRetained <-
        all_scores_2 > window_filtering_threshold & (HD$filter != "excluded")

    # if a window filterering threshold is used, output also a summary of the
    # results obtained, makes it convenient to check the results of a run

    number_snp_retained_window   <- sum(HD$window_filter_TrueIsRetained)
    number_snp_retained_noWindow <- sum(HD$filter != "excluded")
    prop_snp_retained_window     <- (number_snp_retained_window / nrow(HD))
    prop_snp_retained_noWindow   <- (number_snp_retained_noWindow / nrow(HD))

    cat(paste(basename(input_file),
              window_size,
              window_filtering_threshold,
              number_snp_retained_window,
              number_snp_retained_noWindow,
              prop_snp_retained_window,
              prop_snp_retained_noWindow,
              sep = "\t"),
        "\n")

}

write.table(HD, file = output_file, sep = "\t", row.names = F, quote = F)
