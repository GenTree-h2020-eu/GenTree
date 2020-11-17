#!/usr/bin/env Rscript

# vcf2sfs_resample.R

# Author Chedly Kastally <ckastall@gmail.com>
# Copyright (C) 2020 Chedly Kastally <ckastall@gmail.com>
# SPDX-License-Identifier: GPL-3.0-only
# Modified On 2020-05-15 12:32
# Created  2020-04-23 15:26

# Use AN and AC fields from a vcf file to compute the SFS by resampling n
# samples.

# NOTE:
#
# - This is meant to be used in a pipe with bcftools; and it's only with
# using bcftools that this script works atm
#
# eg:
#
# bcftools query -Hf "%CHROM\t%POS[\t%AN\t%AC]\n" input.vcf | \
#     ./vcf2sfs_resample.R -n <sample_size> > output_sfs.tsv

# TODO:
#
# - implement the option to read a vcf file? I don't think it's worth it; I
# would not use R

# library ----

library("dplyr")

# Optparse ----

posterior_n        <- 10
ascertainment_bias <- 0
folded_mode        <- FALSE

library("optparse")

option_list <-
    list(
         make_option(c("-i", "--input_file"),
                     type    = "character",
                     default = "stdin",
                     help    = "input_file [Default: %default]"),
         make_option(c("-a", "--ascertainment_bias"),
                     type    = "numeric",
                     default = ascertainment_bias,
                     help    = "Posterior sample size [Default: %default]"),
         make_option(c("-n", "--posterior_n"),
                     type    = "numeric",
                     default = posterior_n,
                     help    = "Posterior sample size [Default: %default]"),
         make_option(c("-f", "--folded_mode"),
                     action  = "store_true",
                     default = FALSE,
                     help    = "compute the folded SFS [Default: %default]"),
         make_option(c("-R", "--seed"),
                     type    = "numeric",
                     default = NULL,
                     help    = "set seed [Default: %default]")
    )

parser <-
    OptionParser(usage = "%prog -i input_file -o output_file",
                 option_list = option_list)

opt <- parse_args(parser)

input_file         <- opt$input_file
posterior_n        <- opt$posterior_n
ascertainment_bias <- opt$ascertainment_bias
seed               <- opt$seed
folded_mode        <- opt$folded_mode

# setup ----

if (!is.null(seed)) {
    set.seed(seed)
}

# main ----

con <- file(input_file, open = "r")

## header ----

header <- readLines(con, n = 1)

header <- unlist(strsplit(header, "\t"))
header <- gsub("#?\\s?\\[[0-9]+\\]", "", header)
header <- gsub(":.*", "", header)

i_an <- which(header == "AN")
i_ac <- which(header == "AC")

## Process the file line by line ----

res <- vector(mode = "list")
line_count <- 0
i <- 1

while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
        break
    }
    line_count <- line_count + 1
    vec <- unlist(strsplit(line, "\t"))[c(i_ac, i_an)]

    # cat(sprintf("current AC: %s, AN: %s\n", vec[1], vec[2]))
    if (!grepl(",", vec[1])) {

        ac <- as.numeric(vec[1])
        an <- as.numeric(vec[2])

        # Silently excludes sites for which posterior_n < an
        if (an >= posterior_n) {

            res[i] <- sum(sample(c(rep(0, an - ac), rep(1, ac)),
                              posterior_n,
                              replace = F))
            i <- i + 1
        } else {
            # Be vocal about those? At least for gentree it makes sense to be.
            cat(sprintf("# WARNING: dropping line %s with AN = %s\n",
                        line_count,
                        an))
        }

    }
}

close(con)

## output ----

sfs_table <- as.data.frame(table(as.numeric(res)),
                           stringAsFactors = F)

colnames(sfs_table) <- c("bin", "count")

sfs_table[sfs_table$bin == 0, "count"] <-
    sfs_table[sfs_table$bin == 0, "count"] + ascertainment_bias

cat(sprintf("# Parameters\n# Posterior sample size: %s\n", posterior_n))

if (folded_mode) {

    cat("# Folded mode\n")

    bins <- as.numeric(as.character(sfs_table$bin))

    sfs_table$folded_bins <-
        factor(ifelse(bins > posterior_n / 2,
                      (posterior_n) - bins,
                      bins))

    sfs_table <- sfs_table %>%
        dplyr::group_by(folded_bins) %>%
        dplyr::summarize(count = sum(count)) %>%
        as.data.frame(.)

    colnames(sfs_table) <- c("bin", "count")

}

write.table(sfs_table, sep = "\t", row.names = F, file = "", quote = F)

