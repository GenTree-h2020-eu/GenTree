# Copyright (C) 2021, 2022 Katharina B. Budde
# SPDX-License-Identifier: GPL-3.0-only
# The use and distribution terms for this software are covered by the
# GNU General Public License v3.0 only (https://www.gnu.org/licenses/gpl-3.0-standalone.html)
# which can be found in the file LICENSE at the root of this distribution.
# By using this software in any fashion, you are agreeing to be bound by
# the terms of this license.
# You must not remove this notice, or any other, from this software.
#
# Description: The purpose of this script is to estimate nucleotide diversity (pi) per population based on a set of folded site frequency spectra (fsfs) derived from resampling stored in tsv files. In our experiments 1000 such files were utilized, hence the name of the script. 
#
# This script expects your data organized as follows,
# species -|
#          |- populations  -|
#                           |- sitetypes -|
#                                         |- FSFS Files
# e.g., ES/ES_PP_01/4fold/file1.fsfs. Heads up: this script does not check if there are in fact fsfs files, so please keep the content of the directories pure. 
#
# To generate the documentation of this file just run `document("pi_1000resamples.R", check_package=FALSE)' in R, provided the package `document' was loaded. 

library(gmodels)
library(data.table)

#' Compute pi from single fsfs file
#' 
#' Given a \code{file} fsfs file, as well as the currently computed \code{sitetype} and \code{species}
#' 
#' @param file The fsfs file containing the observed folded site frequency spectrum. 
#' @param sitetype The sitetype this the data from \code{file} belongs to (used for referencing TODO@Katha Number)
#' @param species The species the data from \code{file} belongs to (used for referencing TODO@Katha Number)
#' @return The pi value computed
#' @export
#' @examples process_file(xyz.fsfs,"4fold","Ppinaster")
process_file <- function(file, sitetype,  species){
    # read data table from transparent unzipped folded site frequency spectrum file
    freq.data <- read.table(gzfile(file), skip = 3, header =TRUE, sep ='\t') # extracts data frame
    # extract sample size number from 2nd line in file
    sample.size <- read.table(gzfile(file), skip = 1, # extract sample size from file
                              nrows= 1,
                              header =FALSE, sep =':',comment.char = "")$V2
    # compute allele frequency by deviding column entries in bin by sample size
    freq.data$Freq1 <- freq.data$bin / sample.size
    # compute pi per site: 2 * p * (1 - p) * bessel correction * count-per-allele-freq
    freq.data$pi.site <- 2 * freq.data$Freq1 * (1 - freq.data$Freq1) *
        (sample.size) / (sample.size - 1) * freq.data$count
    #  estimate average pi
    pi <- sum(freq.data$pi.site) / (sum(freq.data$count) + monomorphic_sites[Species==species,get(sitetype)])
    return(pi)
}

#' Process whole sitetype for species
#' 
#' Given a \code{species}, a \code{population}, and a \code{sitetype}, this function calls \link{process_file} for any file found in the species/population/sitetype directory and stores the resulting pi value (per fsfs file) into a data.table data structure, which will be returned. 
#' 
#' @param sitetype The sitetype to be computed refers to the annotation of the site: all sites (all); intergenic, intron, and 4fold (putative neutral); 0fold (non-synonymous), 4fold (synonymous)
#' @param species The species the data from \code{file} belongs to (e.g. Ppinaster)
#' @param population The population the data from \code{file} belongs to (e.g. ES_PP_01)
#' @return \code{data.table} containing pi values for all fsfs files in species/population/sitetype directory
#' @export
#' @examples process_sitetype("4fold","ES_PP_01", "Ppinaster")
process_sitetype <- function(sitetype, population, species){
    # Find fsfs files in path species/population/sitetype
    thefiles <- list.files(file.path(species, population, sitetype),full.names=TRUE)
    colnames <- c("file","pi")
    # Allocate data table
    data <- data.table(replication=1:length(thefiles), pi=-1)
    # Iterate through all (fsfs) files (stored in thefiles)
    for (f in thefiles){
        # Extract file number by 1st removing start characters and 2nd removing end chars
        filenumber <- as.numeric(sub(".fsfs.gz","", sub(".*replicate_","", f)))
        data[filenumber,2] <- process_file(f, sitetype, species)
    }
    return(data)
}

#' Estimate average pi values per population per site type
#' 
#' Given a \code{population}, a \code{species}, and \code{ourcolumns}, this function calls \link{process_sitetype} and averages pi over all replicate pi estimates and calculates the ratio of pi0fold/pi4fold for the population being processed. It also counts the number of pi estimates the average is based on.
#' 
#' @param population The population the data from \code{file} belongs to (e.g. ES_PP_01)
#' @param species The species the data from \code{file} belongs to (e.g. Ppinaster)
#' @param ourcolumns Name of the columns indicating the type of estimate.
#' @return \code{data.table} with average pi values, the ratio pi0/pi4 and the number of replicate pi estimates the average per population is based on.
#' @export
#' @examples process_population("ES_PP_01", "Ppinaster", "0fold, ..., pi0div4")
process_population <- function(population, species, ourcolumns ){
    # Allocate data table for result 
    statdata <-data.table(thesitetype = ourcolumns, vals=-1)
    savepivectors <- list() # create empty list
    for (sitetype in sitetypelist){ # iterate through sitetypes, i.e., 4fold, 0fold, etc.
        sitetypetable <- process_sitetype(sitetype, population, species) # table: replicate(file);pivalue
        # Compute simple statistics; 
        stats <- ci(sitetypetable$pi, confidence=0.95, na.rm=TRUE) 
        statdata[thesitetype==sitetype, "vals"] <- as.numeric(stats[1])
        savepivectors[[ sitetype ]] <- sitetypetable$pi
        number<-sum(!is.na(sitetypetable$pi))
        statdata[thesitetype==paste(sitetype,"-count",sep=""), "vals"] <-number
    }
    sitetypesdataframe <- data.frame(sapply(savepivectors,c))
    # Add pi0fold/pi4fold
    sitetypesdataframe$pi0div4 <- sitetypesdataframe$X0fold / sitetypesdataframe$X4fold
    # Write out intermediate result for documentation. 
    write.table(sitetypesdataframe,
                file=file.path(outputdir, paste(species,"-",population,sep="")),
                row.names = T)
    stats <- ci(sitetypesdataframe$pi0div4, confidence=0.95,alpha=1-0.95, na.rm=TRUE) 
    statdata[thesitetype=="pi0div4", "vals"] <- as.numeric(stats[1])
    # Return table with pi result values. 
    return(statdata)
}


#' Obtain population names
#' 
#' Helper function that, given a \code{species}, returns the list of directory names (population names).
#' 
#' @param species The species currently being processed (e.g. Ppinaster)
#' @return \code{list} containing a list of file (population) names
#' @export
#' @examples get_pop_names("Ppinaster")
get_pop_names <- function(species){
    #'Get the population(folder)names inside species(dir)
    pop_names <- list.files(species)
    return(pop_names)
}

#' Compute pi values for species
#'
#' Given a \code{species} and output column names (\code{ourcolums}), this function infers the names of the respective populations by calling \link{get_pop_names} and processes every population using \link{process_population}. 
#' 
#' @param species The species currently being processed (e.g. Ppinaster)
#' @param ourcolumns Name of the columns indicating the type of estimate.
#' @return boolean, if successful \code{TRUE} 
#' @export
#' @examples process_species("Ppinaster","0fold, ..., pi0div4")
process_species <- function(species, ourcolumns){
    speciesdata <-data.table(thesitetype = ourcolumns)          
    population_names <- get_pop_names(species) 
    for (population in population_names){
        result <- process_population(population, species, ourcolumns)
        setnames(result, "vals", population)
        speciesdata <- cbind(speciesdata,result[,2])
    }
    write.table(t(speciesdata), paste(species,"_pi_results",sep=""), sep=" ", col.names=TRUE, quote=FALSE)
    return(TRUE)
}

#' Process all species directories and print results
#'
#'
#' This function loops over the species names, defines the column names of the output table ("ourcolumns"), then it calls \link{process_species} and prints the results table. 
#' @export
process_all_species <- function(){
    for (species in specieslist){
        print(paste("Processing species:",species))
		ourcolumns <- c("0fold","0fold-count",
                        "4fold","4fold-count",
                        "all", "all-count",
                        "putative_neutral","putative_neutral-count",
                        "pi0div4")
        result <- process_species(species, ourcolumns)
        print(paste("Results:",result))
    }
}


# Define runtime parameters by provding species and site type names, as well as output directory. Additionally the number of monormophic sites per population and site type is defined.
specieslist <- c('Bpendula', 'Fsylvatica', 'Pabies', 'Pnigra', 'Ppinaster', 'Psylvestris', 'Qpetraea')
sitetypelist <- c("0fold","4fold","all","putative_neutral")
outputdir <- "replicate-pi-values"
monomorphic_sites <- data.table(
    "Species" = specieslist,
    "0fold" = c(1585956,1553756,1101532,1568165,623085,63214,1138138),
    "4fold" = c(353155,347814,243641,353048,141005,14316,240324),
    "all"   = c(4140170,4911524,3003470,4307690,1787835,282605,4705012),
    "putative_neutral" = c(2049096,2870440,1550060,1883003,968528,200383,2982473))
# access example: monomorphic_sites[Species=="Ppinaster","4fold"]

process_all_species()
