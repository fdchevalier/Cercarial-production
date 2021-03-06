#!/usr/bin/env Rscript
# Title: gt2rqtl.R
# Version: 0.7
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2016-02-18
# Modified in: 2020-08-02



#==========#
# Comments #
#==========#

# Transform GT table from vcftools to R/qtl csvsr table



#==========#
# Versions #
#==========#

# v0.7 - 2020-08-02: physical distance info added
# v0.6 - 2020-06-24: na.string option added / simplified unassembled scaffold option removed / allele reordering improved
# v0.5 - 2016-07-16: adapt script to the Sm genome v6
# v0.4 - 2016-06-21: simplified unassembled scaffold ordered to have SC at the end
# v0.3 - 2016-03-13: simplified unassembled scaffold name improved
# v0.2 - 2016-03-07: simplified unassembled scaffold name option added
# v0.1 - 2016-03-01: use of list for looking a different crosses at once added / new method for counting variants that passed tests
# v0.0 - 2016-02-18: creation



#===================#
# Packages required #
#===================#

# library()



#===========#
# Functions #
#===========#

gt2rqtl <- function(x, parents.cln, F1.cln, F2.cln, alleles=c("A","B"), out.fmt="csvs", out.name, na.string=NA, add.distance=FALSE) {

    #-------#
    # Usage #
    #-------#

    # x                     dataframe of the genotype. This should typically be the GT table generated from vcftools.
    # parents.cln           vector or list of vectors containing the column numbers or names of the parents.
    # F1.cln                vector or list of vectors containing the column numbers or names of the F1 individuals. (optional)
    # F2.cln                vector or list of vectors containing the column numbers or names of the F2 individuals.
    # alleles               alleles used to code genotypes. They are the same as the default of the qtl package.
    # out.fmt               file format used to export data. Might be csvs or csvsr. For more details, see R/qtl manuel and book.
    # out.name              base of the filename for exporting data. (optional)
    # simplify.unass.name   if TRUE, simplify the name of unassembled scaffolds to get only one unassembled scaffold. This is specific of S. manosni genome.
    # na.string             what NA value should be (character or NA).
    # add.distance          add physical distance (in Mb) to the table.


    #---------------------#
    # Checking parameters #
    #---------------------#

    # Check if mendatory arguments present
    if (missing("x")) {stop("x is missing.")}
    if (missing("parents.cln")) {stop("parents.cln is missing.")}
    if (missing("F2.cln")) {stop("F2.cln is missing.")}
    if (missing("out.name")) {stop("out.name is missing.")}

    # Check nature of object
    if (! is.data.frame(x)) {stop("x must be a dataframe.")}
    if (! is.vector(parents.cln) | length(parents.cln) != 2) {stop("parents.cln must be a vector of two elements.")}
    if (exists("F1.cln") & ! is.vector(F1.cln)) {stop("F1.cln must be a vector.")}
    if (! is.vector(F2.cln)) {stop("F2.cln must be a vector.")}
    if (! is.vector(out.fmt) | length(out.fmt) != 1) {stop("out.fmt must be a vector of one element.")}
    if (! is.vector(out.name) | length(out.name) != 1) {stop("out.name must be a vector of one element.")}
    if (! (is.character(na.string) | is.numeric(na.string) | is.na(na.string))) {stop("out.name must be a vector of one element.")}
    if (! is.logical(add.distance)) {stop("add.dsitance must be logical.")}

    # Check if out.fmt is OK
    if (length(grep("^csvs$|^csvsr$", out.fmt)) != 1)  {stop("out.fmt must have the value \"csvs\" or \"csvsr\".")}

    # Check list lengths if list used
    if (any(is.list(parents.cln) | is.list(F1.cln) | is.list(F2.cln))) {
        if (exists("F1.cln") & any(! is.list(parents.cln) | ! is.list(F1.cln) | ! is.list(F2.cln)) & any(is.vector(parents.cln) | is.vector(F1.cln) | is.vector(F2.cln))) {
            stop("parents.cln, F1.cln and F2.cln must be either vectors or list of vectors but not a mix.")
        } else if (any(! is.list(parents.cln) | ! is.list(F2.cln)) & any(is.vector(parents.cln) | is.vector(F2.cln))) {
            stop("parents.cln and F2.cln must be either vectors or list of vectors but not a mix.")
        } else if (length(parents.cln) != length(F2.cln)) {
            stop("parents.cln and F2.cln must have the same length.")
        }
        
        # Check for overlaps
        if (length(unlist(parents.cln)) != length(unique(unlist(parents.cln)))) {
            warning("The different vectors in the parents.cln have common columns.", immediate.=TRUE)
        }
        if (exists("F1.cln") & length(unlist(F1.cln)) != length(unique(unlist(F1.cln)))) {
            warning("The different vectors in the F1.cln have common columns.", immediate.=TRUE)
        }
        if (length(unlist(F2.cln)) != length(unique(unlist(F2.cln)))) {
            warning("The different vectors in the F2.cln have common columns.", immediate.=TRUE)
        }
    }

    if (!is.list(parents.cln)) {parents.cln <- list(parents.cln)}
    if (exists("F1.cln") & ! is.list(F1.cln)) {F1.cln <- list(F1.cln)}
    if (! is.list(F2.cln)) {F2.cln <- list(F2.cln)}


    #---------------#
    # Genotype data #
    #---------------#
    
    # Convert cell content in character string. This step is crucial for the replacement of the genotype.
    x <- apply(x, 2, as.character)
    
    # Replace real NA by string
    if (is.na(na.string)) {
        na.string <- "./."
        x[is.na(x)] <- na.string
    }


    for (i in 1:length(parents.cln)) {
        # Store column number of the given cross
        parents.cln.tmp <- parents.cln[[i]]
        if (exists("F1.cln")) {F1.cln.tmp <- F1.cln[[i]]}
        F2.cln.tmp <- F2.cln[[i]]
        
        cat(paste("Allele report for cross ",colnames(x)[parents.cln.tmp[1]]," x ",colnames(x)[parents.cln.tmp[2]],":\n", sep=""))
        nb.alleles <- sum(rowSums(x[,F2.cln.tmp] == na.string) != length(F2.cln.tmp))
        cat(paste("\t- Initial nunber present in the GT dataframe:", nb.alleles, "\n"))
        
        # Transform non alternative fixed alles from the parents in NA values
        x[ ! ((x[,parents.cln.tmp[1]] == "0/0" & x[,parents.cln.tmp[2]] == "1/1") | (x[,parents.cln.tmp[1]] == "1/1" & x[,parents.cln.tmp[2]] == "0/0")) , c(parents.cln.tmp,F2.cln.tmp) ] <- na.string
        nb.alleles <- sum(rowSums(x[,F2.cln.tmp] == na.string) != length(F2.cln.tmp))  
        cat(paste("\t- Alternative fixed alleles:", nb.alleles, "\n"))
        
        # Check for Mendelian segregation (ie, heterozygous state in F1)
        if (exists("F1.cln.tmp") & ! (any(is.na(F1.cln.tmp)) | any(is.null(F1.cln.tmp)))) {
            x[ ! (x[,F1.cln.tmp[1]] == "0/1" & x[,F1.cln.tmp[2]] == "0/1") , c(parents.cln.tmp,F2.cln.tmp) ] <- na.string
            nb.alleles <- sum(rowSums(x[,F2.cln.tmp] == na.string) != length(F2.cln.tmp))  
            cat(paste("\t- Alleles following Mendelian segregation:", nb.alleles, "\n"))
        }

        # Search and replace GT with corresponding alleles 
        x[x[,parents.cln.tmp[1]] == "0/0", F2.cln.tmp] <- apply(x[x[,parents.cln.tmp[1]] == "0/0", F2.cln.tmp], 2, function(y) {gsub("0", alleles[1], y)})
        x[x[,parents.cln.tmp[1]] == "0/0", F2.cln.tmp] <- apply(x[x[,parents.cln.tmp[1]] == "0/0", F2.cln.tmp], 2, function(y) {gsub("1", alleles[2], y)})
        x[x[,parents.cln.tmp[1]] == "1/1", F2.cln.tmp] <- apply(x[x[,parents.cln.tmp[1]] == "1/1", F2.cln.tmp], 2, function(y) {gsub("1", alleles[1], y)})
        x[x[,parents.cln.tmp[1]] == "1/1", F2.cln.tmp] <- apply(x[x[,parents.cln.tmp[1]] == "1/1", F2.cln.tmp], 2, function(y) {gsub("0", alleles[2], y)})
    
        # Remove multi-allelic sites
        x.lines <- apply(x[, F2.cln.tmp], 1, function(y) {any(grepl("[[:digit:]]", y))})
        x[x.lines, c(parents.cln.tmp,F2.cln.tmp)] <- na.string
        nb.alleles <- sum(rowSums(x[,F2.cln.tmp] == na.string) != length(F2.cln.tmp))  
        cat(paste("\t- Bi-allelic sites:", nb.alleles, "\n"))
    }

    cat("Processing table...\n")

    # Compile all F2 column
    F2.cln <- unique(unlist(F2.cln))
    
    # Remove uninformative variants
    x <- x[ rowSums(x[,F2.cln] == na.string) != length(F2.cln) , ]
    cat(paste("\t- The final table will contain",nrow(x),"variants.\n"))

    # Replace na.string with real NA
    x[,F2.cln] <- apply(x[,F2.cln], 2, function(y) {gsub(na.string, NA, y, fixed=TRUE)})
    
    # Remove /
    x[,F2.cln] <- apply(x[,F2.cln], 2, function(y) {gsub("/", "", y, fixed=TRUE)})
    
    # Reorder genotype (BA -> AB)
    ## Create all genotypes from alleles
    genos <- expand.grid(rep(list(alleles), length(alleles)))
    ## Remove homozygous
    genos <- genos[apply(genos, 1, function(x) unique(x) %>% length()) > 1, ]
    ## Create string genotypes
    genos <- apply(genos, 1, function(x) c(paste0(x, collapse = ''), sort(x) %>% paste0(., collapse = ''))) %>% t()
    ## Keep only those that need to be reordered
    genos <- genos[apply(genos, 1, function(x) unique(x) %>% length()) > 1, , drop = FALSE]
    ## Replace genotypes that need to be reorder
    for (g in 1:nrow(genos)) { x[,F2.cln] <- apply(x[,F2.cln], 2, function(y) {gsub(genos[g,1], genos[g,2], y, fixed=TRUE)}) }
    
    # Replace missing data with NA
    x[!nzchar(x)] <- "NA"

    # Remove spaces 
    x[,c(1:2,F2.cln)] <- apply(x[, c(1:2,F2.cln)], 2, function(y) {gsub(" *","", y)})

    # Final table
    dist   <- NULL
    if (add.distance) {
        dist    <- as.numeric(x[,2]) / 1e6
    }

    x <- cbind(paste(1:nrow(x),x[,1], x[,2], sep="_"), x[,1], dist, x[,F2.cln])

    if (add.distance) {
        colnames(x)[1:3] <- c("id", "", "")
    } else {
        colnames(x)[1:2] <- c("id", "")
    }


    if (out.fmt == "csvsr") {
        myrow.nm <- FALSE
        mycln.nm <- TRUE
    } else {
        x <- t(x)
        myrow.nm <- TRUE
        mycln.nm <- FALSE
    }

    # Write file
    cat("Writting output table...\n")
    write.table(x, file=paste(out.name,out.fmt,sep="."), row.names=myrow.nm, col.names=mycln.nm, quote=FALSE, sep=",")
}
