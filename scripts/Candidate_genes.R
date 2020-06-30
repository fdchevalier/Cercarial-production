#!/usr/bin/env Rscript
# Title: 
# Version: 
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 
# Modified in: 



#==========#
# Comments #
#==========#

# Getting the 1.8 LOD coordinates of the chromosomes around each QTL peak



#==========#
# Versions #
#==========#

# v0.0 - 2016-02-18: creation

#===================#
# Packages required #
#===================#

library("qtl")      # For R/qtl commands
library("vcfR")
library("magrittr")


##-------------------
## Packages
##-------------------

#library(latticeExtra)

#===========#
# Functions #
#===========#

# Working directory
#setwd(file.path(getwd(), "scripts"))



#===========#
# Variables #
#===========#

# Folders
data_fd    <- "../data/"
graph_fd   <- "../graphs/"
result_fd1 <- "../results/1-QTL/"
result_fd2 <- "../results/2-Candidate genes/"

myvcf_f <- paste0(data_fd, "calling/cerc_prod_snpEff_1-3-5.vcf.gz")

myqtl.ls <- paste0(result_fd1, "myqtl.ls.RData")

myfmt <- matrix(c(
            "GT", FALSE,
            "GQ", TRUE,
            "DP", TRUE),
            ncol=2, byrow=TRUE)


# # Object contaning the input file name
# filename <- basename(myvcf_f) %>% strsplit(., ".vcf.gz")

# # PT data file
# myF2.ptf <- paste0(data_fd, "phenotyping/F2.csv")

# # Output name of the R/qtl GT table without any filename extension
# myF2.gt <- "F2_geno"

#-----------------#
# Cross variables #
#-----------------#

# F2A variables
mypA <- c("SmBRE4_m", "SmLE19_f")   # Be careful at the order of the parents (cf. Allele code below)
myF1A <- c("F1A")
myF2A <- c("F2A")

#mypB <- c("SmLE15_m", "SmBRE2_f")
mypB <- c("SmBRE2_f", "SmLE15_m")   # Be careful at the order of the parents (cf. Allele code below)
myF1B <- c("F1B")
myF2B <- c("F2B")

myF2.list <- c("F2A", "F2B")

# GQ and read depth (rd) and missing data thresholds
## Assign NULL to skip to the variable to skip the corresponding filtering step
gq.trsh <- 30
dp.trsh <- 6

na.string <- NA # Could be "./."

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")


#--------------#
# QTL analysis #
#--------------#


# Number of permutationq to get threshold
mylod.trsh <- 0.05



#=================#
# Data processing #
#=================#

#---------------#
# Sanity checks #
#---------------#

if (mymissing.data < 0 | mymissing.data >= 1) { stop("mymissing.data must be in [0;1[") }
if (! any(grepl(pheno.cln, myprefixes[,1]))) {stop("The pheno.cln value is unknown in myprefixes table.", call.=FALSE)}

if (! dir.exists(result_fd)) { dir.create(result_fd, recursive = TRUE)

#--------------#
# Data Loading #
#--------------#

# Load QTL data
load(myqtl.ls)

mypeaks    <- summary(myqtl.ls[[3]][["mr"]]$lod, perm=myqtl.ls[[3]][["mr"]]$perm, alpha=mylod.trsh)

myqtl.chr <- mypeaks[rownames(mypeaks[grep("_[0-9ZW]$", mypeaks[,1], perl=TRUE), ]), 1]

lodint.ls <- lapply(myqtl.chr, function(x) lodint(myqtl.ls[[3]][[2]]$lod, chr = x, drop = 1.8))


# Load VCF data
myvcf <- read.vcfR(myvcf_f)
myvcf.bk <- myvcf

# Retains only the cross parents
myvcf@gt <- myvcf@gt[, colnames(myvcf@gt) %in% c("FORMAT", mypA, mypB)]


# Select markers within QTLs
myvec <- sapply(lodint.ls, function(x) getFIX(myvcf)[, 1] == x[1,1] & (getFIX(myvcf)[,2] >= x[1,2] & getFIX(myvcf)[,2] <= x[3,2])) %>% rowSums(.) > 0
myvcf <- myvcf[myvec, ]

# Select site that are all NAs
myvec.na <- extract.gt(myvcf, "GQ") %>% apply(., 1, function(x) all(is.na(x)))

# Select sites that are all below DP threshold
myvec.dp <- (extract.gt(myvcf, "DP") < dp.trsh) %>% apply(., 1, function(x) all(x, na.rm=TRUE))

# Select sites that are all below GQ threshold
myvec.gq <- (extract.gt(myvcf, "GQ") < gq.trsh) %>% apply(., 1, function(x) all(x, na.rm=TRUE))

# Select uninformative markers (same as reference or between all parents)
myvec.mrkr <- extract.gt(a) %>% gsub("|", "/", ., fixed = TRUE) %>% apply(., 1, function(x) all(x %in% "1/1", na.rm=TRUE) | all(x %in% "0/0", na.rm=TRUE))

# Remove sites that are below thresholds
myvec <- rowSums(cbind(myvec.na, myvec.dp, myvec.gq, myvec.mrkr)) > 0
myvcf <- myvcf[! myvec, ]




#---------------------
# Loading data table
#---------------------

mydata <- read.table("tables/combination.em.tsv", dec=".", sep="\t", header=T)

#------------------------------
# subseting the initial table
#------------------------------

# Order LOD table
mydata.ordered <- mydata[order(mydata[,3], decreasing=T),]

# QTL on x fist chr
x <- 4
mychroms <- unique(mydata.ordered[,1])[1:x]

for (c in mychroms) {
	# New data tables corresponding to each chromosome of interest
	mychrom <- mydata[ mydata[,1] == c, ]

	# Finding the maximum LOD value for each chromosome of interest
	mychrom.max <- max(mychrom$lod)

	# Calculate the max lod value - 1.8 lod
	mychrom.1.8 <- mychrom.max - 1.8
#	mychrom.1.8 <- mychrom.max

	# Searching all the row of the chromosome table that have a lod value >= to the max lod -1.8 lod
	mychrom.1.8 <- mychrom[ mychrom[,3] >= mychrom.1.8 , ]
	mychrom.1.8 <- mychrom.1.8[c(1,nrow(mychrom.1.8)),]
	print(mychrom.1.8)
	
	# Subtract position to get segment length
	myseg <- diff(as.numeric(as.vector(mychrom.1.8[,2])))
	
	# Average number of gene by bp
	mygenes <- 10000 / 364500000
	
	# Number genes expected in segment
	mygenes.nb <- myseg * mygenes
	
	# Print final message
	print(paste0("Number of expected genes under QTL of ", c, ": ", mygenes.nb))

}
