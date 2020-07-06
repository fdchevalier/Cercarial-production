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

library("parallel")
library("doParallel")

library("rtracklayer")

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

myvcf_f  <- paste0(data_fd, "calling/cerc_prod_snpEff_1-3-5.vcf.gz")
mygff_f  <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
myann_f  <- "~/data/sm_Gene_table/Sm_v7.1_ann/Sm_v7.1_transcript_table_gff-hhpred.tsv" 
# Expression data
myexpr_f <- "~/data/sm_Gene_table/TPM_isoforms_Sm_v7.1.tsv"
myexpr.cln <- 3:4    # Column containing juvenil and adult expression data
myexpr.nm  <- paste0("TPM ", c("sp 48h", "cercariae"))


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

if (! dir.exists(result_fd2)) { dir.create(result_fd2, recursive = TRUE) }

#--------------#
# Data Loading #
#--------------#

mygff <- readGFF(mygff_f)
myann <- read.csv(myann_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Expression data
myexpr <- read.csv(myexpr_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Load QTL data
load(myqtl.ls)

mylod     <- myqtl.ls[["combination"]][["em"]]$lod
myperm    <- myqtl.ls[["combination"]][["em"]]$perm

mypeaks   <- summary(mylod, perm=myperm, alpha=mylod.trsh)

myqtl.chr <- mypeaks[rownames(mypeaks[grep("_[0-9ZW]$", mypeaks[,1], perl=TRUE), ]), 1]

lodint.ls <- lapply(myqtl.chr, function(x) lodint(myqtl.ls[[3]][[2]]$lod, chr = x, drop = 1.8))


# Load VCF data
myvcf <- read.vcfR(myvcf_f)
myvcf.bk <- myvcf

# Retains only the cross parents
myvcf@gt <- myvcf@gt[, colnames(myvcf@gt) %in% c("FORMAT", mypA, mypB)]


# GT, GQ and DP data file
mydata <- vector("list", length(myfmt[,1])+1)
names(mydata) <- c("FIX", myfmt[,1])

mydata[["FIX"]] <- getFIX(myvcf, getINFO = TRUE)

for (i in myfmt[,1]) {
    mynum <- myfmt[ myfmt[,1] %in% i, 2]
    mydata[[i]] <- extract.gt(myvcf, i, as.numeric = mynum)
}

# GQ
mat.gq <- mydata$GQ < gq.trsh

# DP
mat.dp <- mydata$DP < dp.trsh

mat.sel <- apply(mat.gq * mat.dp, 2, as.logical)

mydata$GQ[mat.sel] <- NA
mydata$GT[mat.sel] <- NA
mydata$DP[mat.sel] <- NA


# Select markers within QTLs
myvec.mrkr <- sapply(lodint.ls, function(x) getFIX(myvcf)[, 1] == x[1,1] & as.numeric(getFIX(myvcf)[,2]) >= as.numeric(x[1,2]) & as.numeric(getFIX(myvcf)[,2]) <= as.numeric(x[3,2])) %>% rowSums(.) > 0

# Select site that are all NAs
myvec.na <- mydata$GT %>% apply(., 1, function(x) all(is.na(x)))

# Remove sites that are below thresholds
# myvec <- rowSums(cbind(! myvec.mrkr, myvec.na, myvec.dp, myvec.gq, myvec.uninfo)) > 0
myvec <- rowSums(cbind(! myvec.mrkr, myvec.na)) > 0

# Clean the lists
mydata <- lapply(mydata, function(x) x[! myvec,])

# Select uninformative markers (same as reference or between all parents)
myvec.uninfo <- mydata$GT %>% gsub("|", "/", ., fixed = TRUE) %>% apply(., 1, function(x) strsplit(x, "[/]") %>% unlist) %>% lapply(., function(x) unique(x) %>% length(.) == 1) %>% unlist()

# Remove sites that are below thresholds
# myvcf <- myvcf[! myvec.uninfo, ]
mydata <- lapply(mydata, function(x) x[! myvec.uninfo,])

# Find LOD score when available
mylod.cln <- apply(mydata$FIX[,1:2], 1, function(x) mylod[mylod[,1] == x[1] & mylod[,2] == x[2], 3]) %>% lapply(., function(x) if (length(x) == 0 ) { NA } else { x }) %>% unlist()

# Genotype score
mysc.cln <- matrix(NA, ncol=2, nrow=nrow(mydata$GT))
colnames(mysc.cln) <- c("A", "B")

myp <- list(A=mypA, B=mypB)

for(i in 1:length(myp)) {
    mydata.tmp <- mydata$GT[, myp[[i]]]
    mysc.cln[, i] <- apply(mydata.tmp, 1, function(x) {
                   
                                       myna  <- is.na(x)
                                       myhet <- is_het(as.matrix(x), na_is_false = FALSE)

                                       if (all(myna)) { # If all is NA
                                           y <- 0
                                       } else if (any(myna)) { # If one NA
                                           if (all(na.omit(myhet))) { # If remaining genotype heterozygous
                                               y <- 0
                                           } else if (any(grepl("0[/|]0", x))) { # If homozygous by deduction is all homozygous for the reference
                                               y <- 0
                                           } else { # If homozygous for the alternative
                                               y <- 5
                                           }
                                       } else if (all(! myna)) { # If no NA
                                           if (all(myhet)) { # If all heterozygous
                                               y <- 0
                                           } else if (any(myhet)) { # If one homozygous
                                               if (any(grepl("0[/|]0", x))) { # If homozygous site is the same as reference
                                                   y <- 0
                                               } else {
                                                   y <- 0 # 10
                                               }
                                           } else if (all(! myhet)) { # If no heterozygous
                                               if (strsplit(x, "[/|]") %>% unlist %>% unique %>% length() == 1) { # Same allele everywhere
                                                   y <- 0
                                               } else {
                                                   y <- 20
                                               }
                                           }
                                       }

                                       return(y)
                                })

}



# Build table
mydata.tb <- data.frame(mydata$FIX, mydata$GT[, c(mypA, mypB)], mysc.cln, lod=mylod.cln)

mygff.genes <- mygff[ mygff[,3] == "gene", ]

tb.ls <- vector("list", length(myqtl.chr))
names(tb.ls) <- myqtl.chr

for (i in myqtl.chr) {
    tb.ls[[i]]$qtl.tb <- qtl.tb(mydata.tb, chr = i, myann, myexpr, myexpr.cln, expr.nm=myexpr.nm, cln.print=9:15)
    tb.ls[[i]]$qtl.tb.sum <- qtl.tb.sum(tb.ls[[i]]$qtl.tb, myann, myexpr, myexpr.cln, myexpr.nm, mygff.genes,baits.bed=NULL)

    write.table(tb.ls[[i]]$qtl.tb, paste0(result_fd2, "/", i, "_QTL_table_details.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(tb.ls[[i]]$qtl.tb.sum, paste0(result_fd2, "/", i, "_QTL_table_summary.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}
