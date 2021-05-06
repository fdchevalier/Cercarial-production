#!/usr/bin/env Rscript
# Title: Candidate_genes.R
# Version: 0.1
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2021-01-12
# Modified in: 2021-05-06



#==========#
# Comments #
#==========#

# Generate list of candidate gene under the 1.8 LOD confidence interval of each QTL peak.



#==========#
# Versions #
#==========#

# v0.1 - 2021-05-06: reshape code / improve chromosome selection
# v0.0 - 2021-01-12: creation



#===================#
# Packages required #
#===================#

cat("Loading packages...\n")

suppressMessages({
    library("qtl")          # For R/qtl commands
    library("vcfR")
    library("magrittr")

    library("parallel")
    library("doParallel")

    library("rtracklayer")
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

source("functions/Candidate_genes_func.R")


#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd    <- "../data/"
graph_fd   <- "../graphs/"
result_fd1 <- "../results/2-QTL/"
result_fd2 <- "../results/3-Candidate genes/"

myvcf_f  <- paste0(data_fd, "calling/cerc_prod_snpEff_1-2-3-4-5.vcf.gz")
mygff_f  <- paste0(data_fd, "genome/schistosoma_mansoni.PRJEA36577.WBPS14.annotations.gff3")
myann_f  <- paste0(data_fd, "genome/Sm_transcript_table_gff-hhpred_2020-09-21.tsv")

# Expression data
myexpr_f <- paste0(data_fd, "genome/TPM_isoforms_Sm_2020-08-07.tsv")
myexpr.cln <- 6:7    
myexpr.nm  <- paste0("TPM ", c("sp shedder", "cercariae"))

mypheno.qtl.ls <- paste0(result_fd1, "mypheno.qtl.ls.RData")
myqtl.i <- paste0(result_fd1, "scantwo_qtl.i.RData")
st.perm <- paste0(result_fd1, "scantwo_perm.RData")

myfmt <- matrix(c(
            "GT", FALSE,
            "GQ", TRUE,
            "DP", TRUE),
            ncol=2, byrow=TRUE)


##-----------------#
## Cross variables #
##-----------------#

# F2A variables
mypA <- c("SmBRE4_m", "SmLE19_f")   # Be careful at the order of the parents (cf. Allele code below)
mypB <- c("SmBRE2_f", "SmLE15_m")   # Be careful at the order of the parents (cf. Allele code below)

# GQ and read depth (dp) thresholds
## Assign NULL to skip to the variable to skip the corresponding filtering step
gq.trsh <- 30
dp.trsh <- 6


#--------------#
# QTL analysis #
#--------------#

# Number of permutationq to get threshold
mylod.trsh <- 0.05

# Drop
mydrop <- 1.8



#=================#
# Data processing #
#=================#

#---------------#
# Sanity checks #
#---------------#

if (! dir.exists(result_fd2)) { dir.create(result_fd2, recursive = TRUE) }


#--------------#
# Data loading #
#--------------#

cat("Loading data...\n")

mygff <- readGFF(mygff_f)
myann <- read.csv(myann_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Expression data
myexpr <- read.csv(myexpr_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)

# Load QTL data
load(mypheno.qtl.ls)
myqtl.ls <- mypheno.qtl.ls[["average"]]

mylod    <- myqtl.ls[["combination"]][["em"]]$lod

load(myqtl.i)
load(st.perm)


#---------------#
# LOD intervals #
#---------------#

cat("Loading data...\n")

# Get the thresholds as a vector
thres <- summary(st.perm) %>% as.data.frame()
colnames(thres) <- (summary(st.perm) %>% attributes())$"names"
thres_vec <- thres[1,]

# Chromosomes of interest
myqtl.chr <- summary(myqtl.i, thresholds = as.numeric(thres_vec[1:5])) %>% as.data.frame() %>% .[,1:2] %>% unlist() %>% unique() %>% grep("_[1-7]$|_ZW$", ., value=TRUE) %>% sort()

# LOD intervals
lodint.ls <- lapply(myqtl.chr, function(x) lodint(myqtl.ls[[3]][[2]]$lod, chr = x, drop = mydrop))


#----------------#
# VCF processing #
#----------------#

cat("Processing VCF file...\n")

# Load VCF data
myvcf <- read.vcfR(myvcf_f)

# Retain only the cross parents
myvcf@gt <- myvcf@gt[, colnames(myvcf@gt) %in% c("FORMAT", mypA, mypB)]

# GT, GQ and DP data file
mydata <- vector("list", length(myfmt[,1])+1)
names(mydata) <- c("FIX", myfmt[,1])

mydata[["FIX"]] <- getFIX(myvcf, getINFO = TRUE)

for (i in myfmt[,1]) {
    mynum <- myfmt[ myfmt[,1] %in% i, 2]
    mydata[[i]] <- extract.gt(myvcf, i, as.numeric = mynum)
}

# Filtering on GQ
mat.gq <- mydata$GQ < gq.trsh

# Filtering on DP
mat.dp <- mydata$DP < dp.trsh

mat.sel <- apply(mat.gq * mat.dp, 2, as.logical)

mydata$GQ[mat.sel] <- NA
mydata$GT[mat.sel] <- NA
mydata$DP[mat.sel] <- NA


#------------------#
# Preparing tables #
#------------------#

cat("\nPreparing input tables...\n")

# Select markers within QTLs
myvec.mrkr <- sapply(lodint.ls, function(x) getFIX(myvcf)[, 1] == x[1,1] & as.numeric(getFIX(myvcf)[,2]) >= as.numeric(x[1,2]) & as.numeric(getFIX(myvcf)[,2]) <= as.numeric(x[3,2])) %>% rowSums(.) > 0

# Select sites that are all NAs
myvec.na <- mydata$GT %>% apply(., 1, function(x) all(is.na(x)))

# Remove sites that are below thresholds
myvec <- rowSums(cbind(! myvec.mrkr, myvec.na)) > 0

# Clean the lists
mydata <- lapply(mydata, function(x) x[! myvec,])

# Select uninformative markers (same as reference or between all parents)
myvec.uninfo <- mydata$GT %>% gsub("|", "/", ., fixed = TRUE) %>% apply(., 1, function(x) strsplit(x, "[/]") %>% unlist) %>% lapply(., function(x) unique(x) %>% length(.) == 1) %>% unlist()

# Remove sites that are below thresholds
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

mydata.tb <- data.frame(mydata$FIX, mydata$GT[, c(mypA, mypB)], mysc.cln, lod=mylod.cln)


#-----------------------#
# Candidate gene tables #
#-----------------------#

cat("Generating candidate gene tables...\n")

tb.ls <- vector("list", length(myqtl.chr))
names(tb.ls) <- myqtl.chr

for (i in myqtl.chr) {
    tb.ls[[i]]$qtl.tb <- qtl.tb(mydata.tb, chr = i, myann, myexpr, myexpr.cln, expr.nm = myexpr.nm, cln.print = 9:15)
    tb.ls[[i]]$qtl.tb.sum <- qtl.tb.sum(tb.ls[[i]]$qtl.tb, myann, myexpr, myexpr.cln, myexpr.nm, mygff, baits.bed=NULL)

    write.table(tb.ls[[i]]$qtl.tb, paste0(result_fd2, "/", i, "_QTL_table_details.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
    write.table(tb.ls[[i]]$qtl.tb.sum, paste0(result_fd2, "/", i, "_QTL_table_summary.tsv"), row.names = FALSE, quote = FALSE, sep = "\t")
}

