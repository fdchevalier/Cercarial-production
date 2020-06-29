#!/usr/bin/env Rscript
# Title: Cercaire_intensity_qtl.R
# Version: 0.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2016-02-18
# Modified in: 2016-07-15



#==========#
# Comments #
#==========#

# Transform GT table from vcftools to R/qtl csvr table
# Derived from Chronobio_qtl.R



#==========#
# Versions #
#==========#

# v0.2 - 2016-07-15: missing data bug correction
# v0.1 - 2016-05-18: combined analysis included
# v0.0 - 2016-02-18: creation



#===================#
# Packages required #
#===================#

library("qtl")      # For R/qtl commands
library("vcfR")
library("magrittr")
library("snow")     # For parallel permutations



#===========#
# Functions #
#===========#

# Working directory
#setwd(file.path(getwd(), "scripts"))



#--------------------#
# External functions #
#--------------------#

source("Cercariae_intensity_qtl_func.R")

# Function specific to R/qtl
source("gt2rqtl.R")

# Function to plot data on exome
source("Sm.matplot.data.R")



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/1-QTL/"

myvcf_f <- paste0(data_fd, "calling/cerc_prod_snpEff_reduced.vcf.gz")

myfmt <- matrix(c(
            "GT", FALSE,
            "GQ", TRUE,
            "DP", TRUE),
            ncol=2, byrow=TRUE)


# Object contaning the input file name
filename <- basename(myvcf_f) %>% strsplit(., ".vcf.gz")

# PT data file
myF2.ptf <- paste0(data_fd, "phenotyping/F2.csv")

# Output name of the R/qtl GT table without any filename extension
myF2.gt <- "F2_geno"


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
rd.trsh <- 10
# mymis.data must be in [0;1[
mymissing.data <- 0.2

na.string <- NA # Could be "./."

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")


#--------------#
# QTL analysis #
#--------------#

n.cluster <- 3  # ATTENTION: increasing this number carefully (make sure you have sufficient RAM) otherwise server might be unusable

# Number of permutationq to get threshold
my.n.perm <- 1000   # ATTENTION: increasing number will increase computation time substantially
mylod.trsh <- 0.05

# Model for QTL scan
pheno.cln <- 8

# Prefix for graph
myprefixes <- matrix(c(
#    #pheno.cln, graph_prefix, scan_model
#        2,      "no_class",   "np",
#        3,      "2_classes",  "binary",
#        4,      "3_classes",  "np"), ncol=3, byrow=T)       # The pheno.cln cannot start at one because it is the ID column
    #pheno.cln, graph_prefix, scan_model
        7,      "sum",      "normal",
        8,      "average",  "normal"), ncol=3, byrow=T)       # The pheno.cln cannot start at one because it is the ID column
#        8,      "average",  "np"), ncol=3, byrow=T)       # The pheno.cln cannot start at one because it is the ID column

g.prefix <- myprefixes[ myprefixes[,1] == pheno.cln, 2 ]
mymodel <- myprefixes[ myprefixes[,1] == pheno.cln, 3]




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

myvcf <- read.vcfR(myvcf_f)

# Remove mitochondrial markers
myvcf <- myvcf[! grepl("mito", getFIX(myvcf)[,1], ignore.case = TRUE), ]

# Remove indels
myvcf <- extract.indels(myvcf)

# Keep only biallelic sites
myvcf <- myvcf[is.biallelic(myvcf), ]

# GT, GQ and DP data file
mydata <- vector("list", length(myfmt[,1]))
names(mydata) <- myfmt[,1]

for (i in myfmt[,1]) {
    mynum <- myfmt[ myfmt[,1] %in% i, 2]
    mydata[[i]] <- data.frame(getFIX(myvcf)[, 1:2], extract.gt(myvcf, i, as.numeric = mynum), stringsAsFactors = FALSE)
}

# Clean any phasing
mydata$GT[, -(1:2)] <- apply(mydata$GT[, -(1:2)], 2, function(x) gsub("|", "/", x, fixed = TRUE))

# Reverse phased data to default
if (any(grepl("1/0", as.matrix(mydata$GT[, -(1:2)]), fixed = TRUE))) {
    mydata$GT[, -(1:2)] <- apply(mydata$GT[, -(1:2)], 2, function(x) gsub("1/0", "0/1", x, fixed = TRUE))
}

# Remove uninformative markers (same as reference or between all parents)
myvec  <- apply(mydata$GT[, colnames(mydata$GT) %in% c(mypA, mypB)], 1, function(x) all(x == "0/0", na.rm=TRUE) | all(x == "1/1", na.rm=TRUE))
mydata <- lapply(mydata, function(x) x[! myvec, ])





#-----------------------------#
# Data filtering on DP and GQ #
#-----------------------------#
cat("\nFiltering data based on GQ and DP...\n")

# Converting data in numeric
cat("\t- Initial number of alleles before filtering in F2A and F2B:", nrow(mydata$GT), "\n")

for (i in myF2.list) {

    # Data filtering removing SNPs with low GQ
    if (! is.null(gq.trsh)) {        # gq.trsh is in the section variables
        gq.sign <- ">="
        mydata.gt.tmp.flt.gq <- mydata.fltr(mydata$GT, mydata$GQ, i, missing.data = mymissing.data, sign = gq.sign, trsh = gq.trsh)

        if (is.na(na.string)) {
            na.mt <- is.na(mydata.gt.tmp.flt.gq[,grep(i,colnames(mydata.gt.tmp.flt.gq))])
        } else {
            na.mt <- mydata.gt.tmp.flt.gq[,grep(i,colnames(mydata.gt.tmp.flt.gq))] == na.string
        }

        nb.alleles <- nrow(mydata.gt.tmp.flt.gq[(rowSums(na.mt)/length(grep(i,colnames(mydata.gt.tmp.flt.gq))) <= mymissing.data),])
        cat("\t- Remaining alleles after GQ filtering in ", i, ": ", nb.alleles, "\n", sep = "")
    } else {
        mydata.gt.tmp.flt.gq <- mydata$GT
    }

    # Data filtering removing SNPs with low read.depth
    if (! is.null(rd.trsh)) {
        rd.sign <- ">="
        mydata.gt.tmp.flt.gq.dp <- mydata.fltr(mydata.gt.tmp.flt.gq, mydata$DP, i, missing.data = mymissing.data, sign = rd.sign, trsh = rd.trsh, na.string = na.string)

        if (is.na(na.string)) {
            na.mt <- is.na(mydata.gt.tmp.flt.gq.dp[,grep(i,colnames(mydata.gt.tmp.flt.gq.dp))])
        } else {
            na.mt <- mydata.gt.tmp.flt.gq.dp[,grep(i,colnames(mydata.gt.tmp.flt.gq.dp))] == na.string
        }

        nb.alleles <- nrow(mydata.gt.tmp.flt.gq.dp[(rowSums(na.mt)/length(grep(i,colnames(mydata.gt.tmp.flt.gq.dp))) <= mymissing.data),])
        cat("\t- Remaining alleles after DP filtering in ", i,": ", nb.alleles, "\n", sep = "")
    } else {
        mydata.gt.tmp.flt.gq.dp <- mydata$GT
    }

#    # Rename data
#    assign(paste("mydata.gt.F2",sep=""), mydata.gt.tmp.flt.gq.dp)
#    assign(paste("my",i,".gt",sep=""), name.gt)

    mydata$GT <- mydata.gt.tmp.flt.gq.dp
}

# Suffix for output genotype file generated during next step
myflt.suffix <- ""
if (! is.null(gq.trsh) | ! is.null(rd.trsh)) {
    myflt.suffix <- ".flt"
    if (! is.null(gq.trsh)) {myflt.suffix <- paste(myflt.suffix,"-gt", gq.trsh, sep="")}
    if (! is.null(rd.trsh)) {myflt.suffix <- paste(myflt.suffix,"-dp", rd.trsh, sep="")}
}


#----------------------------------#
# Conversion of GQ in R/qtl format #
#----------------------------------#
cat("\nConverting GQ table in R/qtl format...\n")

# File name
myF2.gt  <- paste0(result_fd, myF2.gt, myflt.suffix)
myF2.gtf <- paste0(myF2.gt, ".csvs")


# Column associated with name
mypA <- sapply(mypA, function(x) grep(x, colnames(mydata$GT)))
myF1A <- grep(myF1A, colnames(mydata$GT))
myF2A <- grep(myF2A, colnames(mydata$GT))

mypB <- sapply(mypB, function(x) grep(x, colnames(mydata$GT)))
myF1B <- grep(myF1B, colnames(mydata$GT))
myF2B <- grep(myF2B, colnames(mydata$GT))

myp <- list(mypA, mypB)
myF1 <- list(myF1A, myF1B)
myF2 <- list(myF2A, myF2B)


# Table conversion
# if (file.exists("tables/") == FALSE) {dir.create("tables")}
gt2rqtl(mydata$GT, parents.cln = myp, F1.cln = myF1, F2.cln = myF2, out.fmt = "csvs", out.name = myF2.gt, alleles = my.alleles, na.string = na.string)


#-------------------------#
# Loading data with R/qtl #
#-------------------------#
cat("\nLoading cross data with R/qtl...\n")

#mydata.qtl <- read.cross("csvs", ".", genfile=myF2.gtf, phefile=myF2.ptf, estimate.map=FALSE, genotypes = c("AA", "AB", "BB")) # Change my alleles
mydata.qtl <- read.cross("csvs", ".", genfile = myF2.gtf, phefile = myF2.ptf, estimate.map = FALSE, genotypes = c("LL", "HL", "HH"), alleles = my.alleles) # Change my alleles

cross.cln <- grep("cross", colnames(mydata.qtl$pheno), ignore.case=TRUE, value=TRUE)

# Check if cross ID can be find in the qtl data
if (length(cross.cln) != 1) {stop("Can't identify the cross column.")}
if (! all(sort(unique(mydata.qtl$pheno[,cross.cln])) == myF2.list)) {
    if (length(myF2.list) > 1) {myF2.msg <- myF2.list} else {myF2.msg <- paste(myF2.list, collapse=" and ")}
    stop(paste("Can't find", myF2.msg,"in the cross column."))
}

# Removing uninformative markers
myuninfo.mrkr <- vector("list", length(myF2.list)) 
names(myuninfo.mrkr) <- myF2.list
for (i in myF2.list) {
    mydata.qtl.tmp     <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
    myuninfo.mrkr[[i]] <- unlist(findDupMarkers(mydata.qtl.tmp, exact.only=FALSE, adjacent.only=TRUE))
}


# Methods used for the QTL analysis
mymethods    <- c("em")
mymethods.nm <- c("EM")
if (mymodel != "np") {
    mymethods    <- c(mymethods, "mr")
    mymethods.nm <- c(mymethods.nm, "marker regression")
}

mycomp.ls <- c(myF2.list,"combination")
myqtl.ls  <- vector("list", length(mycomp.ls))
names(myqtl.ls) <- mycomp.ls

for (i in mycomp.ls) {
    cat(paste(" --Cross",i,"\n"))

    # Subset cross of interest
    if (i != "combination") {
        in.qtl <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
        in.qtl <- drop.markers(in.qtl, myuninfo.mrkr[[i]])
    } else {
        mylist.mkr <- myuninfo.mrkr %>% unlist() %>% unique()
        in.qtl     <- drop.markers(mydata.qtl, mylist.mkr)
    }

    # Genotype probability calculation (needed for some analysis)
    cat("\t -Computing genotype probability...\n")
    in.qtl <- calc.genoprob(in.qtl)
    

    #~~~~~~~~~~~~~~#
    # Scan for QTL #
    #~~~~~~~~~~~~~~#

    out.ls <- vector("list", length(mymethods))
    names(out.ls) <- mymethods

    for (m in mymethods) {
        cat("\t -Performing ", mymethods.nm[mymethods %in% m], " analysis and permutation computation (", my.n.perm, ")...\n", sep="")
        out <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=mymodel)
        out.perm <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=mymodel, n.perm=my.n.perm, n.cluster = n.cluster, verbose=FALSE)
        out.trsh <- as.numeric(sort(out.perm)[round(my.n.perm-my.n.perm*mylod.trsh)])
        if(max(out[,3]) < out.trsh) {warning(m, " method: all LOD scores are under LOD threshold.", immediate.=TRUE, call.=FALSE)}

        # Rename positions
        out[,2] <- unlist(lapply(rownames(out), function(x) rev(strsplit(x,"_")[[1]])[1]))

        # Check for spurious results
        out <- out[ ! is.na(out[,2]),]
        if (any(is.infinite(out[,3]))) {
            warning("Infinite values present. They will be replaced by the maximum value.", immediate.=TRUE, call.=FALSE)
            out[is.infinite(out[,3]),3] <- max(out[!is.infinite(out[,3]),3])
        }

        # Store results
        out.ls[[m]] <- list(lod=out, perm=out.perm, trsh=out.trsh)
    }

    # Creating QTL object
    myqtl.ls[[i]]         <- out.ls
    myqtl.ls[[i]]$genopob <- in.qtl
}

# Save QTL analysis for other scripts to use
save(myqtl.ls, file = paste0(result_fd, "myqtl.ls.RData"))

# QTL identification
mypeaks    <- summary(myqtl.ls[[3]][["mr"]]$lod, perm=myqtl.ls[[3]][["mr"]]$perm, alpha=mylod.trsh)
myqtl.mrkr <- rownames(mypeaks[grep("_[0-9ZW]$", mypeaks[,1], perl=TRUE), ])
myqtl.nb <- length(myqtl.mrkr)


# Genotypes at QTL peaks
mygeno.tb <- pull.geno(myqtl.ls[[3]]$genoprob)
myAF.pheno <- data.frame( 
                "pheno" = pull.pheno(mydata.qtl.combination)[,pheno.cln],
                "AF" = rowSums(mygeno.tb[ , grep(paste(myqtl.mrkr, collapse="|"), colnames(mygeno.tb)) ]) - myqtl.nb   # Normalized number of "alternative" alleles regarding the number of QTLs
                )

# Exporting data
for (i in mycomp.ls) {
    for (j in unique(mymethods)) {
        write.table(get(paste(i,j,sep=".")), paste0("tables/",i,".",j,".tsv"), row.names=FALSE, quote=FALSE, sep="\t")
    }
}





#=========#
# Figures #
#=========#

cat("\nDrawing graphs...\n")

# hist(apply(mydata.gt[,grepl("F2A", colnames(mydata.gt))], 1, function(y) {length(grep("./.",y,fixed=T))}))
#boxplot(mydata.gq[,3:ncol(mydata.gq)])
#boxplot(mydata.dp[,3:ncol(mydata.dp)])

if (dir.exists(graph_fd) == FALSE) {dir.create(graph_fd)}

g.prefix <- paste(g.prefix, colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")


for (i in mycomp.ls) {

    for (m in mymethods) {
    
        myqtl.data <- myqtl.ls[[i]][[m]]
        
        # Determine the max of y axis
        if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
            my.ylim <- ceiling(myqtl.data$trsh)
        } else {
            my.ylim <- ceiling(max(myqtl.data$lod[,3]))
        }

        mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
        mylod <- mylod[! grepl("loc", mylod[,2]), ]

        pdf(paste0(graph_fd,g.prefix,"_",i,"_",m,".pdf"), width=15)
            matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=my.ylim, abline.h=myqtl.data$trsh, data.order=TRUE)
        dev.off()
    }
}


# Phenotype regarding proportion of alleles from QTLs
pdf(paste0(graph_fd,"pheno-AF.pdf"), useDingbats=FALSE)
boxplot(myAF.pheno[,1] ~ myAF.pheno[,2], ylab="Sum of cercariae shed", xlab="Number of high shedder alleles at QTLs", frame=FALSE)
dev.off()


pdf(paste0(graph_fd,"pheno-AF.pdf"), width = 15, useDingbats = FALSE)
plotPXG(mydata.qtl, marker=myqtl.mrkr, pheno.col = pheno.cln)
dev.off()

