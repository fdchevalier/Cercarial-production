#!/usr/bin/env Rscript
# Title: Cercarial_production_qtl.R
# Version: 0.6
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>, Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2016-02-18
# Modified in: 2021-05-06



#==========#
# Comments #
#==========#

# Test the presence of QTL related to production of cercariae in Schistosoma mansoni parasites.



#==========#
# Versions #
#==========#

# v0.6 - 2021-05-06: update path to phenotype files / add genotype variable
# v0.5 - 2021-05-05: add line2user function
# v0.4 - 2021-05-03: move QTL interaction test in a new script / add missing function / correct typos
# v0.3 - 2021-05-03: script renamed and code reshaped
# v0.2 - 2016-07-15: missing data bug correction
# v0.1 - 2016-05-18: combined analysis included
# v0.0 - 2016-02-18: creation



#===================#
# Packages required #
#===================#

cat("Loading packages...\n")

suppressMessages({
    library("qtl")          # For R/qtl commands
    library("vcfR")
    library("magrittr")
    library("snow")         # For parallel permutations
    library("multcompView")	# For multcompLetters
})



#===========#
# Functions #
#===========#

# Working directory
setwd(file.path(getwd(), "scripts"))

#--------------------#
# External functions #
#--------------------#

source("functions/Cercarial_production_qtl_func.R")

# Function specific to R/qtl
source("functions/gt2rqtl.R")

# Function to rename chromosomes
source("functions/rename_chr.R")

# Function to plot data on exome
source("functions/Sm.matplot.data.R")
source("functions/line2user.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/2-QTL/"

myvcf_f <- paste0(data_fd, "calling/cerc_prod_snpEff_reduced.vcf.gz")

# Numeric format of the different data type
myfmt <- matrix(c(
            "GT", FALSE,
            "GQ", TRUE,
            "DP", TRUE),
            ncol=2, byrow=TRUE)


# Object containing the input file name
filename <- basename(myvcf_f) %>% strsplit(., ".vcf.gz")

# PT data file
myF2.ptf <- paste0(data_fd, "phenotypes/F2.csv")

# PT data file
myF2.sex.f <- paste0(data_fd, "phenotypes/sex.tsv")

# Output name of the R/qtl GT table without any filename extension
myF2.gt <- "F2_geno"

# Seed for replicating random procedures
myseed <- 1582


#-----------------#
# Cross variables #
#-----------------#

# F2A variables
mypA  <- c("SmBRE4_m", "SmLE19_f")   # Be careful at the order of the parents (cf. Allele code below)
myF1A <- c("F1A")
myF2A <- c("F2A")

#mypB <- c("SmLE15_m", "SmBRE2_f")
mypB  <- c("SmBRE2_f", "SmLE15_m")   # Be careful at the order of the parents (cf. Allele code below)
myF1B <- c("F1B")
myF2B <- c("F2B")

myF2.list <- c("F2A", "F2B")

# GQ and read depth (rd) and missing data thresholds
## Assign NULL to skip to the variable to skip the corresponding filtering step
gq.trsh <- 30
rd.trsh <- 10
## mymissing.data must be in [0;1[
mymissing.data <- 0.2

na.string <- NA # Could be "./."

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")
my.geno    <- c("LL", "HL", "HH")


#--------------#
# QTL analysis #
#--------------#

nb.markers <- 400 # NULL to include all markers

n.cluster <- 3  # ATTENTION: increasing this number could make your server unusable (make sure you have sufficient RAM)

# Number of permutations to get threshold
my.n.perm <- 1000   # ATTENTION: increasing number will increase computation time substantially
mylod.trsh <- 0.05

# # Model for QTL scan
# pheno.cln <- 8

# Prefix for graph
mypheno.mt <- matrix(c(
    #pheno.cln, graph_prefix, scan_model
        3,      "d1",         "normal",
        4,      "d2",         "normal",
        5,      "d3",         "normal",
        6,      "d4",         "normal",
        7,      "sum",        "normal",
        8,      "average",    "normal",
        9,      "PO",         "normal",
        10,     "Hb",         "normal"), ncol=3, byrow=TRUE)       # The pheno.cln cannot start at one because it is the ID column

# g.prefix <- mypheno.mt[ mypheno.mt[,1] == pheno.cln, 2 ]
# mymodel <- mypheno.mt[ mypheno.mt[,1] == pheno.cln, 3]



#=================#
# Data processing #
#=================#

#---------------#
# Sanity checks #
#---------------#

if (mymissing.data < 0 | mymissing.data >= 1) { stop("mymissing.data must be in [0;1[") }
# if (! any(grepl(pheno.cln, mypheno.mt[,1]))) {stop("The pheno.cln value is unknown in mypheno.mt table.", call.=FALSE)}

if (! dir.exists(result_fd)) { dir.create(result_fd, recursive = TRUE) }


#--------------#
# Data Loading #
#--------------#

cat("Loading data...\n")

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

cat("Filtering data based on GQ and DP...\n")

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
myF2.gt <- paste0(result_fd, myF2.gt, myflt.suffix)


# Column associated with name
mypA  <- sapply(mypA, function(x) grep(x, colnames(mydata$GT)))
myF1A <- grep(myF1A, colnames(mydata$GT))
myF2A <- grep(myF2A, colnames(mydata$GT))

mypB  <- sapply(mypB, function(x) grep(x, colnames(mydata$GT)))
myF1B <- grep(myF1B, colnames(mydata$GT))
myF2B <- grep(myF2B, colnames(mydata$GT))

myp  <- list(mypA, mypB)
myF1 <- list(myF1A, myF1B)
myF2 <- list(myF2A, myF2B)


# Table conversion
gt2rqtl(mydata$GT, parents.cln = myp, F1.cln = myF1, F2.cln = myF2, out.fmt = "csvs", out.name = myF2.gt, alleles = my.alleles, na.string = na.string, add.distance = TRUE)


#-------------------------#
# Loading data with R/qtl #
#-------------------------#

cat("\nLoading cross data with R/qtl...\n")

mydata.qtl <- read.cross("csvs", ".", genfile = paste0(myF2.gt, ".csvs"), phefile = myF2.ptf, estimate.map = FALSE, genotypes = my.geno, alleles = my.alleles)

cross.cln <- grep("cross", colnames(mydata.qtl$pheno), ignore.case=TRUE, value=TRUE)
    
# Check if cross ID can be find in the qtl data
if (length(cross.cln) != 1) {stop("Can't identify the cross column.")}
if (! all(sort(unique(mydata.qtl$pheno[,cross.cln])) == myF2.list)) {
    if (length(myF2.list) > 1) {myF2.msg <- myF2.list} else {myF2.msg <- paste(myF2.list, collapse=" and ")}
    stop(paste("Can't find", myF2.msg,"in the cross column."))
}

# Identifying uninformative markers
myuninfo.mrkr <- vector("list", length(myF2.list)) 
names(myuninfo.mrkr) <- myF2.list
for (i in myF2.list) {
    mydata.qtl.tmp     <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
    myuninfo.mrkr[[i]] <- unlist(findDupMarkers(mydata.qtl.tmp, exact.only=FALSE, adjacent.only=TRUE))
}

# List to store QTL and related information
mycomp.ls <- c(myF2.list,"combination")

# Sex covariate
mysex <- read.table(myF2.sex.f, header = TRUE, stringsAsFactors = FALSE)
mysex <- mysex[grepl("F2", mysex[,1]), ]
mysex[is.na(mysex[, 3]), 3] <- 0.75
mysex[mysex[, 3] > 0.7 & mysex[, 3] < 0.9, 4] <- NA
mysex <- mysex[, c(1,4)]
mysex[! is.na(mysex[, 2]) & mysex[, 2 ] == "female", 2] <- 0    # Coding sex as recommended per the book
mysex[! is.na(mysex[, 2]) & mysex[, 2 ] == "male", 2]   <- 1    # Coding sex as recommended per the book
mysex[, 2] <- as.numeric(mysex[, 2])
mysex <- merge(pull.pheno(mydata.qtl), mysex, by = 1, all = TRUE)
mysex <- mysex[, ncol(mysex), drop = FALSE]

# Lists to store results
mypheno.qtl.ls  <- vector("list", nrow(mypheno.mt)) ## For simple scanone
names(mypheno.qtl.ls) <- mypheno.mt[,2]
mypheno.qtl.ls.ac <- mypheno.qtl.ls ## For scanone with covariate

# Analyze all phenotypes
for (p in 1:nrow(mypheno.mt)) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])

    # Phenotype name
    pheno.nm <- pull.pheno(mydata.qtl) %>% colnames() %>% .[3]

    # Methods used for the QTL analysis
    mymethods    <- c("em")
    mymethods.nm <- c("EM")
    if (mymodel != "np") {
        mymethods    <- c(mymethods, "mr")
        mymethods.nm <- c(mymethods.nm, "marker regression")
    }

    # QTL list
    myqtl.ls  <- vector("list", length(mycomp.ls))
    names(myqtl.ls) <- mycomp.ls

    myqtl.ls.ac <- myqtl.ls

    for (i in mycomp.ls) {
        cat(" --Cross ", i, " (", pheno.nm, ")\n", sep="")

        # Subset cross of interest
        if (i != "combination") {
            in.qtl <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
            in.qtl <- drop.markers(in.qtl, myuninfo.mrkr[[i]])
            mysex.tmp <- mysex[ pull.pheno(mydata.qtl)[, 2] == i, , drop = FALSE ]
        } else {
            mylist.mkr <- myuninfo.mrkr %>% unlist() %>% unique()
            in.qtl     <- drop.markers(mydata.qtl, mylist.mkr)
            mysex.tmp  <- mysex
        }

        # Genotype probability calculation (needed for some analysis)
        cat("\t -Computing genotype probability...\n")
        in.qtl <- calc.genoprob(in.qtl)
        

        #~~~~~~~~~~~~~~#
        # Scan for QTL #
        #~~~~~~~~~~~~~~#

        myqtl.ls[[i]]          <- myscanone(in.qtl, pheno.cln, mymethods, mymodel, my.n.perm, n.cluster)
        myqtl.ls[[i]]$genoprob <- in.qtl
        
        myqtl.ls.ac[[i]]          <- myscanone(in.qtl, pheno.cln, mymethods, mymodel, my.n.perm, n.cluster, addcovar=mysex.tmp)
        myqtl.ls.ac[[i]]$genoprob <- in.qtl
    }

    # Store in the designated slot
    mypheno.qtl.ls[[mypheno.mt[p,2]]]    <- myqtl.ls
    mypheno.qtl.ls.ac[[mypheno.mt[p,2]]] <- myqtl.ls.ac

}

# Save QTL analysis for other scripts to use
save(mypheno.qtl.ls, file = paste0(result_fd, "mypheno.qtl.ls.RData"))
save(mypheno.qtl.ls.ac, file = paste0(result_fd, "mypheno.qtl.ls.ac.RData"))

# QTL identification
cat("\nSummary of the scanone analysis (EM method):\n")
for (p in 1:nrow(mypheno.mt)) {
    summary(mypheno.qtl.ls[[p]][[3]][["em"]]$lod, perm=mypheno.qtl.ls[[p]][[3]][["em"]]$perm, alpha=mylod.trsh, pvalues=TRUE) %>% print()
}



#=========#
# Figures #
#=========#

cat("\nGenerating graphs...\n")

if (dir.exists(graph_fd) == FALSE) {dir.create(graph_fd)}

# Analyze all phenotypes
for (p in 1:nrow(mypheno.mt)) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])
    # g.prefix <- mypheno.mt[p, 2 ]
    # g.prefix <- paste(g.prefix, colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")
    g.prefix <- paste(colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")

    myqtl.ls <- mypheno.qtl.ls[[p]]

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
            mylod[,2] <- as.numeric(mylod[,2])

            pdf(paste0(graph_fd,g.prefix,"_",i,"_",m,".pdf"), width=15)
                matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=my.ylim, abline.h=myqtl.data$trsh, data.order=TRUE, by.pos=TRUE)
            dev.off()
        }
    }
}


# Superimposed LOD scores for cercarial shedding time points
p.ls  <- c(3:6,8)-2
# myclr <- viridisLite::viridis(length(p.ls), end=0.8, direction=-1)
myclr <- c(rgb(1, 0, 0, alpha=0.5), gray.colors(length(p.ls) - 1, start=0.05, end=0.75, gamm=1)) %>% rev()
# myclr <- viridisLite::magma(length(p.ls))
# mylwd <- seq(1, 6, length.out=length(p.ls))
mylwd <- c(rep(1, length(p.ls) - 1), 3)
# mylty <- 5:1
for (i in mycomp.ls) {
    for (m in mymethods) {
		
		g.prefix <- paste(i, m, paste0("md-", mymissing.data), sep=".")

		mymax <- lapply(mypheno.qtl.ls[p.ls], function(x) unlist(x[[i]][[m]][["lod"]][,3]) %>% max()) %>% unlist() %>% max()

        pdf(paste0(graph_fd,g.prefix,".pdf"), width=15)
	
        for (p in rev(p.ls)) {
    
			# Phenotype column
			pheno.cln <- as.numeric(mypheno.mt[p, 1])

			myqtl.ls <- mypheno.qtl.ls[[p]]
			
			myqtl.data <- myqtl.ls[[i]][[m]]
            
            # Determine the max of y axis
            if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
                my.ylim <- ceiling(myqtl.data$trsh)
            } else {
                my.ylim <- ceiling(max(myqtl.data$lod[,3]))
            }

            mylod     <- rename_chr_SmV7(myqtl.data$lod, 1) %>% as.data.frame()
            mylod     <- mylod[! grepl("loc", mylod[,2]), ]
            mylod[,2] <- as.numeric(mylod[,2])

            p.i <- match(p, p.ls)
            if (p.i == length(p.ls)) {
                myadd <- FALSE
            } else {
                myadd <- TRUE
            }

            matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=mymax, col=rep(myclr[p.i], 2), data.order=TRUE, by.pos=TRUE, lwd=mylwd[p.i], chr.bg="white", add=myadd)
        }

        legend("topright", legend=c("Week 4 PI", "Week 5 PI", "Week 6 PI", "Week 7 PI", "Average"), col=myclr, bty="n", lty=1, lwd=mylwd, xjust=0) # Need work on the adjustment
        
        dev.off()
    }
}


# Separate LOD scores for separate shedding time points
myrow <- 4
pdf(paste0(graph_fd,"Supp. Fig. 2 - Shedding_panel_weeks.pdf"), width=10, height=4*myrow)

par(mar=c(5,4.5,3,1.5) + 0.1)

layout(matrix(1:myrow, ncol=1))

for (p in 1:nrow(mypheno.mt[1:4,])) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])

    myqtl.ls <- mypheno.qtl.ls[[p]]

    # for (i in mycomp.ls) {
    for (i in "combination") {
		m <- "em"

        my.ylim <- lapply(mypheno.qtl.ls[p.ls], function(x) unlist(x[[i]][[m]][["lod"]][,3]) %>% max()) %>% unlist() %>% max() %>% signif(., 2)

        # for (m in mymethods) {
        for (m in m) {
        
            myqtl.data <- myqtl.ls[[i]][[m]]

            mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
            mylod <- mylod[! grepl("loc", mylod[,2]), ]
            mylod[,2] <- as.numeric(mylod[,2])

            matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.min=0, ylim.max=my.ylim, abline.h=myqtl.data$trsh, data.order=TRUE, by.pos=TRUE)
            mtext(paste0(LETTERS[p], "."), side=3, line=0.75, at=line2user(par("mar")[2],2), cex=par("cex")*2.5, adj=0)
            title(paste("Shedding", p), cex.main=par("cex")*3)
        }
    }
}
dev.off()


# LOD scores for each cross
myrow <- 2
pdf(paste0(graph_fd,"Supp. Fig. 3 - LOD score by cross.pdf"), width=10, height=4*myrow)

par(mar=c(5,4.5,3,1.5) + 0.1)

layout(matrix(1:myrow, ncol=1))

myqtl.ls <- mypheno.qtl.ls[["average"]]

m <- "em"

for (i in myF2.list) {

	myletter <- LETTERS[match(i, myF2.list)]

	my.ylim <- lapply(mypheno.qtl.ls[p.ls], function(x) unlist(x[[i]][[m]][["lod"]][,3]) %>% max()) %>% unlist() %>% max() %>% signif(., 2)

	myqtl.data <- myqtl.ls[[i]][[m]]

	mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
	mylod <- mylod[! grepl("loc", mylod[,2]), ]
	mylod[,2] <- as.numeric(mylod[,2])

	matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.min=0, ylim.max=my.ylim, data.order=TRUE, by.pos=TRUE, cex.axis=0.75)
	mtext(paste0(myletter, "."), side=3, line=0.75, at=line2user(par("mar")[2],2), cex=par("cex")*2.5, adj=0)
	title(paste("Cross", myletter), cex.main=par("cex")*2)
}

dev.off()
