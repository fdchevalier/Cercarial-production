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

library(qtl)    # For R/qtl commands



#===========#
# Functions #
#===========#

#------------------------#
# Function to save plots #
#------------------------#

save.myplot <- function (filename, width=8, height=8) {
	# Folder creation for graphic output (if does not exist)
	if (file.exists("Graphs/") == FALSE) {dir.create("Graphs")}

	# Graphic saving (png and eps formats)
	dev.copy(png,paste("Graphs/",filename,".png",sep=""), width=width*72, height=height*72, bg="white")
	dev.off()
	dev.copy2pdf(file=paste("Graphs/",filename,".pdf",sep=""), width=width, height=height)
}


#-------------------------#
# Function to save tables #
#-------------------------#

save.mytable <- function (x, filename, ext=".tsv", sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE){
	# Folder creation for table output (if does not exist)
	if (file.exists("Results/") == FALSE) {dir.create("Results")}
	
	# Data export
	write.table(x, file=paste("Results/",filename, ext, sep=""), col.names=col.names, row.names=row.names, append=append, sep=sep)
}


#-------------------#
# Dynamic cross tab #
#-------------------#

mydata.fltr <- function(data.x, data.slct, pattern.data, sign, trsh, missing.data=0.2, na.string=NA) {

    # Usage
    ## data.x: dataframe to filter
    ## data.slct: dataframe used for row selection. Can't be the same as data.
    ## pattern.data: pattern to select column(s) in the dataframe to filter and used for selection
    ## sign: sign for the comparison
    ## trsh: numerical threshold value used for selection
    ## missing.data: proportion of missing data allowed
    ## na.rm: remove any row that contain NA in any selected columns

    # Check steps
    if (missing(data.x) || ! is.data.frame(data.x)) {stop("data.x must be present and must be a dataframe.")}
    if (missing(data.slct) || ! is.data.frame(data.slct)) {stop("data.slct must be present and must be a dataframe.")}
    if (missing(sign) || ! is.character(sign)) {stop("sign must be present and must be a character.")}
    if (missing(trsh) || ! is.numeric(trsh)) {stop("trsh must be present and must be numeric.")}

    # check if same number of row
    # check if same number of column selected
    if (! is.numeric(missing.data) | missing.data < 0 | missing.data >= 1) {stop("missing.data must be numeric and in [0;1[")}    # check if missing data is between 0 and 1

    # Tmp objects
    data.cln <- grep(pattern.data, colnames(data.x))
    data.x.tmp <- data.x[,data.cln]
    data.slct.tmp <- data.slct[,data.cln]

    # Put NA in tmp objects
    data.x.tmp[ ! eval(parse(text=paste("data.slct.tmp",sign,trsh))) ] <- na.string     # Warning: we want the opposite value (!) of the selection to be NA because the expression select the cells that pass the selection criteria

    # Replace values in the original data
    data.x[,data.cln] <- data.x.tmp
    
    # Remove rows with missing data over the threshold
    data.rw.md <- (rowSums(data.x[,data.cln] == na.string)/length(data.cln)) <= missing.data
    #data.x <- data.x[data.rw.md,]  # Keep only rows that pass the missing data test
    data.x[! data.rw.md,data.cln] <- na.string
    return(data.x)
}


#-------------------------------#
# Genotyping data quality check #
#-------------------------------#

quality.check <- function(x, est.rf=FALSE) {
    cat("Segregation distortion.\n")
    gt <- geno.table(x)
    gti <- nrow(gt[gt$P.value < 1e-7,])
    cat(paste("Number of sites showing segregation distortion (chi-sq < 1e-7): ", gti, ". These sites should be removed.", sep=""))

    cat("Pairwise genotype comparisons.\n")
    cg <- comparegeno(x)
    cgi <- nrow(which(cg >0.9, arr.ind=T))/2
    cat(paste("Number of individuals showing 0.9 indentity in their genotypes: ", cgi, ". This may indicate a mix-up or a very very closely related individuals.", sep=""))

    # est.rf takes quite a lot of time so it is why it is off by default
    if(est.rf) {
        cat("Pairwise recombinaison fraction (it will take a while...).\n")
        x <- est.rf(x)
        checkAlleles(x)
    }
}


#--------------------#
# External functions #
#--------------------#

# Function specific to R/qtl
source("~/scripts/0-R_functions_repository/gt2rqtl.R")

# Function to plot data on exome - v2.4
source("~/scripts/0-R_functions_repository/Sm.matplot.data.R_v3.0")



#===========#
# Variables #
#===========#

#--------------#
# Data Loading #
#--------------#

# GT, GQ and DP data file
myfile <- "data/Cerc_prod.ugtpr.snps_indels.raw.wi-baits.nuc.flt-wo_indls"
mydata.gt <- read.csv(paste(myfile,"GT.FORMAT",sep="."), header=TRUE, sep="\t")
mydata.gq <- read.csv(paste(myfile,"GQ.FORMAT",sep="."), header=TRUE, sep="\t", na.strings=".")
mydata.dp <- read.csv(paste(myfile,"DP.FORMAT",sep="."), header=TRUE, sep="\t", na.strings=".")

# Object contaning the input file name
filename <- strsplit(myfile,".csv")

# PT data file
myF2.ptf <- "data/F2.csv"

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
gq.trsh <- NULL
rd.trsh <- 10
# mymis.data must be in [0;1[
mymissing.data <- 0.2

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")


#--------------#
# QTL analysis #
#--------------#

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

if(! any(grepl(pheno.cln, myprefixes[,1]))) {stop("The pheno.cln value is unknown in myprefixes table.", call.=FALSE)}



#=================#
# Data processing #
#=================#

if (mymissing.data < 0 | mymissing.data >= 1) { stop("mymissing.data must be in [0;1[") }


#-----------------------------#
# Data filtering on DP and GQ #
#-----------------------------#
cat("\nFiltering data based on GQ and DP...\n")

# Converting data in numeric
cat(paste("\t- Initial number of alleles before filtering in F2A and F2B:", nrow(mydata.gt), "\n"))

for (i in myF2.list) {

    # Data filtering removing SNPs with low GQ
    if (! is.null(gq.trsh)) {        # gq.trsh is in the section variables
        gq.sign <- ">="
        mydata.gt.tmp.flt.gq <- mydata.fltr(mydata.gt, mydata.gq, i, missing.data=mymissing.data, sign=gq.sign, trsh=gq.trsh)

        nb.alleles <- nrow(mydata.gt.tmp.flt.gq.dp[(rowSums(mydata.gt.tmp.flt.gq.dp[,grep(i,colnames(mydata.gt.tmp.flt.gq.dp))] == "./.")/length(grep(i,colnames(mydata.gt.tmp.flt.gq.dp))) <= mymissing.data),])
        cat(paste("\t- Remaining alleles after GQ filtering in ", i, ":", nb.alleles, "\n", sep=""))
    } else {
        mydata.gt.tmp.flt.gq <- mydata.gt
    }

    # Data filtering removing SNPs with low read.depth
    if (! is.null(rd.trsh)) {
        rd.sign <- ">="
        mydata.gt.tmp.flt.gq.dp <- mydata.fltr(mydata.gt.tmp.flt.gq, mydata.dp, i, missing.data=mymissing.data, sign=rd.sign, trsh=rd.trsh, na.string="./.")

        nb.alleles <- nrow(mydata.gt.tmp.flt.gq.dp[(rowSums(mydata.gt.tmp.flt.gq.dp[,grep(i,colnames(mydata.gt.tmp.flt.gq.dp))] == "./.")/length(grep(i,colnames(mydata.gt.tmp.flt.gq.dp))) <= mymissing.data),])
        cat(paste("\t- Remaining alleles after DP filtering in ", i,": ", nb.alleles, "\n", sep=""))
    } else {
        mydata.gt.tmp.flt.gq.dp <- mydata.gt
    }

#    # Rename data
#    assign(paste("mydata.gt.F2",sep=""), mydata.gt.tmp.flt.gq.dp)
#    assign(paste("my",i,".gt",sep=""), name.gt)

    mydata.gt <- mydata.gt.tmp.flt.gq.dp
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
myF2.gt <- paste("tables/", myF2.gt, myflt.suffix, sep="")
myF2.gtf <- paste(myF2.gt, "csvs", sep=".")


# Column associated with name
mypA <- sapply(mypA, function(x) grep(x, colnames(mydata.gt)))
myF1A <- grep(myF1A, colnames(mydata.gt))
myF2A <- grep(myF2A, colnames(mydata.gt))

mypB <- sapply(mypB, function(x) grep(x, colnames(mydata.gt)))
myF1B <- grep(myF1B, colnames(mydata.gt))
myF2B <- grep(myF2B, colnames(mydata.gt))

myp <- list(mypA, mypB)
myF1 <- list(myF1A, myF1B)
myF2 <- list(myF2A, myF2B)


# Table conversion
if (file.exists("tables/") == FALSE) {dir.create("tables")}
gt2rqtl(mydata.gt, parents.cln=myp, F1.cln=myF1, F2.cln=myF2, out.fmt="csvs", out.name=myF2.gt, simplify.unass.name=TRUE, alleles=my.alleles)


#-------------------------#
# Loading data with R/qtl #
#-------------------------#
cat("\nLoading cross data with R/qtl...\n")

#mydata.qtl <- read.cross("csvs", ".", genfile=myF2.gtf, phefile=myF2.ptf, estimate.map=FALSE, genotypes = c("AA", "AB", "BB")) # Change my alleles
mydata.qtl <- read.cross("csvs", ".", genfile=myF2.gtf, phefile=myF2.ptf, estimate.map=FALSE, genotypes=c("LL", "HL", "HH"), alleles=my.alleles) # Change my alleles

cross.cln <- grep("cross", colnames(mydata.qtl$pheno), ignore.case=TRUE, value=TRUE)

# Check if cross ID can be find in the qtl data
if (length(cross.cln) != 1) {stop("Can't identify the cross column.")}
if (! all(sort(unique(mydata.qtl$pheno[,cross.cln])) == myF2.list)) {
    if (length(myF2.list) > 1) {myF2.msg <- myF2.list} else {myF2.msg <- paste(myF2.list, collapse=" and ")}
    stop(paste("Can't find", myF2.msg,"in the cross column."))
}

# Removing uninformative markers
for (i in myF2.list) {
    mydata.qtl.tmp <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
#    mymap.tmp <- est.map(mydata.qtl.tmp)
#
#    mylist.tmp <- unlist(lapply(mymap.tmp, function(x) names(x)[! duplicated(round(x))]))
    mylist.tmp <- unlist(findDupMarkers(mydata.qtl.tmp, exact.only=FALSE, adjacent.only=TRUE))
    assign(paste0(i,".inf.mrkr"), mylist.tmp)
}

mymethods <- c()

for (i in c(myF2.list,"combination")) {
    cat(paste(" --Cross",i,"\n"))

    # Subset cross of interest
    if (i != "combination") {
        in.qtl <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
#        in.qtl <- pull.markers(in.qtl, get(paste0(i,".inf.mrkr")))
        in.qtl <- drop.markers(in.qtl, get(paste0(i,".inf.mrkr")))
    } else {
        mylist.mkr <- mydata.gt[,1]
        for (j in myF2.list) {
            if (length(get(paste0(j,".inf.mrkr"))) < length(mylist.mkr)) {
                mylist.mkr <- get(paste0(j,".inf.mrkr"))
            }
    }
        in.qtl <- mydata.qtl
#        in.qtl <- pull.markers(in.qtl, mylist.mkr)
        in.qtl <- drop.markers(in.qtl, mylist.mkr)
    }

    # Genotype probability calculation (needed for some analysis)
    cat("\t -Computing genotype probability...\n")
    in.qtl <- calc.genoprob(in.qtl)
    
    myobjects <- c("myobjects")

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Scan for QTL using Marker Regression #
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    if (mymodel != "np") {
        cat(paste("\t -Performing marker regresion analysis and permutation computation (", my.n.perm, ")...\n", sep=""))
        out.mr <- scanone(in.qtl, pheno.col=pheno.cln, method="mr", model=mymodel)
        out.mr.perm <- scanone(in.qtl, pheno.col=pheno.cln, method="mr", model=mymodel, n.perm=my.n.perm, verbose=FALSE)
        out.mr.trsh <- as.numeric(sort(out.mr.perm)[round(my.n.perm-my.n.perm*mylod.trsh)])
        if(max(out.mr[,3]) < out.mr.trsh) {warning("MR method: all LOD scores are under LOD threshold.", immediate.=TRUE, call.=FALSE)}

        # Rename positions
        #out.mr[,2] <- unlist(lapply(rownames(out.mr), function(x) strsplit(x,"_")[[1]][5]))
        out.mr[,2] <- unlist(lapply(rownames(out.mr), function(x) rev(strsplit(x,"_")[[1]])[1]))
        
        # Check for spurious results
        if (any(is.infinite(out.mr[,3]))) {
            warning("Infinite values present. They will be replaced by the maximum value.", immediate.=TRUE, call.=FALSE)
            out.mr[is.infinite(out.mr[,3]),3] <- max(out.mr[!is.infinite(out.mr[,3]),3])
        }

        # Store methods
        mymethods <- c(mymethods, "mr")

        # Create object
        assign(paste(i,".mr",sep=""), out.mr)
        assign(paste(i,".mr.trsh",sep=""), out.mr.trsh)
        myobjects <- c(myobjects, "out.mr", "out.mr.trsh")
    }


    #~~~~~~~~~~~~~~~~~~~~~~~#
    # Scan for QTL using EM #
    #~~~~~~~~~~~~~~~~~~~~~~~#

    cat(paste("\t -Performing EM analysis and permutation computation (", my.n.perm, ")...\n", sep=""))
    out.em <- scanone(in.qtl, pheno.col=pheno.cln, method="em", model=mymodel)
    out.em.perm <- scanone(in.qtl, pheno.col=pheno.cln, method="em", model=mymodel, n.perm=my.n.perm, verbose=FALSE)
    out.em.trsh <- as.numeric(sort(out.em.perm)[round(my.n.perm-my.n.perm*mylod.trsh)])
    if(max(out.em[,3]) < out.em.trsh) {warning("EM method: all LOD scores are under LOD threshold.", immediate.=TRUE, call.=FALSE)}

    # Rename positions
    #out.em[,2] <- unlist(lapply(rownames(out.em), function(x) strsplit(x,"_")[[1]][5]))
    out.em[,2] <- unlist(lapply(rownames(out.em), function(x) rev(strsplit(x,"_")[[1]])[1]))

    # Check for spurious results
    out.em <- out.em[ ! is.na(out.em[,2]),]
    if (any(is.infinite(out.em[,3]))) {
        warning("Infinite values present. They will be replaced by the maximum value.", immediate.=TRUE, call.=FALSE)
        out.em[is.infinite(out.em[,3]),3] <- max(out.em[!is.infinite(out.em[,3]),3])
    }

    # Store methods
    mymethods <- c(mymethods, "em")

    # Create object
    assign(paste(i,".em",sep=""), out.em)
    assign(paste(i,".em.trsh",sep=""), out.em.trsh)
    myobjects <- c(myobjects, "out.em", "out.em.trsh")

    # Creating QTL object
    assign(paste("mydata.qtl.",i,sep=""), in.qtl)
    myobjects <- c(myobjects, "in.qtl")

    # Remove tmp objects
    rm(list=myobjects)
}


# QTL identification
ce.ordered <- combination.em[ order(combination.em[,3], decreasing=T), ]
myqtl.mrkr <- rownames( ce.ordered[! duplicated(ce.ordered[1]) & ce.ordered[,3] > combination.em.trsh,] )
myqtl.nb <- length(myqtl.mrkr)

# Genotypes at QTL peaks
mygeno.tb <- pull.geno(mydata.qtl.combination)
myAF.pheno <- data.frame( 
                "pheno" = pull.pheno(mydata.qtl.combination)[,pheno.cln],
                "AF" = rowSums(mygeno.tb[ , grep(paste(myqtl.mrkr, collapse="|"), colnames(mygeno.tb)) ]) - myqtl.nb   # Normalized number of "alternative" alleles regarding the number of QTLs
                )

# Exporting data
for (i in c(myF2.list,"combination")) {
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

if (file.exists("graphs/") == FALSE) {dir.create("graphs")}

g.prefix <- paste(g.prefix, colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste("md-",mymissing.data, sep=""), sep=".")


for (i in c(myF2.list,"combination")) {
    for (j in unique(mymethods)) {

        # Determine the max of y axis
        if (get(paste(i,j,"trsh",sep=".")) > max(get(paste(i,j,sep="."))[,3])) {
            my.ylim <- ceiling(get(paste(i,j,"trsh",sep=".")))
        } else {
            my.ylim <- ceiling(max(get(paste(i,j,sep="."))[,3]))
        }

        pdf(paste("graphs/",g.prefix,"_",i,"_",j,".pdf",sep=""), width=15)
            matplot.data(get(paste(i,j,sep=".")), 3, datatype="freq", ylab="LOD score", ylim.max=my.ylim, abline.h=get(paste(i,j,"trsh",sep=".")), data.order=FALSE)
        dev.off()
    }
}


# Phenotype regarding proportion of alleles from QTLs
pdf(paste0("graphs/","pheno-AF.pdf"), useDingbats=FALSE)
    boxplot(myAF.pheno[,1] ~ myAF.pheno[,2], ylab="Sum of cercariae shed", xlab="Number of high shedder alleles at QTLs", frame=FALSE)
dev.off()

