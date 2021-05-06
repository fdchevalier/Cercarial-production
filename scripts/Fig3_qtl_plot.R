#!/usr/bin/env Rscript
# Title: Fig3_qtl_plot.R
# Version: 0.2
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-12-22
# Modified in: 2021-05-06



#==========#
# Comments #
#==========#

# Produce figure 3 of the paper showing the QTL plot and the timing of the QTL appearance.



#==========#
# Versions #
#==========#

# v0.2 - 2021-05-06: update path to phenotype files
# v0.1 - 2021-05-06: reshape code
# v0.0 - 2020-12-22: creation



#===================#
# Packages required #
#===================#

cat("Loading packages...\n")

suppressMessages({
    library("qtl")          # For R/qtl commands
    library("magrittr")
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


#-----------------#
# Cross variables #
#-----------------#

mycomp.ls <- "combination"

# mymissing.data must be in [0;1[
mymissing.data <- 0.2

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")
my.geno    <- c("LL", "HL", "HH")

# QTL analysis methods
mymethods <- "em"


#--------------#
# QTL analysis #
#--------------#

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

myshed <- c(1:4)
myQTL <- list("", "Chr_5" ,"Chr_3", c("Chr_1", "Chr_2", "Chr_3"))
mypad <- 0.25

chr.names.cplt <- c("Chr_1", "Chr_2", "Chr_3", "Chr_4", "Chr_5", "Chr_6", "Chr_7", "Chr_W")
# xlab.axis <- c("1", "2", "3", "4", "5", "6", "7", "Z", "Unassembled scaffolds")
xlab.axis <- c("1", "2", "3", "4", "5", "6", "7", "Z")



#=================#
# Data processing #
#=================#

#---------------#
# Sanity checks #
#---------------#

# Check presence of mandatory files
if (! file.exists(paste0(result_fd, "mypheno.qtl.ls.RData"))) { stop(paste0(result_fd, "mypheno.qtl.ls.RData"), " is missing. Exiting...") }
if (! file.exists(paste0(result_fd, "mypheno.qtl.ls.ac.RData"))) { stop(paste0(result_fd, "mypheno.qtl.ls.ac.RData"), " is missing. Exiting...") }


#--------------#
# Data Loading #
#--------------#

cat("Loading data...\n")

# Recreate R/qtl object
mydata.qtl <- read.cross("csvs", ".", genfile = paste0(result_fd, "F2_geno.flt-gt30-dp10.csvs"), phefile = paste0(data_fd, "phenotypes/F2.csv"), estimate.map = FALSE, genotypes = my.geno, alleles = my.alleles)

# Load QTL analysis from singe QTL analysis
load(paste0(result_fd, "mypheno.qtl.ls.RData"))
load(paste0(result_fd, "mypheno.qtl.ls.ac.RData"))



#=========#
# Figures #
#=========#

cat("Generating graphs...\n")

if (dir.exists(graph_fd) == FALSE) {dir.create(graph_fd)}

for (p in 6) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])
    g.prefix <- paste(colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")

    myqtl.ls <- mypheno.qtl.ls[[p]]

    for (j in mycomp.ls) {

        for (m in mymethods) {
        
            myqtl.data <- myqtl.ls[[j]][[m]]
            
            # Determine the max of y axis
            if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
                my.ylim <- ceiling(myqtl.data$trsh)
            } else {
                my.ylim <- ceiling(max(myqtl.data$lod[,3]))
            }

            mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
            mylod <- data.order(mylod)
            mylod <- mylod[! grepl("loc", mylod[,2]), ]
            mylod <- mylod[ mylod[,1] %in% chr.names.cplt, ]
            mylod[,2] <- as.numeric(mylod[,2])

            pdf(paste0(graph_fd,g.prefix,"_",j,"_",m,"_Fig3.pdf"), width=10, height=7)
                layout(matrix(1:2, ncol=1))
                par(mar=c(2,4.5,1,1) + 0.1)

                # LOD score plot
                matplot.data(mylod, 3, datatype="freq", xlab="", ylab="LOD score", ylim.min=0, ylim.max=my.ylim, lwd=2, abline.lwd=2, abline.h=myqtl.data$trsh, data.order=TRUE, by.pos=TRUE)
                mtext(paste0(LETTERS[1], "."), side=3, line=-0.5, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

                # Rectangle plot
                chr.names.true.cplt <- xlab.axis

                chr.names      <- NULL
                chr.names.true <- NULL
                chr.names.all  <- unique(as.vector(mylod[,1]))
                
                for (i in 1:length(chr.names.cplt)) {
                    if (any(grepl(chr.names.cplt[i], chr.names.all))) {
                        chr.names <- c(chr.names, chr.names.cplt[i])
                        chr.names.true <- c(chr.names.true, chr.names.true.cplt[i])
                    }
                }

                # Get the last pb position of each chr
                chr.length <- c(0, sapply(chr.names.all, function (x) tail(mylod[ mylod[,1] == x, 2],1)))
                # Transform the vector to have cumulative positions
                chr.length <- c(0, sapply(2:length(chr.length), function(i) sum(chr.length[1:(i-1)])+chr.length[i]))
                # Get the positions of the chr + scaffold associated
                chr.SNP.length <- c(0, sapply(chr.names, function(x) tail(chr.length[grep(x, names(chr.length))],1)))

                par(mar=c(5,4.5,1,1) + 0.1)
                plot(1, type="n", xlab="Chromosomes", ylab="Shedding week", xlim=c(0, max(chr.length)), ylim=rev(range(c(myshed-mypad, myshed+mypad))), bty = 'n', axes=FALSE)


                for (i in 1:length(myshed)) { 
                    myQTLcol <- rep(NA, length(chr.length[-1]))
                    names(myQTLcol) <- names(chr.length[-1])
                    myQTLcol[myQTL[[i]]] <- "red"
                    rect(chr.length[-length(chr.length)], myshed[i]-mypad, chr.length[-1], myshed[i]+mypad, col=myQTLcol)
                }

                # Axis 1
                mypos <- NULL
                for (i in 1:length(chr.names)) {
                    mypos[i] <- chr.SNP.length[i]+(chr.SNP.length[i+1]-chr.SNP.length[i])/2
                }
                axis(1, at=mypos, labels=xlab.axis, tick=FALSE)

                # Axis 2
                axis(2, at=myshed, labels=myshed, tick=FALSE, las=2)
                
                mtext(paste0(LETTERS[2], "."), side=3, line=-0.5, at=line2user(par("mar")[2],2), cex=par("cex")*2, adj=0)

            dev.off()
        }
    }
}
