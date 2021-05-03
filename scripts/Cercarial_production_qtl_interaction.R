#!/usr/bin/env Rscript
# Title: Cercarial_production_qtl_interaction.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>, Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2021-05-03
# Modified in:



#==========#
# Comments #
#==========#

# Test QTL interactions. The Cercarial_production_qtl.R must be run at least once before running this script.
# This script takes quite a long time to run.



#==========#
# Versions #
#==========#

# v0.0 - 2021-05-03: creation



#===================#
# Packages required #
#===================#

cat("Loading packages...\n")

suppressMessages({
    library("qtl")          # For R/qtl commands
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

# Function to plot data on exome
source("functions/Sm.matplot.data.R")



#===========#
# Variables #
#===========#

cat("Setting variables...\n")

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"
result_fd <- "../results/2-QTL/"

# Seed for replicating random procedures
myseed <- 1582


#-----------------#
# Cross variables #
#-----------------#

# Allele code
## Be careful to the position of the parents in the vector mypX as the first allele will be given to the first parent. In other words, the first parent needs to correspond to the first allele of my.alleles vector.
my.alleles <- c("L", "H")


#--------------#
# QTL analysis #
#--------------#

nb.markers <- 400 # NULL to include all markers

n.cluster <- 3  # ATTENTION: increasing this number could make your server unusable (make sure you have sufficient RAM)

# Number of permutations to get threshold
my.n.perm <- 1000   # ATTENTION: increasing number will increase computation time substantially
mylod.trsh <- 0.05

# Phenotype to analyze
pheno.cln <- 8  # average



#=================#
# Data processing #
#=================#

#---------------#
# Sanity checks #
#---------------#

# if (! any(grepl(pheno.cln, mypheno.mt[,1]))) {stop("The pheno.cln value is unknown in mypheno.mt table.", call.=FALSE)}

if (! dir.exists(result_fd)) { dir.create(result_fd, recursive = TRUE) }

# Check presence of mandatory files
if (! file.exists(paste0(result_fd, "mypheno.qtl.ls.RData"))) { stop(paste0(result_fd, "mypheno.qtl.ls.RData"), " is missing. Exiting...") }
if (! file.exists(paste0(result_fd, "mypheno.qtl.ls.ac.RData"))) { stop(paste0(result_fd, "mypheno.qtl.ls.ac.RData"), " is missing. Exiting...") }


#--------------#
# Data Loading #
#--------------#

cat("Loading data...\n")

# Load QTL analysis from singe QTL analysis
load(paste0(result_fd, "mypheno.qtl.ls.RData"))
load(paste0(result_fd, "mypheno.qtl.ls.ac.RData"))


#---------------------------#
# Scan for QTL interactions #
#---------------------------#

# Work only on average and combine crosses
myqtl.ls <- mypheno.qtl.ls[["average"]]

#~~~~~~~~~~~~~~~~~~~#
# Subsample markers #
#~~~~~~~~~~~~~~~~~~~#

cat("Subsampling markers...\n")

# This is needed to reduce the computation time for testing interaction

# Chromosome ID
mylod.tb <- myqtl.ls[["combination"]][["em"]]$lod
mychr <- unique(mylod.tb[,1])

# Create an empty list
lod.sub.markers <- vector("list", length(mychr))
names(lod.sub.markers) <- mychr

# Subsample markers proportionally to the number of markers per chromosome
j <- 0
for (i in mychr) {

    # Isolate chromosome and count markers
    chr      <- mylod.tb[mylod.tb[,1] == i, , drop=FALSE]
    chr.mrkr <- nrow(chr)

    # Generate proportion for sampling
    prop <- (nb.markers * chr.mrkr) / nrow(mylod.tb)
    prop <- round(prop)

    # Sample
    set.seed(myseed + j)
    lod.sub.markers[[i]] <- chr[sort(sample(chr.mrkr, prop)),]

    # Increment j
    j <- j + 1
}

# Unlist
myrow.nm <- lapply(lod.sub.markers, row.names) %>% unlist()

if (length(myrow.nm) != nb.markers) {warning(length(myrow.nm), " markers have been subsampled. This is different from the ", nb.markers, " markers requested because the subsampling is porportional to the initial number of markers per chromosome.", call.=FALSE)}

# Drop all but a selected number of markers 
genoprob      <- myqtl.ls[["combination"]][["genoprob"]]
myqtl.reduced <- pull.markers(genoprob, myrow.nm)

#~~~~~~~~~~~~~~~~~~~#
# Test interactions #
#~~~~~~~~~~~~~~~~~~~#

cat("\nTesting interactions. This while take a while...\n")

# Reorder chromosomes (for plotting purpose)
## Select and reorder main chromosomes
mychr.o <- names(myqtl.reduced$geno) %>% grep("_[0-9]$|W$", ., value=TRUE) %>% sort() %>% match(., names(myqtl.reduced$geno))
## Add the others
mychr.o <- c(names(myqtl.reduced$geno)[mychr.o], names(myqtl.reduced$geno)[-mychr.o])
## Reorder
myqtl.reduced$geno <- myqtl.reduced$geno[mychr.o]

# Compute first interaction
myqtl.reduced <- calc.genoprob(myqtl.reduced)
myqtl.i       <- scantwo(myqtl.reduced, pheno.col=pheno.cln)

# Compute permutation thresholds
myperm <- c(1:10)
st.perm  <- vector("list", length(myperm))

for (i in myperm) {
    set.seed(myseed+i)
    st.perm[[i]] <- scantwo(myqtl.reduced, pheno.col=pheno.cln, n.perm=100, n.cluster=10)
}

## combined the permutation files together (c.scantwoperm function)
st.perm <- do.call("c", st.perm)

## Save scantwo and permutation objects as files
save(myqtl.i, file = paste0(result_fd, "scantwo_qtl.i.RData"))
save(st.perm, file = paste0(result_fd, "scantwo_perm.RData"))

## Get the thresholds as a vector
thres <- summary(st.perm) %>% as.data.frame()
colnames(thres) <- (summary(st.perm) %>% attributes())$"names"
thres_vec <- thres[1,]

## QTL interaction/additive effect summary table (complement to heatmap)
cat("\nResult of the interaction test:\n")
summary(myqtl.i, thresholds = as.numeric(thres_vec[1:5]))

# Chromosomes of interest
## The grep command allows for selecting the main scaffolds only
mychr      <- summary(myqtl.i, thresholds = as.numeric(thres_vec[1:5])) %>% as.data.frame() %>% .[,1:2] %>% unlist() %>% unique() %>% grep("_[1-7]$|_ZW$", ., value=TRUE) %>% sort()
myqtl.tb   <- summary(myqtl.ls[["combination"]][["em"]]$lod, perm=myqtl.ls[["combination"]][["em"]]$perm)
myqtl.mrkr <- myqtl.tb[ myqtl.tb[,1] %in% mychr, ]

# Reorder markers
myqtl.mrkr[,1] <- factor(myqtl.mrkr[,1], levels(myqtl.mrkr[,1]) %>% sort())
myqtl.mrkr     <- myqtl.mrkr[ order(myqtl.mrkr[,1]), ]

mygeno.sim <- sim.geno(myqtl.ls[[3]]$genoprob, n.draw=500)

mypos  <- sapply(rownames(myqtl.mrkr), function(x) find.markerpos(mygeno.sim, x)) %>% .[2,] %>% unlist()

myqtl.imp <- makeqtl(mygeno.sim, chr=myqtl.mrkr[,1], pos=mypos)

myqtl.fit <- fitqtl(mygeno.sim, qtl=myqtl.imp, formula=y~Q1+Q2+Q3+Q4+Q5, pheno.col=pheno.cln)

cat("\nResult of the fitting model for additivity:\n")
summary(myqtl.fit)


#------------------#
# Allele dominance #
#------------------#

cat("Analyzing allele dominance...\n")

mygeno.cd <- expand.grid(rep(list(my.alleles), length(my.alleles))) %>% apply(., 1, function(x) x[order(match(x,my.alleles))] %>% paste(., collapse='')) %>% unique()

mygeno.pheno        <- vector("list", nrow(myqtl.mrkr))
names(mygeno.pheno) <- myqtl.mrkr[,1]
for (m in 1:nrow(myqtl.mrkr)) {
    mygeno.pheno[[m]] <- data.frame(geno = pull.geno(myqtl.ls[["combination"]]$genoprob, chr=myqtl.mrkr[m,1])[,rownames(myqtl.mrkr)[m]] %>% sapply(., function(x) mygeno.cd[x]),
                                    pheno = pull.pheno(myqtl.ls[["combination"]]$genoprob, pheno.cln))

	# Reorder factor
	mygeno.pheno[[m]]$geno <- factor(mygeno.pheno[[m]]$geno, mygeno.cd)
}

# Number of genotypes
mygeno.count <- lapply(mygeno.pheno, function(x) data.frame(count = sapply(mygeno.cd, function(y) { z <- x[,1] ; sum(z == y, na.rm=TRUE) }),
                                                            poportion = sapply(mygeno.cd, function(y) { z <- x[,1] ; sum(z == y, na.rm=TRUE) / sum(!is.na(z)) })) )
names(mygeno.count) <- rownames(myqtl.mrkr)

# Test normality of phenotype data
cat("\t- Test normlality of phenotype data:\n")
pull.pheno(myqtl.ls[["combination"]]$genoprob, pheno.cln) %>% shapiro.test() %>% print()

# Test overall differences at each locus
cat("\t- Test overall differences at each locus:\n")
lapply(mygeno.pheno, function(x) kruskal.test(x$pheno ~ x$geno))

# Test differences between genotypes at each locus
cat("\t- Test differences between gentoypes at each locus:\n")
lapply(mygeno.pheno, function(x) pairwise.wilcox.test(x$pheno, x$geno, p.adjust="none")) # p.adjust="bon"

mya.test <- vector("list", length(mygeno.pheno))
mya.letters <- vector("list", length(mygeno.pheno))
for (i in 1:length(mygeno.pheno)) {
    mya.test[[i]] <- pairwise.wilcox.test(mygeno.pheno[[i]]$pheno, mygeno.pheno[[i]]$geno, p.adjust.method="none")
    d <- as.vector(mya.test[[i]]$p.value)
    names(d) <- sapply(colnames(mya.test[[i]]$p.value), function(x) paste0(x, "-", rownames(mya.test[[i]]$p.value))) %>% as.vector()
    d <- d[!is.na(d)]
    mya.letters[[i]] <- multcompLetters(d)
}



#=========#
# Figures #
#=========#

cat("Generating graphs...\n")

if (dir.exists(graph_fd) == FALSE) {dir.create(graph_fd)}

# Interaction heatmap
## Rename chromosomes
mylevels                <- levels(myqtl.i$map$chr) %>% gsub("SM_V7_", "", .)
levels(myqtl.i$map$chr) <- mylevels
## Select chromosomes to plot
chr.o <- mylevels[1:8]
## Plot
pdf(paste0(graph_fd,"scantwo_plot.pdf"), useDingbats = FALSE)
    plot(myqtl.i, chr=chr.o, lower="int", upper="add", col.scheme="viridis")
dev.off()

# Effect plot
pdf(paste0(graph_fd,"effect_plot.pdf"), width = 15, useDingbats = FALSE)
    layout(matrix(1:5, nrow=1))

    # Phenotype against genotype
    for (m in 1:nrow(myqtl.mrkr)) {
        effectplot(myqtl.ls[["combination"]]$genoprob, pheno.col=pheno.cln, rownames(myqtl.mrkr[m,]))
    }
dev.off()
