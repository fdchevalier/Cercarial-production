#!/usr/bin/env Rscript
# Title: Supp_Fig1_cercarial_shedding_sex.R
# Version: 0.2
# Author: Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2020-06
# Modified in: 2021-05-05



#-------------------
# Loading packages
#-------------------

suppressMessages({
    library("gplots")
    library("plotrix")
    library("magrittr")
})

#-----------------------------
# Loading dataset
#-----------------------------

# Working directory
setwd(file.path(getwd(), "scripts"))

# Functions
source("functions/line2user.R")

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"

mydataF1 <- read.table(paste0(data_fd, "phenotyping/F1.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")

#---------F1 progeny----------#
mydataA <- mydataF1[mydataF1[,2] == "F1A",]
mydataB <- mydataF1[mydataF1[,2] == "F1B",]

dataF1Am <- mydataA[mydataA[,11] == "m",]
dataF1Af <- mydataA[mydataA[,11] == "f",]

dataF1Bm <- mydataB[mydataB[,11] == "m",]
dataF1Bf <- mydataB[mydataB[,11] == "f",]

#---------F2 progeny----------#

pheno <- read.csv(paste0(data_fd, "phenotyping/F2.csv"))

b <- read.table(paste0(data_fd, "phenotyping/sex.tsv"), header = TRUE, stringsAsFactors = FALSE)
names(b) <- c("id", "read_depth", "ratio", "sex")

c <- merge(pheno, b ,by = 1)

# Mask uncertain sexing
c <- c[! is.na(c[, "sex"]), ]
c <- c[! (c[,"ratio"] > 0.7 & c[,"ratio"] < 0.9), ]

cA <- c[ c[,2] == "F2A", ]
cB <- c[ c[,2] == "F2B", ]

dataF2Am <- cA[cA[,13] == "male",]
dataF2Af <- cA[cA[,13] == "female",]

dataF2Bm <- cB[cB[,13] == "male",]
dataF2Bf <- cB[cB[,13] == "female",]


#----------
# Variables
#----------


F1_2A <- c("gray50")
F1_2B <- c("white")
myborders <- c("royalblue1", "deeppink3")
mynames <- c("male","female")

c.type <- c("m","f")
c.typeF2 <- c("male","female")

#mycex.text <- 1
mycex.axis <- 1
mycex.text <- 1
my.x.pp.text <- 1.9
shifts <- c(-0.2,0.2)

#---------
# Figures
#---------

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

pdf(file=paste0(graph_fd, "Supplementary_Figure1.pdf"), width=6, height=6)

layout(matrix(c(1,2,3,4),2,2))

#------------------#
#----Graph F1A-----#

par(mar=c(3,5,4,2))
boxplot(dataF1Am[,8], dataF1Af[,8], ylab="", ylim=c(0,3000), boxwex=0.6, cex.axis=1, main="F1A", col=F1_2A, border=myborders, names=mynames, frame=FALSE, lwd=2, las=1, outline=FALSE)

mtext(side=2, text="Average cercariae produced", line=3.5, cex=0.9)


#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.type)) {   
       inds <- which(mydataA$sex == c.type[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(mydataA$average[inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.1, "**", xpd=TRUE, cex=1.5)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.03, 1.90)
}

#---Add letters on the figure panel---#
mtext("A.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------#
#----Graph F2A-----#

par(mar=c(3,5,4,2))
boxplot(dataF2Am[,8], dataF2Af[,8], ylab="", ylim=c(0,4000),boxwex=0.6, cex.axis=1, main="F2A", col=F1_2A, border=myborders, names=mynames, las=1, frame=FALSE, lwd=2, outline=FALSE)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.typeF2)) {   
       inds <- which(cA$sex == c.typeF2[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(cA$average[inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.1, "NS", xpd=TRUE, cex=1)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.03, 1.90)
}

mtext(side=2, text="Average cercariae produced", line=3.5, cex=0.9)

#---Add letters on the figure panel---#
mtext("C.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------#
#----Graph F1B-----#

par(mar=c(3,3,4,2))
boxplot(dataF1Bm[,8], dataF1Bf[,8], ylab="", ylim=c(0,3000), boxwex=0.6, cex.axis=1, main="F1B", col=F1_2B, border=myborders, names=mynames, las=1, frame=FALSE, lwd=2, outline=FALSE)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.type)) {   
       inds <- which(mydataB$sex == c.type[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(mydataB$average[inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.15, "NS", xpd=TRUE, cex=1)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.035, 1.90)
}

#---Add letters on the figure panel---#
mtext("B.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------#
#----Graph F2B-----#

par(mar=c(3,3,4,2))
boxplot(dataF2Bm[,8], dataF2Bf[,8], ylab="", ylim=c(0,4000), boxwex=0.6, cex.axis=1, main="F2B", col=F1_2B, border=myborders, names=mynames, las=1, frame=FALSE, lwd=2, outline=FALSE)

#---Add statistical analysis on the plot---#

y.maxs <- c()  
    for (i in 1:length(c.typeF2)) {   
       inds <- which(cB$sex == c.typeF2[i]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(cB$average[inds], plot=FALSE)[[1]][5])
     
}

y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1.5, y.maxs.grp*1.1, "*", xpd=TRUE, cex=1.5)

for (t.id in 1:length(c.type)) {
	segments(1.10, y.maxs.grp[t.id]*1.035, 1.90)
}

#---Add letters on the figure panel---#
mtext("D.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#-----------------------------------
# Plot saving
#-----------------------------------

dev.off()
