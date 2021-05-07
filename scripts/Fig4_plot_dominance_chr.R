#!/usr/bin/env Rscript
# Title: Fig4_plot_dominance_chr.R
# Version: 0.2
# Author: Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2020-03
# Modified in: 2021-05-07



#-------------------
# Loading packages
#-------------------

suppressMessages({
    library("gplots")
    library("plotrix")
    library("magrittr")
    library("multcompView")
})

#----------------
# Loading datasets
#----------------

# Working directory
setwd(file.path(getwd(), "scripts"))

# Folders
graph_fd  <- "../graphs/"
result_fd <- "../results/2-QTL/"

load(paste0(result_fd, "mygeno.pheno.RData"))

# Print the summary of the element of the list
str(mygeno.pheno)

# Extract dataframe from the list object
#chr1 <- mygeno.pheno[[1]]

# Compute mean and standard error of the phenotype for each loci per genotype
mydata <- lapply(mygeno.pheno, function(x) aggregate(. ~ geno, x, function(y) c(mean(y), sd(y)/sqrt(length(y)))))



#-----------
# Variables
#-----------

myorder <- c("LL", "LH", "HH")

# matrix of colors		
mytb <- matrix(c("LL", "darkorange",
				 "LH", "darkgoldenrod", 
				 "HH", "darkolivegreen3",
				 "NA", "white"), ncol=2, byrow=TRUE) 

myab <- c(4,8,12,16)

# matrix of segments (first column: start / second column: end)
myseg <- matrix(c(1,3,
				5,7,
				9,11,
				13,15,
				17,19), ncol=2, byrow=TRUE)
				

#--------------------------
# Create final data table
#--------------------------

# Create a dynamique table with chr, geno, mean, se
mytable <- NULL

for (i in 1:length(mygeno.pheno)) {
	# Add a column to the table with chr names 
	mychr   <- rep(i, nrow(mydata[[i]]))
	mydata.tmp <- cbind(mychr, mydata[[i]])
	
	# Reoder rows
	mydata.tmp[,2] <- as.character(mydata.tmp[,2])
	mydata.tmp     <- mydata.tmp[ order(match(mydata.tmp[,2], myorder)), ]
	
	# Add NA lines to separate the chr (except last row of the final table)
	if (i < length(mygeno.pheno)) { 
		myNA <- rep(NA, ncol(mydata.tmp))
		mydata.tmp <- rbind(mydata.tmp, myNA)
	}
	# Flatten the matrix into columns
	mydata.tmp <- do.call(data.frame, mydata.tmp)
	
	# Concatenate all the tables together
	mytable <- rbind(mytable, mydata.tmp)

}

#---------
# Figures
#---------

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

#pdf(file="Figure 4_plot.pdf", width=9, height=6)
pdf(file=paste0(graph_fd, "Figure 4_barplot.pdf"), width=9, height=6)

mymean <- mytable[,3]
myse <- mytable[,4]
x <- mytable[,2]
mymax <- mymean+myse*1.3

#---Stats---#
# Generating the letters for the posthoc test 

mya.test <- vector("list", length(mygeno.pheno))
mya.letters <- vector("list", length(mygeno.pheno))
for (i in 1:length(mygeno.pheno)) {
    mya.test[[i]] <- pairwise.wilcox.test(mygeno.pheno[[i]]$pheno, mygeno.pheno[[i]]$geno, p.adjust.method="none")
    d <- as.vector(mya.test[[i]]$p.value)
    names(d) <- sapply(colnames(mya.test[[i]]$p.value), function(x) paste0(x, "-", rownames(mya.test[[i]]$p.value))) %>% as.vector()
    d <- d[!is.na(d)]
    mya.letters[[i]] <- multcompLetters(d)
}

#------barplot-----#
mybarplot<- barplot2(mymean, col=mytb[,2], ylim=c(1200,2300), beside=TRUE, ylab="Average cercariae produced", xlab="Parasite genotype", names.arg=mytable[,2] , space=0.5, cex.axis=1, cex.lab=1, cex=1, plot.ci=TRUE, ci.u=mymean+myse, ci.l=mymean-myse, ci.width=0.1, ci.lwd=1.4, lwd=1.4, xpd=FALSE)


text(c(2.5,8.5,14.5,20.5,26.5),2400,labels=c("Chr.1", "Chr.2", "Chr.3", "Chr.4", "Chr.5"),cex=1.5,xpd=TRUE)

##----Plotting stats posthoc test-----#
text(mybarplot,mymax+30,labels=c("a", "a", "b", "", "a", "b", "c", "", "a", "b", "c", "", "a", "a", "b", "", "a", "b", "c"),cex=0.8,xpd=TRUE)



#--------------
# Plot saving
#--------------

dev.off()
