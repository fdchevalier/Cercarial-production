#!/usr/bin/env Rscript
# Title: Castle_Wright_estimator_script.R
# Version: 0.1
# Author: Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2020-10-30
# Modified in: 2021-01-12



#-------------------
# Comments
#-------------------

# Compute the Castle-Wright estimator to evaluate the number of loci involved


#-------------------
# Versions
#-------------------

# v0.1 - 2021-01-12: clean code / integrate within local environment
# v0.0 - 2020-10-30: creation


#-------------------
# Loading packages
#-------------------

suppressMessages({
    library("magrittr")
})


#-----------------------------
# Loading data
#-----------------------------

# Working directory
setwd(file.path(getwd(), "scripts"))

mydataF0 <- read.csv("../data/phenotyping/F0_parental_populations.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF1 <- read.csv("../data/phenotyping/F1.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF2 <- read.csv("../data/phenotyping/F2.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")


#-----------------------------
# Castle-Wright estimator Ne
#-----------------------------
# Combine the tables F1/F2

mydataF1F2 <- rbind(mydataF1[,1:8], mydataF2[,1:8])

# Make a list with the id of each parents of the crosses
## Be careful to the order of the parents in the list
F0 <- list("A" = c("SmLE_15","SmBRE_4"), "B" = c("SmLE_19","SmBRE_2"))
mycomp <- c("A", "B", "combined") 

# Ne estimated for Cross A and B (independantly)
#------------------------------------------------

av_F0 <- NULL
var_F <- NULL
Ne    <- NULL


for (i in mycomp) {
	# if "A" slot of the list take F1A or F2A, if "B" slot on the list take F1B or F2B
	if (i == "A") {
		F0_temp <- F0[[i]]
		F1_2 <- c("F1A", "F2A")
	} else if (i == "B") {
		F0_temp <- F0[[i]]
		F1_2 <- c("F1B", "F2B")
	} else if (i == "combined") {
		F0_temp <- c(paste(c(F0[["A"]][[1]], F0[["B"]][[1]]), collapse = "$|"), paste(c(F0[["A"]][[2]], F0[["B"]][[2]]), collapse = "$|"))
		F1_2 <- c("F1", "F2")
	}	
	# Extract the average (col 8 of the table) nb of cercariae produced for each parents of the crosses
	for (j in F0_temp){
		a <- grep(paste0(j, "$"), mydataF0[,1])  
		av_F0[ match(j, F0_temp) ] <- mean(mydataF0 [a, 8])
	}
	
	# Compute the variance on the average of cercariae produced for each generation (F1 and F2) for each crosses
	for (j in F1_2){ 
		b <- grep(j, mydataF1F2[,2]) 
		var_F[ match(j, F1_2) ] <- var(mydataF1F2 [b, 8])
	} 
	
	# Compute the Castle-Wright estimator (Ne)
	Ne <- ((av_F0[1]-av_F0[2])^2)/((var_F[2]-var_F[1])*8)
    cat("Castle-Wright estimator for cross ", i, ": ", Ne, "\n", sep="")
}

