#-------------------
# Loading packages
#-------------------

library("gtools")
library("gplots")
library("plotrix")
#library("Hmisc")


#-----------------------------
# Loading data
#-----------------------------

mydataF0A <-read.csv("CrossA_F0_F1.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF0B <-read.csv("CrossB_F0_F1.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF1 <-read.csv("F1.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF2 <-read.csv("F2.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")

#-----------------------------
# Castle-Wright estimator Ne
#-----------------------------
# Combine the tables F0A/B and F1/F2

mydataF1F2 <- rbind(mydataF1,mydataF2[,1:8])
mydataF0 <- rbind(mydataF0A[1:2,], mydataF0B[1:2,])

# Make a list with the id of each parents of the crosses

F0 <- list("A" = c("LE_15_f","BRE_4_m"), "B" = c("LE_19_m","BRE_2_f"))
mycomp <- c("A", "B", "combined") 

# Ne estimated for Cross A and B (independantly)
#------------------------------------------------

av_F0 <- NULL
var_F <- NULL
Ne <- NULL


for(i in mycomp){
	# if "A" slot of the list take F1A or F2A, if "B" slot on the list take F1B or F2B
	#if(grep("A$", names(F0), value = T) == "A") {
	if (i == "A") {
		F0_temp <- F0[[i]]
		F1_2 <- c("F1A", "F2A")
	} else if (i == "B") {
		F0_temp <- F0[[i]]
		F1_2 <- c("F1B", "F2B")
	} else if (i == "combined") {
		F0_temp <- F0 %>% unlist() %>% strsplit(., "_") %>% lapply(., function(x) x[[1]]) %>% unlist() %>% unique()
		F1_2 <- c("F1", "F2")
	}	
	# Extract the average (col 8 of the table) nb of cercariae produced for each parents of the crosses
	for (j in F0_temp){
		a <- grep(j, mydataF0[,1])  
		av_F0[ match(j, F0_temp) ] <- mean(mydataF0 [ a,8 ] )
	}
	
	# Compute the variance on the average of cercariae produced for each generation (F1 and F2) for each crosses
	for (j in F1_2){ 
		b <- grep(j, mydataF1F2[,2]) 
		var_F[ match(j, F1_2)] <- var(mydataF1F2 [ b, 8 ] )
	} 
	
	# Compute the Castle-Wright estimator (Ne)
	Ne <- ((av_F0[1]-av_F0[2])^2)/((var_F[2]-var_F[1])*8)
	print(Ne)
}



## Average shedding parent LE (#id: LE_15_f)

#LEA <- mydataF0A[mydataF0A[,1]=="LE_15_f",]
#av_LEA <- rowMeans(LEA[,3:6], na.rm=TRUE)

## Average shedding parent BRE (#id: BRE_4_m)
#BREA <- mydataF0A[mydataF0A[,1]=="BRE_4_m",]
#av_BREA <- rowMeans(BREA[,3:6], na.rm=TRUE)

## Variation in shedding F1A
#F1A <- mydataF1[mydataF1[,2] == "F1A",]
#varF1A <- var(F1A[,8])

## Variation in shedding F2A
#F2A <- mydataF2[mydataF2[,2] == "F2A",]
#varF2A <- var(F2A[,8])

## Ne estimation
#Ne_A <- ((av_LEA-av_BREA)^2)/((varF2A-varF1A)*8)


## Average shedding parent LE (#id: LE_19_m)

#LEB <- mydataF0B[mydataF0B[,1]=="LE_19_m",]
#av_LEB <- rowMeans(LEB[,3:6], na.rm=TRUE)

## Average shedding parent BRE (#id: BRE_2_f)
#BREB <- mydataF0B[mydataF0B[,1]=="BRE_2_f",]
#av_BREB <- rowMeans(BREB[,3:6], na.rm=TRUE)

## Variation in shedding F1A
#F1B <- mydataF1[mydataF1[,2] == "F1B",]
#varF1B <- var(F1B[,8])

## Variation in shedding F2A
#F2B <- mydataF2[mydataF2[,2] == "F2B",]
#varF2B <- var(F2B[,8])

## Ne estimation
#Ne_B <- ((av_LEB-av_BREB)^2)/((varF2B-varF1B)*8)


# Ne estimated for combined crosses (A+B)
#----------------------------------------

for(i in ){
	# Grab all the LE and all the BRE for both crosses
	
	
	# Compute the average of nb of cercariae produced for parents LE (LE F0A+F0B) and for BRE (BRE F0A+F0B) parents
	for (j in F0[[i]]){ 
		av_F0[ match(j, F0[[i]]) ] <- mydataF0[mydataF0[,1] == j,8]
	}
	
	# Compute the variance on the average of cercariae produced for each generation (F1A+B and F2A+B) 
	for (j in F1_2){ 
		var_F[ match(j, F1_2)] <- var(mydataF1F2[mydataF1F2[,2] == j, 8])
	}
	
	# Compute the Castle-Wright estimator (Ne)
	Ne <- ((av_F0[1]-av_F0[2])^2)/((var_F[2]-var_F[1])*8)
	print(Ne)
}
