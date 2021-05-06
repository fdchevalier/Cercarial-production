#!/usr/bin/env Rscript
# Title: Fig2_cercarial_shedding_pheno.R
# Version: 0.2
# Author: Winka LE CLEC'H <winkal@txbiomed.org>
# Created in: 2020-03
# Modified in: 2021-05-05


#-------------------
# Loading packages
#-------------------

suppressMessages({
    library("gtools")
    library("gplots")
    library("plotrix")
    library("reshape")
})

#-----------------------------
# Loading data
#-----------------------------

# Suppress warning messages
options(warn=-1)

# Working directory
setwd(file.path(getwd(), "scripts"))

# Functions
source("functions/line2user.R")

# Folders
data_fd   <- "../data/"
graph_fd  <- "../graphs/"

mydataF0 <-read.csv(paste0(data_fd, "phenotyping/F0_parental_populations.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF1 <-read.csv(paste0(data_fd, "phenotyping/F1.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")
mydataF2 <-read.csv(paste0(data_fd, "phenotyping/F2.csv"), header = TRUE, sep = ",", dec = ".", na.strings = "NA")


#------------
# Variables
#------------

mycol <- c("white", "darkorange", "darkolivegreen3")
mynames <- c("Control", "LS", "HS")

c.type <- c("SmLE", "SmBRE")
c.typeF1 <- c("F1A", "F1B")
c.typeF2 <- c("F2A", "F2B")
mytime <- c(3,4,5,6,7,8)
mytime_cerc <- c("Shed.1", "Shed.2", "Shed.3", "Shed.4")
week_mytime_cerc <- c(4,5,6,7)
mycolor <- c("darkolivegreen3","darkorange")
mycolorF1F2 <- c("gray50","white")

mycex.axis <- 1
mycex.text <- 1
my.x.pp.text <- 1.9
 
cx <- 1.5  
x.ttl <- "Week post-infection" 
y.ttl <- "Cercariae produced"  
x.lim <- c(0.5, (length(mytime_cerc)+0.5))
shifts <- c(-0.15,0.15)


#---Matrix of variables---#
			
mymat <- matrix(c("BRE", "darkorange",19, NA,
				"LE","darkolivegreen3",19, NA,
				"F1A", "black",21, "gray50", 
				"F1B", "black",21, "white",
				"F2A", "black",21, "gray50",
				"F2B", "black",21, "white"), ncol=4, byrow=TRUE)

#-----------------------------
# mydata table (Fig A)
#-----------------------------

dataBRE <- mydataF0[mydataF0[,2] == "SmBRE",]
dataLE <- mydataF0[mydataF0[,2] == "SmLE",]

dataF1A <- mydataF1[mydataF1[,2] == "F1A",]
dataF1B <- mydataF1[mydataF1[,2] == "F1B",]

dataF2A <- mydataF2[mydataF2[,2] == "F2A",]
dataF2B <- mydataF2[mydataF2[,2] == "F2B",]

#-----------------------------
# mydata table (Fig B)
#-----------------------------
colnames(mydataF0)[2] <- c("cross")
cerc_dataF0 <- melt(mydataF0[,1:6], id=c("id", "cross"))
cerc_dataF1 <- melt(mydataF1[,1:6], id=c("id", "cross"))
cerc_dataF2 <- melt(mydataF2[,1:6], id=c("id", "cross"))

colnames(cerc_dataF0)[3] <- "time"
colnames(cerc_dataF1)[3] <- "time"
colnames(cerc_dataF2)[3] <- "time"

colnames(cerc_dataF0)[4] <- "nb_cerc"
colnames(cerc_dataF1)[4] <- "nb_cerc"
colnames(cerc_dataF2)[4] <- "nb_cerc"

#-----------------------------
# mydata table (Fig C and D)
#-----------------------------

dataF2A <- dataF2A[order(dataF2A[,8]),]
dataF2B <- dataF2B[order(dataF2B[,8]),]

# Extract parents from crosses for F0 and F1 tables
pF0A <- mydataF0[which(mydataF0[,10] == "A"),]
pF0B <- mydataF0[which(mydataF0[,10] == "B"),]

pF1A <- mydataF1[which(mydataF1[,12] == "A"),]
pF1B <- mydataF1[which(mydataF1[,12] == "B"),]

# Concatenate pF0, pF1 and F2 progeny tables
colnames(pF0A)[2] <- "cross"
colnames(pF0B)[2] <- "cross"

dataA <- rbind(pF0A[,1:8], pF1A[,1:8])
dataB <- rbind(pF0B[,1:8], pF1B[,1:8])

dataF0_F1_F2A <- rbind(dataA,dataF2A[,1:8])
dataF0_F1_F2B <- rbind(dataB,dataF2B[,1:8])

dataA <- NULL

for (i in unique(dataF0_F1_F2A[,2])) {
	# Identify the rows with F0, F1 and F2 
	data.tmp <- dataF0_F1_F2A[dataF0_F1_F2A[,2] == i,]

	# Add NA lines to separate the parasite generations (F0 F1 and F2: except last row of the final table)
	myNA <- NULL
	if (i != "F2A") { 
		myNA <- rep(NA, ncol(data.tmp))
	}
	dataA <- rbind(dataA, data.tmp, myNA)
}


dataB <- NULL

for (i in unique(dataF0_F1_F2B[,2])) {
	# Identify the rows with F0, F1 and F2 
	data.tmp <- dataF0_F1_F2B[dataF0_F1_F2B[,2] == i,]

	# Add NA lines to separate the parasite generations (F0 F1 and F2: except last row of the final table)
	myNA <- NULL
	if (i != "F2B") { 
		myNA <- rep(NA, ncol(data.tmp))
	}
	dataB <- rbind(dataB, data.tmp, myNA)
}

# matrix of segments (first column: start / second column: end)
myseg <- matrix(c(1,2,
				4,5,
				7, nrow(dataA)), ncol=2, byrow=TRUE)
				
#---------
# Figures
#---------

if (! dir.exists(graph_fd)) { dir.create(graph_fd, recursive = TRUE) }

pdf(file=paste0(graph_fd, "Figure2.pdf"), width=15, height=15)

layout(matrix(c(1,2,3,4,5,6,3,4,7,8,3,4),4,3), heights = c(0.5,0.5,0.5,0.5), widths= c(0.35,0.35,0.35))

#--------------------------------------------#
#---Figure A F0: Parental shedding average---#

par(mar=c(4,5,3,1))
boxplot(dataBRE[,8], dataLE[,8], ylab= "Average cercariae produced", main="F0", ylim=c(0,5000), boxwex=0.5,cex.lab=1.5, col=c("darkorange", "darkolivegreen3"),names=c("BRE","LE"),cex.names=1.5, cex.main=2,axes=FALSE) 

axis(1,at=c(1,2),labels=c("SmBRE (LS)","SmLE (HS)"), las=1, tcl=-0.5,cex.axis=1.5)
axis(2, tcl=-0.5)

y.pos <- par("usr")[4]+par("usr")[4]*0.05

# Add statistical analysis on the plot
text(c(1,2), c(max(dataBRE[,8], na.rm=TRUE)+250, max(dataLE[,8],na.rm=TRUE)+250), c("a","d"), xpd=TRUE,cex=1.5)

# Add letters on the panel
mtext("A.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#----------------------------------------------#
#----Figure B F0: Parental shedding per week---#

par(mar=c(4,5,3,1))

y.maxs <- c()
for (t.id in 1:length(mytime_cerc)) {  
    for (i in 1:length(c.type)) {   
       inds <- which(cerc_dataF0$cross == c.type[i] & cerc_dataF0$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(cerc_dataF0$nb_cerc[inds], plot=FALSE)[[1]][5])
    }  
}

y.max <- max(y.maxs)
y.max <- ceiling(y.max/10^floor(log10(y.max)))*10^floor(log10(y.max))
y.lim <- c(min(cerc_dataF0$nb_cerc), y.max)


plot(x=0, y=0, xlim=x.lim, ylim=c(0,6000), main="F0", col='white', xlab=x.ttl, ylab=y.ttl, cex.axis=cx, cex.lab=cx, cex.main=2, axes=FALSE)  
axis(side=1, at=1:length(mytime_cerc), labels=week_mytime_cerc, cex.axis=cx, cex.lab=cx)  # put on the x-axis labels  
axis(2)    #tcl=-0.3

for (t.id in 1:length(mytime_cerc)) {  
  for (i in 1:length(c.type)) {   
   inds <- which(cerc_dataF0$cross == c.type[i] & cerc_dataF0$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
   boxplot(cerc_dataF0$nb_cerc[inds], at=t.id+shifts[i], col=mycolor[i], add=TRUE, bty='n', boxwex=0.5, cex=1, axes=FALSE, outline=FALSE)  
  }  
 }  

# Add statistical analysis on the plot
y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1:4, y.maxs.grp*1.09, "***", xpd=TRUE, cex=1.5)

for (t.id in 1:length(mytime_cerc)) {
	segments(t.id-shifts[2], y.maxs.grp[t.id]*1.025, t.id+shifts[2])
}

# Add letters on the panel
mtext("B.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#--------------------------------------------------------------#				
#---Figure C: Cercarial production per F2A indiv (rank plot)---#	
				
par(mar=c(4,5,3,1))

myse <- apply(dataA[,3:6], 1, function(x) sd(x, na.rm=TRUE) / sqrt(length(na.omit(x))))

plot(1,1, xlim=c(1,nrow(dataA)), ylim=range(dataA[,8], na.rm=TRUE)+range(myse, na.rm=TRUE),type = "n", xpd=TRUE, bty="n", ylab="Cercariae produced", xlab="", axes=F, cex.lab=1.5)

for (i in 1:nrow(myseg)){
	segments(myseg[i,1], par("usr")[3], myseg[i,2], par("usr")[3], xpd=TRUE)
	for (j in myseg[i,]){
		segments(j, y0=par("usr")[3], j, par("usr")[3]+par("cxy")[2]*par("tcl"), xpd=TRUE)
	}
}

axis(2)

# Add F0, F1 and F2 to x axis
text(rowMeans(myseg), par("usr")[3] + par("cxy")[2] / par("tcl"), c("F0", "F1A", "F2A"), xpd=TRUE, cex=1)

# Create an empty vector
mycolvec <- rep(NA, nrow(dataA))
mypchvec <- rep(NA, nrow(dataA))
mybgvec <- rep(NA, nrow(dataA))   
 
# grep values in dataA[,1] and assign the right color following the matrix mycolmat
for (i in 1:nrow(mymat)){
	mycolvec[grep(mymat[i,1], dataA[,1])] <- mymat[i,2]
	mypchvec[grep(mymat[i,1], dataA[,1])] <- as.numeric(mymat[i,3])
	mybgvec[grep(mymat[i,1], dataA[,1])] <- mymat[i,4]
} 
# Replace NA values by "" (no color)
mycolvec[is.na(mycolvec)] <- "white"

for (i in 1:nrow(dataA)){
	myse <- sd(dataA[i,3:6], na.rm=TRUE) / sqrt(length(na.omit(dataA[i,3:6])))
	plotCI(i, dataA[i,8], add = TRUE, uiw = myse, sfrac = 0.003, pch = mypchvec[i], pt.bg = mybgvec[i], col = mycolvec[i], xpd = TRUE)	
}

y.pos <- par("usr")[4]+par("usr")[4]*0.05
mtext("C.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

				
# Add stats on the graph
text(c(3, 5.5), c(3800, 3300), c("*", "NS"), xpd=TRUE,cex=1)


#--------------------------------------------------------------#				
#---Figure D: Cercarial production per F2B indiv (rank plot)---#	
				
par(mar=c(5,5,3,1))

myse <- apply(dataB[,3:6], 1, function(x) sd(x, na.rm=TRUE) / sqrt(length(na.omit(x))))

plot(1,1, xlim=c(1,nrow(dataB)), ylim=range(dataB[,8], na.rm=TRUE)+range(myse, na.rm=TRUE),type = "n", xlab="", ylab="Cercariae produced", xpd=TRUE, bty="n", axes=F, cex.lab=1.5)

for (i in 1:nrow(myseg)){
	segments(myseg[i,1], par("usr")[3], myseg[i,2], par("usr")[3], xpd=TRUE)
	for (j in myseg[i,]){
		segments(j, y0=par("usr")[3], j, par("usr")[3]+par("cxy")[2]*par("tcl"), xpd=TRUE)
	}
}

axis(2)

# Add F0, F1 and F2 to x axis
text(rowMeans(myseg), par("usr")[3] + par("cxy")[2] / par("tcl"), c("F0", "F1B", "F2B"), xpd=TRUE, cex=1)

# Create an empty vector
mycolvec <- rep(NA, nrow(dataB))
mypchvec <- rep(NA, nrow(dataB))
mybgvec <- rep(NA, nrow(dataB))   
 
# grep values in dataA[,1] and assign the right color following the matrix mycolmat
for (i in 1:nrow(mymat)){
	mycolvec[grep(mymat[i,1], dataB[,1])] <- mymat[i,2]
	mypchvec[grep(mymat[i,1], dataB[,1])] <- as.numeric(mymat[i,3])
	mybgvec[grep(mymat[i,1], dataB[,1])] <- mymat[i,4]
} 
# Replace NA values by "" (no color)
mycolvec[is.na(mycolvec)] <- "white"
#mypchvec[is.na(mypchvec)] <- ""
#mybgvec[is.na(mybgvec)] <- "white"

for (i in 1:nrow(dataB)){
	myse <- sd(dataB[i,3:6], na.rm=TRUE) / sqrt(length(na.omit(dataB[i,3:6])))
	plotCI(i, dataB[i,8], add=TRUE, uiw=myse, sfrac= 0.003, pch = mypchvec[i], pt.bg = mybgvec[i], col=mycolvec[i], xpd=TRUE)
	
}

y.pos <- par("usr")[4]+par("usr")[4]*0.05

mtext("D.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

text(c(3, 5.5), c(5000, 3100), c("**", "NS"), xpd=TRUE,cex=1)

#----------------------------------------------#
#---Figure A F1: F1 progeny shedding average---#

par(mar=c(4,1,3,1))
boxplot(dataF1A[,8], dataF1B[,8],ylim=c(0,6000), boxwex=0.5,cex.axis=0.6, col=c("gray50", "white"), names=c("F1A","F1B"),main="F1", cex.main=2, cex.names=1.5,axes=FALSE)


axis(1,at=c(1,2),labels=c("F1A","F1B"), las=1, tcl=-0.5,cex.axis=1.5)
axis(2, tcl=-0.5)

y.pos <- par("usr")[4]+par("usr")[4]*0.05

text(c(1,2), c(max(mydataF1[,8], na.rm=TRUE)+250, max(mydataF1[,8],na.rm=TRUE)+250), c("bc","b"), xpd=TRUE,cex=1.5)

#-----------------------------------------------#
#---Figure B F1: F1 progeny shedding per week---#

par(mar=c(4,1,3,1))

y.maxs <- c()
for (t.id in 1:length(mytime_cerc)) {  
    for (i in 1:length(c.typeF1)) {   
       inds <- which(cerc_dataF1$cross == c.typeF1[i] & cerc_dataF1$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(cerc_dataF1$nb_cerc[inds], plot=FALSE)[[1]][5])
    }  
}
y.max <- max(y.maxs)
y.max <- ceiling(y.max/10^floor(log10(y.max)))*10^floor(log10(y.max))
y.lim <- c(min(cerc_dataF1$nb_cerc, na.rm=TRUE), y.max)

plot(x=0, y=0, xlim=x.lim, ylim=c(0,6000), col='white', xlab=x.ttl, main="F1", cex.main=2, cex.axis=cx, cex.lab=cx, axes=FALSE)  
axis(side=1, at=1:length(mytime_cerc), labels=week_mytime_cerc, cex.axis=cx, cex.lab=cx)  # put on the x-axis labels  
axis(2, tcl=-0.5)   #tcl=-0.3

for (t.id in 1:length(mytime_cerc)) {  
  for (i in 1:length(c.typeF1)) {   
   inds <- which(cerc_dataF1$cross == c.typeF1[i] & cerc_dataF1$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
   boxplot(cerc_dataF1$nb_cerc[inds], at=t.id+shifts[i], col=mycolorF1F2[i], add=TRUE, xaxt='n', yaxt='n', bty='n', boxwex=0.5, cex=1, axes=FALSE, outline=FALSE)  
  }  
 }  
 
# Add statistical analysis on the plot
y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1:4, y.maxs.grp*1.09, c("***", "*", "NS", "*"), xpd=TRUE, cex=1.5)

for (t.id in 1:length(mytime_cerc)) {
	segments(t.id-shifts[2], y.maxs.grp[t.id]*1.025, t.id+shifts[2])
}

#----------------------------------------------#
#---Figure A F2: F2 progeny shedding average---#

par(mar=c(4,1,3,1))
boxplot(dataF2A[,8], dataF2B[,8],ylim=c(0,5000), boxwex=0.5,cex.axis=0.6, col=c("gray50", "white"), names=c("F2A","F2B"),main="F2", cex.main=2, cex.names=1.5,axes=FALSE)
#border= c("darkolivegreen", "darkgoldenrod")

axis(1,at=c(1,2),labels=c("F2A","F2B"), las=1, tcl=-0.5,cex.axis=1.5)
axis(2, tcl=-0.5)

y.pos <- par("usr")[4]+par("usr")[4]*0.05

text(c(1,2), c(max(mydataF2[,8], na.rm=TRUE)+250, max(mydataF2[,8],na.rm=TRUE)+250), c("c","bc"), xpd=TRUE,cex=1.5)



#-----------------------------------------------#
#---Figure B F2: F2 progeny shedding per week---#

par(mar=c(4,1,3,1))

y.maxs <- c()
for (t.id in 1:length(mytime_cerc)) {  
    for (i in 1:length(c.typeF2)) {   
       inds <- which(cerc_dataF2$cross == c.typeF2[i] & cerc_dataF2$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
       y.maxs <- c(y.maxs, boxplot(cerc_dataF2$nb_cerc[inds], plot=FALSE)[[1]][5])
    }  
}
y.max <- max(y.maxs)
y.max <- ceiling(y.max/10^floor(log10(y.max)))*10^floor(log10(y.max))
y.lim <- c(min(cerc_dataF2$nb_cerc, na.rm=TRUE), y.max)

plot(x=0, y=0, xlim=x.lim, ylim=c(0,6000), col='white', xlab=x.ttl, main="F2", cex.main=2, cex.axis=cx, cex.lab=cx, axes=FALSE)  
axis(side=1, at=1:length(mytime_cerc), labels=week_mytime_cerc, cex.axis=cx, cex.lab=cx)  # put on the x-axis labels  
axis(2, tcl=-0.5)  #tcl=-0.3

for (t.id in 1:length(mytime_cerc)) {  
  for (i in 1:length(c.typeF2)) {   
   inds <- which(cerc_dataF2$cross == c.typeF2[i] & cerc_dataF2$time == mytime_cerc[t.id]) # rows for this sub.type and task.type; what we'll plot.  
   boxplot(cerc_dataF2$nb_cerc[inds], at=t.id+shifts[i], col=mycolorF1F2[i], add=TRUE, xaxt='n', yaxt='n', bty='n', boxwex=0.5, cex=1, axes=FALSE, outline=FALSE)  
  }  
 }  

# Add statistical analysis on the plot
y.maxs.grp <- sapply(seq(1, length(y.maxs), length(c.type)), function(x) max(y.maxs[x:(x+length(c.type)-1)]))
text(1:4, y.maxs.grp*1.09, c("***", "**", "**", "NS"), xpd=TRUE, cex=1.5)

for (t.id in 1:length(mytime_cerc)) {
	segments(t.id-shifts[2], y.maxs.grp[t.id]*1.025, t.id+shifts[2])
}
#----------------
# Figure saving
#----------------

dev.off()


