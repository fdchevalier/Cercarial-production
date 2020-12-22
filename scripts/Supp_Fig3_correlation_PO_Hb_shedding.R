#-------------------
# Loading packages
#-------------------

library("gplots")
library("plotrix")

#-----------------------------
# Loading dataset
#-----------------------------

F1 <-read.table("F1.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")
F2 <-read.table("F2.csv", header = TRUE, sep = ",", dec = ".", na.strings = "NA")


#===========#
# Functions #
#===========#

# Line in units
## source: https://stackoverflow.com/a/30835971
line2user <- function(line, side) {
    lh <- par('cin')[2] * par('cex') * par('lheight')
    x_off <- diff(grconvertX(c(0, lh), 'inches', 'npc'))
    y_off <- diff(grconvertY(c(0, lh), 'inches', 'npc'))
    switch(side,
        `1` = grconvertY(-line * y_off, 'npc', 'user'),
        `2` = grconvertX(-line * x_off, 'npc', 'user'),
        `3` = grconvertY(1 + line * y_off, 'npc', 'user'),
        `4` = grconvertX(1 + line * x_off, 'npc', 'user'),
        stop("Side must be 1, 2, 3, or 4", call.=FALSE))
}

#-----------
# Variables
#-----------



#---------
# Figures
#---------

pdf(file="Supplementary_Figure3.pdf", width=6, height=8, useDingbats=FALSE)

layout(matrix(c(1,2,3,4,5,6,7,4),4,2), heights = c(0.5,0.5,0.5,0.1))

#--------------------------------------------------------------------#
#------F1 - Correlation Average cercariae/Laccase activity-----------#

par(mar=c(5,6,3,1))

plot(F1$PO, F1$average, axes=F, xlab="Total laccase activity", ylab="Average cercariae", main="F1", cex.main=2, cex.lab=1.3, ylim=c(0,4000), xlim=c(1,3), frame.plot=FALSE, bg=ifelse(F1$cross=="F1A", "gray50", "white"), pch=ifelse(F1$cross=="F1B", 21, 21), xpd=TRUE)

#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(average~PO, data=F1)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(2,4000, expression(italic("Pearson's r= -0.090; p = 0.42")), cex=1, xpd=TRUE) 

#---Add letters on the figure panel---#
mtext("A.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------------------------------------------------#
#-------F1- Correlation Average cercariae/Hb rate------------#

par(mar=c(5,6,3,1))

plot(F1$Hb, F1$average, xlab="Hemoglobin rate", ylab="Average cercariae",axes=F, cex.lab=1.3, ylim=c(0,4000), xlim=c(0,3),frame.plot=FALSE, bg=ifelse(F1$cross=="F1A", "gray50", "white"), pch=ifelse(F1$cross=="F1B", 21, 21), xpd=TRUE, )


#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(average~Hb, data=F1)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(1.5,4000, expression(italic("Pearson's r= 0.050; p = 0.66")), cex=1, xpd=TRUE) 

#---Add letters on the figure panel---#
mtext("B.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------------------------------------------------#
#-------F1 - Correlation Laccase activity/Hb rate------------#

par(mar=c(5,6,3,1))

plot(F1$PO, F1$Hb, xlab="Total laccase activity", ylab="Hemoglobin rate", xlim=c(1,3), ylim=c(0,3), cex.lab=1.3, axes=F, frame.plot=FALSE, bg=ifelse(F1$cross=="F1A", "gray50", "white"), pch=ifelse(F1$cross=="F1B", 21, 21), xpd=TRUE, )


#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(Hb~PO, data=F1)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(2,3, expression(italic("Pearson's r= 0.72; p < 0.0001")), cex=1, xpd=TRUE) 

#---Add letters on the figure panel---#
mtext("C.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#-----------------------------#
# Legende of the figure panel
#-----------------------------#

par(mar=c(0,0,0,0)) #retravaille les marges du layout
plot(1, 1, type = "n", axes = FALSE, ann = FALSE) #graphe vide

legend("center",c("F1A", "F1B", "F2A", "F2B"),pch=c(21,21,24,24), col="black", pt.bg=c("gray50","white", "gray50","white"), hor=TRUE, bty="n", cex=1.5)

#--------------------------------------------------------------------#
#------F2 - Correlation Average cercariae/Laccase activity-----------#

par(mar=c(5,3,3,1))

ylim <- c(0,max(F2$average))
xlim <- c(1,max(F2$PO, na.rm=T))

plot(F2$PO, F2$average, axes=F, xlab="Total laccase activity", ylab="", main="F2", ylim=ylim, xlim=xlim, cex.main=2, cex.lab=1.3,  frame.plot=FALSE, bg=ifelse(F2$cross=="F2A", "gray50", "white"), pch=ifelse(F2$cross=="F2B", 24, 24), xpd=TRUE)

#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(average~PO, data=F2)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(2.5,4200, expression(italic("Pearson's r= -0.56; p < 0.0001")), cex=1, xpd=TRUE) 

#---Add letters on the figure panel---#
mtext("D.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------------------------------------------------#
#-------F2- Correlation Average cercariae/Hb rate------------#

par(mar=c(5,3,3,1))

ylim <- c(0,max(F2$average))
xlim <- c(0,max(F2$Hb, na.rm=T))

plot(F2$Hb, F2$average, xlab="Hemoglobin rate", ylab="",axes=F, cex.lab=1.3, ylim=ylim, xlim=xlim, frame.plot=FALSE, bg=ifelse(F2$cross=="F2A", "gray50", "white"), pch=ifelse(F2$cross=="F2B", 24, 24), xpd=TRUE)


#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(average~Hb, data=F2)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(1.5,4200, expression(italic("Pearson's r= -0.50; p < 0.0001")), cex=1, xpd=TRUE)

#---Add letters on the figure panel---#
mtext("E.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#------------------------------------------------------------#
#-------F2 - Correlation Laccase activity/Hb rate------------#

par(mar=c(5,3,3,1))

ylim <- c(0,max(F2$Hb, na.rm=T))
xlim <- c(1,max(F2$PO, na.rm=T))

plot(F2$PO, F2$Hb,  xlab="Total laccase activity", ylab="", xlim=xlim, ylim=ylim, cex.lab=1.3, axes=F, frame.plot=FALSE, bg=ifelse(F2$cross=="F2A", "gray50", "white"), pch=ifelse(F2$cross=="F2B", 24, 24), xpd=TRUE)


#---Axis definition---#
axis(1, cex.axis=1, las=1)
axis(2, las=1)

#---Linear regretion---#
myreg <- lm(Hb~PO, data=F2)
abline(myreg, col="red", lty=2)

#---Add statistical analysis on the plot---#
text(2,2.7, expression(italic("Pearson's r= 0.66; p < 0.0001")), cex=1, xpd=TRUE)

#---Add letters on the figure panel---#
mtext("F.", side=3, line=0.5, at=line2user(3,2), cex=par("cex")*2, adj=1)

#-----------------------------------
# Plot saving
#-----------------------------------

dev.off()
