#!/usr/bin/env Rscript
# Title: schisto_sex_graph.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2020-12-01
# Modified in: 



#==========#
# Comments #
#==========#

# Plot read depth of the Z chromosome



#==========#
# Versions #
#==========#

# v0.0 - 2020-12-01: creation



#===========#
# Variables #
#===========#

# Folders
data_fd   <- "../data/libraries/"
graph_fd  <- "../graphs/"
result_fd <- "../results/1-QTL/"

# Samples to plot
myspl <- c("SmBRE4_m", "SmLE19_f", "F2A_117", "F2A_173", "F2B_326", "F2B_344")

myseed <- 1542268



#=================#
# Data processing #
#=================#

cat("Loading and processing data. This may take a while...\n")

# Load data
mydata <- vector("list", length(myspl))
names(mydata) <- myspl
for (i in myspl) {
    mydata[[i]] <- read.delim(paste0(data_fd,  i, "/", i, "_Z.cov"), header = FALSE)
}

# Smoothing data
mydata[1:2]             <- lapply(mydata[1:2], function(x) { x[, 3] <- runmed(x[, 3], 2001) ; return(x) })
mydata[3:length(myspl)] <- lapply(mydata[3:length(myspl)], function(x) { x[, 3] <- runmed(x[, 3], 101) ; return(x) })

# Reducing data by removing low covered regions
mydata <- lapply(mydata, function(x) { x <- x[x[, 3] > 5, ] ; return(x) })

# Sub-sampling
set.seed(myseed)
mydata[1:2] <- lapply(mydata[1:2], function(x) x[sort(sample(1:nrow(x), nrow(x)/2000)), ])

set.seed(myseed)
mydata[3:length(myspl)] <- lapply(mydata[3:length(myspl)], function(x) x[sort(sample(1:nrow(x), nrow(x)/100)), ])

# Max position
x_max <- max(unlist(lapply(mydata, function(x) max(x[,2]))))



#=========#
# Figures #
#=========#

cat("\nDrawing graphs. This may take a while...\n")

png(paste0(graph_fd, "Supp. Fig. 4.png"), width = 5 * 2, height = 3 * length(myspl) / 2, unit = "in", res = 300)
layout(matrix(1:length(myspl), ncol = 2, byrow = TRUE))
for (i in myspl) {
    # Point graph
    plot(mydata[[i]][,3] ~ mydata[[i]][,2], xlab = "Position on chromosome Z (Mb)", ylab = "Read depth", xlim = c(1, x_max), ylim = c(5, 200), main = i, xaxt = "n", pch = 20, col = "grey", log = "y")
    axis(1, at = seq(0, 80, 10) * 1e6, labels = seq(0, 80, 10))

    # Smoothed line
    lw1 <- loess(mydata[[i]][,3] ~ mydata[[i]][,2], span=0.05)
    lines(mydata[[i]][,2], predict(lw1,mydata[[i]][,2]), col = "red", lwd = 3)
}
dev.off()
