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




mypheno.mt <- matrix(c(
    #pheno.cln, graph_prefix, scan_model
        3,      "d1",         "normal",
        4,      "d2",         "normal",
        5,      "d3",         "normal",
        6,      "d4",         "normal",
        7,      "sum",        "normal",
        8,      "average",    "normal",
        9,      "PO",         "normal",
        10,     "Hb",         "normal", 
        # 11,     "peak_wk",    "normal",
        11,     "traj",    "normal"), ncol=3, byrow=TRUE)       # The pheno.cln cannot start at one because it is the ID column


# Add column to the phenotype table
pheno.tb <- pull.pheno(mydata.qtl)

## Add peak week
pheno.tb[, ncol(pheno.tb) + 1] <- apply(pheno.tb[, 3:6], 1, which.max)

## Add angular momentum
# myang <- (apply(pheno.tb[, 3:6], 1, diff) / 1000) %>% atan(.)
myang <- apply(pheno.tb[, 3:6], 1, diff) 
# pheno.tb[, ncol(pheno.tb) + 1] <- apply(myang, 2, function(x) sum(x, na.rm=T) * diff(range(x, na.rm=T)))
pheno.tb[, ncol(pheno.tb) + 1] <- apply(myang, 2, function(x) sum(x, na.rm=T))

colnames(pheno.tb)[(ncol(pheno.tb) - 1):ncol(pheno.tb)] <- c("peak_wk", "ang_mom")

# # Coefficient of variation
# pheno.tb[, ncol(pheno.tb) ] <- apply(pheno.tb[, 3:6], 1, function(x) sd(x, na.rm=T)/mean(x, na.rm=TRUE)*100)
# colnames(pheno.tb)[(ncol(pheno.tb) - 1):ncol(pheno.tb)] <- c("peak_wk", "cv")

# ## Put back phenotype table
# mydata.qtl$pheno <- pheno.tb



# Tajectories and Frechet distance
library("trajectories")
b <- pheno.tb[ rowSums(is.na(pheno.tb[,3:6])) == 0, ]

lt <- vector("list", nrow(b))
names(lt) <- b[,1]
for (i in 1:nrow(b)) { lt[[i]] <- as.Track(1:4, unlist(b[i,3:6]), seq(as.POSIXct("2015-1-1 0:00"), as.POSIXct("2015-1-1 3:00"), by = "hour")) }

e <- b[,1] %>% as.character() %>% combn(., 2)
tmp <- sapply(split(e, col(e)), function(x) frechetDist(lt[[x[1]]], lt[[x[2]]]))

mstat <- matrix(NA, ncol = nrow(b), nrow = nrow(b))
mstat[lower.tri(mstat)] <- tmp
## directly fill lower triangle of matrix
## need to transpose
tmstat <- t(mstat)
## then fill in lower triangle again, to get correct order
tmstat[lower.tri(tmstat)] <- tmp
## transpose back
mstat <- t(tmstat)
## add on identifiers
colnames(mstat) <- rownames(mstat) <- b[,1]

mymeans <- colMeans(mstat, na.rm=T)

# pheno.tb[, ncol(pheno.tb) +1 ] <- mymeans
# colnames(pheno.tb)[ncol(pheno.tb)] <- c("traj")

b[, ncol(b) +1 ] <- mymeans
colnames(b)[ncol(b)] <- c("traj")

## Put back phenotype table
mydata.qtl$pheno <- merge(pheno.tb, b[,c(1,11)], by=1, all=T) 








# mycomp.ls <- c(myF2.list,"combination")
mycomp.ls <- c("combination")


# List to store results
# mypheno.qtl.ls  <- vector("list", nrow(mypheno.mt))
mypheno.qtl.ls  <- vector("list", nrow(mypheno.mt))
names(mypheno.qtl.ls) <- mypheno.mt[,2]

# Analyze all phenotypes
# for (p in 1:nrow(mypheno.mt)) {
for (p in nrow(mypheno.mt)) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])

    # Methods used for the QTL analysis
    mymethods    <- c("em")
    mymethods.nm <- c("EM")
    if (mymodel != "np") {
        mymethods    <- c(mymethods, "mr")
        mymethods.nm <- c(mymethods.nm, "marker regression")
    }

    # QTL list
    myqtl.ls  <- vector("list", length(mycomp.ls))
    names(myqtl.ls) <- mycomp.ls

    for (i in mycomp.ls) {
        cat(paste(" --Cross",i,"\n"))

        # Subset cross of interest
        if (i != "combination") {
            in.qtl <- subset(mydata.qtl, ind=(pull.pheno(mydata.qtl, cross.cln) == i))
            in.qtl <- drop.markers(in.qtl, myuninfo.mrkr[[i]])
        } else {
            mylist.mkr <- myuninfo.mrkr %>% unlist() %>% unique()
            in.qtl     <- drop.markers(mydata.qtl, mylist.mkr)
        }

        # Genotype probability calculation (needed for some analysis)
        cat("\t -Computing genotype probability...\n")
        in.qtl <- calc.genoprob(in.qtl)
        

        #~~~~~~~~~~~~~~#
        # Scan for QTL #
        #~~~~~~~~~~~~~~#

        out.ls <- vector("list", length(mymethods))
        names(out.ls) <- mymethods

        for (m in mymethods) {
            cat("\t -Performing ", mymethods.nm[mymethods %in% m], " analysis and permutation computation (", my.n.perm, ")...\n", sep="")
            out <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=mymodel)
            cat("\t")
            set.seed(myseed + p + match(i, mycomp.ls) + match(m, mymethods))
            # out.perm <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=mymodel, n.perm=my.n.perm, n.cluster = n.cluster, verbose=FALSE)
            out.perm <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=mymodel, n.perm=100, n.cluster = n.cluster, verbose=FALSE)
            # out.trsh <- as.numeric(sort(out.perm)[round(my.n.perm-my.n.perm*mylod.trsh)])
            out.trsh <- as.numeric(sort(out.perm)[round(100-100*mylod.trsh)])
            if(max(out[,3]) < out.trsh) {warning(m, " method: all LOD scores are under LOD threshold.", immediate.=TRUE, call.=FALSE)}

            # Rename positions
            out[,2] <- unlist(lapply(rownames(out), function(x) rev(strsplit(x,"_")[[1]])[1]))

            # Check for spurious results
            out <- out[ ! is.na(out[,2]),]
            if (any(is.infinite(out[,3]))) {
                warning("Infinite values present. They will be replaced by the maximum value.", immediate.=TRUE, call.=FALSE)
                out[is.infinite(out[,3]),3] <- max(out[!is.infinite(out[,3]),3])
            }

            # Store results
            out.ls[[m]] <- list(lod=out, perm=out.perm, trsh=out.trsh)
        }

        # Creating QTL object
        myqtl.ls[[i]]          <- out.ls
        myqtl.ls[[i]]$genoprob <- in.qtl
    }

    # Store in the designated slot
    mypheno.qtl.ls[[mypheno.mt[p,2]]] <- myqtl.ls

}



# Analyze all phenotypes
# for (p in 1:nrow(mypheno.mt[1:4,])) {
for (p in nrow(mypheno.mt)) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])
    # g.prefix <- mypheno.mt[p, 2 ]
    # g.prefix <- paste(g.prefix, colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")
    g.prefix <- paste(colnames(mydata.qtl$pheno)[pheno.cln], mymodel, paste0("md-", mymissing.data), sep=".")

    myqtl.ls <- mypheno.qtl.ls[[p]]

    for (i in mycomp.ls) {

        for (m in mymethods) {
        
            myqtl.data <- myqtl.ls[[i]][[m]]
            
            # Determine the max of y axis
            if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
                my.ylim <- ceiling(myqtl.data$trsh)
            } else {
                my.ylim <- ceiling(max(myqtl.data$lod[,3]))
            }

            mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
            mylod <- mylod[! grepl("loc", mylod[,2]), ]
            mylod[,2] <- as.numeric(mylod[,2])

            pdf(paste0(graph_fd,g.prefix,"_",i,"_",m,".pdf"), width=15)
                matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=my.ylim, abline.h=myqtl.data$trsh, data.order=TRUE, by.pos=TRUE)
            dev.off()
        }
    }
}


# Superimposed LOD scores for cercarial shedding
p.ls  <- c(3:6,8)-2
# myclr <- viridisLite::viridis(length(p.ls), end=0.8, direction=-1)
myclr <- c(rgb(1, 0, 0, alpha=0.5), gray.colors(length(p.ls) - 1, start=0.05, end=0.75, gamm=1)) %>% rev()
# myclr <- viridisLite::magma(length(p.ls))
# mylwd <- seq(1, 6, length.out=length(p.ls))
mylwd <- c(rep(1, length(p.ls) - 1), 3)
# mylty <- 5:1
for (i in mycomp.ls) {
    for (m in mymethods) {
		
		g.prefix <- paste(i, m, paste0("md-", mymissing.data), sep=".")

		mymax <- lapply(mypheno.qtl.ls[p.ls], function(x) unlist(x[[i]][[m]][["lod"]][,3]) %>% max()) %>% unlist() %>% max()

        pdf(paste0(graph_fd,g.prefix,".pdf"), width=15)
	
        for (p in rev(p.ls)) {
    
			# Phenotype column
			pheno.cln <- as.numeric(mypheno.mt[p, 1])

			myqtl.ls <- mypheno.qtl.ls[[p]]
			
			myqtl.data <- myqtl.ls[[i]][[m]]
            
            # Determine the max of y axis
            if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
                my.ylim <- ceiling(myqtl.data$trsh)
            } else {
                my.ylim <- ceiling(max(myqtl.data$lod[,3]))
            }

            mylod     <- rename_chr_SmV7(myqtl.data$lod, 1) %>% as.data.frame()
            mylod     <- mylod[! grepl("loc", mylod[,2]), ]
            mylod[,2] <- as.numeric(mylod[,2])

            p.i <- match(p, p.ls)
            if (p.i == length(p.ls)) {
                myadd <- FALSE
            } else {
                myadd <- TRUE
            }

            matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=mymax, col=rep(myclr[p.i], 2), data.order=TRUE, by.pos=TRUE, lwd=mylwd[p.i], chr.bg="white", add=myadd)
        }

        legend("topright", legend=c("Week 4 PI", "Week 5 PI", "Week 6 PI", "Week 7 PI", "Average"), col=myclr, bty="n", lty=1, lwd=mylwd, xjust=0) # Need work on the adjustment
        
        dev.off()
    }
}



# Analyze all phenotypes
myrow <- 4
pdf(paste0(graph_fd,g.prefix,"panel_weeks.pdf"), width=10, height=4*myrow)

par(mar=c(5,4.5,3,1.5) + 0.1)

layout(matrix(1:myrow, ncol=1))

for (p in 1:nrow(mypheno.mt[1:4,])) {
    
    # Model for QTL scan
    mymodel <- mypheno.mt[p, 3]

    # Phenotype column
    pheno.cln <- as.numeric(mypheno.mt[p, 1])

    myqtl.ls <- mypheno.qtl.ls[[p]]

    # for (i in mycomp.ls) {
    for (i in "combination") {
		m <- "em"

        my.ylim <- lapply(mypheno.qtl.ls[p.ls], function(x) unlist(x[[i]][[m]][["lod"]][,3]) %>% max()) %>% unlist() %>% max() %>% signif(., 2)

        # for (m in mymethods) {
        for (m in m) {
        
            myqtl.data <- myqtl.ls[[i]][[m]]
            
            # # Determine the max of y axis
            # if (myqtl.data$trsh > max(myqtl.data$lod[,3])) {
            #     my.ylim <- ceiling(myqtl.data$trsh)
            # } else {
            #     my.ylim <- ceiling(max(myqtl.data$lod[,3]))
            # }

            mylod <- rename_chr_SmV7(myqtl.data$lod, 1)
            mylod <- mylod[! grepl("loc", mylod[,2]), ]
            mylod[,2] <- as.numeric(mylod[,2])

            # matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.max=my.ylim, abline.h=myqtl.data$trsh, data.order=TRUE, by.pos=TRUE)
            matplot.data(mylod, 3, datatype="freq", ylab="LOD score", ylim.min=0, ylim.max=my.ylim, data.order=TRUE, by.pos=TRUE)
            mtext(paste0(LETTERS[p], "."), side=3, line=0.75, at=line2user(par("mar")[2],2), cex=par("cex")*2.5, adj=0)
            title(paste("Shedding", p), cex.main=par("cex")*3)
        }
    }
}
dev.off()
