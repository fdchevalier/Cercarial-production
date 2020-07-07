#!/usr/bin/env Rscript
# Title: PZQ_analysis_func.R
# Version: 0.0
# Author: Frédéric CHEVALIER <fcheval@txbiomed.org>
# Created in: 2019-05-15
# Modified in: 



#==========#
# Comments #
#==========#

# v0.0 - 2019-05-15: store functions in a separate file from the inital script



#======================#
# Packages and options #
#======================#

library("rtracklayer")
library("magrittr")
library("vcfR")

options(stringsAsFactors = FALSE)



#===========#
# Functions #
#===========#

#------------#
# Save plots #
#------------#

save.myplot <- function (filename, width=8, height=8) {
	# Folder creation for graphic output (if does not exist)
	if (file.exists("Graphs/") == FALSE) {dir.create("Graphs")}

	# Graphic saving (png and eps formats)
	dev.copy(png,paste("Graphs/",filename,".png",sep=""), width=width*72, height=height*72, bg="white")
	dev.off()
	dev.copy2pdf(file=paste("Graphs/",filename,".pdf",sep=""), width=width, height=height)
}


#-------------#
# Save tables #
#-------------#

save.mytable <- function (x, filename, ext=".tsv", sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE){
	# Folder creation for table output (if does not exist)
	if (file.exists("Results/") == FALSE) {dir.create("Results")}
	
	# Data export
	write.table(x, file=paste("Results/",filename, ext, sep=""), col.names=col.names, row.names=row.names, append=append, sep=sep)
}


#----------#
# Filename #
#----------#

# Improvement of file_path_sans_ext from the tools package
filename <- function (x, compression=FALSE, length=TRUE) {
    # Test if compression extension
    if (compression) {x <- sub("[.](gz|bz2|xz)$", "", x)}
    # Split name
    myname <- strsplit(x,"\\.")[[1]]

    # Set length of the filename to keep if several point
    if (isTRUE(length)) {
        length <- length(myname)-1
    } else if (! is.numeric(length)) {
        stop("Length must be numeric.")
    } else if (length > length(myname)) {
        stop("Specified lenght is greater than the actual length.")
    }

    basename(paste(myname[0:length], collapse='.'))
}



#===========================#
# Data processing functions #
#===========================#

#--------------#
# Peaks report #
#--------------#

# Peaks report
pk.rpt <- function ( x, sign, trsh, cln.print=NULL, cln.nm.ext=NULL, res.folder, filename=paste("peaks_report.",x.nm,".",cln.nm.ext,sep="") ) {
    # Usage
    ## x: dataframe to analyse
    ## sign: >, <, >=, <= or ==
    ## trsh: trehsold value
    ## cln.print: vector containing the column number or name to print in the report (default all)
    ## cln.nm.ext: name to be printed in the filename of the report (default, name from cln.print)
    ## res.folder: folder to store tables
    ## filename: filename of the table
    
    # Check steps
    ## Check mandatory objects
    if (missing(x) || missing(sign) || missing(trsh) || missing(res.folder)) {stop("x, sign, trsh and res.folder are mandatory.", call.=FALSE)}
    ## Check class of objects
#    if (! is.data.frame(x)) {stop("x must be a data.frame.", call.=FALSE)}
    if (! is.numeric(trsh) || length(trsh) > 1) {stop("trsh must be a single numerical value.", call.=FALSE)}
    if (! is.null(cln.print) && ! is.vector(cln.print)) {stop("cln.print must be a vector.", call.=FALSE)}
#    if (! is.null(cln.nm.ext) && ! (is.vector(cln.nm.ext) || is.character(cln.nm.ext))) {stop("cln.nm.ext must be a vector.", call.=FALSE)}
    if (! is.null(cln.nm.ext) && ! is.character(cln.nm.ext) && ! length(cln.nm.ext) > 1) {stop("cln.nm.ext must be a single character string.", call.=FALSE)}
    if (! is.character(res.folder) || length(res.folder) > 1) {stop("res.foler must be a character string of length 1.")}

    ## Check signs
    if (! is.character(sign) || ! (sign == ">" || sign == "<" || sign == "<=" || sign == ">=" || sign == "==")) {stop("sign must be >, <, >=, <= or ==.", call.=FALSE)}
    
    ## Check for directory
    if (file.exists(res.folder) == FALSE) {dir.create(res.folder, recursive=TRUE)}

    # x name 
    x.splt <- unlist(strsplit(deparse(substitute(x)), "\\[|,|\\]$"))[1]
    x.nm <- x.splt[1]

    # Original data.frame
    x.data <- get(x.nm)

    # Column range and names
    if (is.null(cln.print)) { cln.print <- colnames(get(x.nm)) }

    # Table generation
    if (! is.data.frame(x)) { x <- as.data.frame(x) }
    myrows <- apply(x, 1, function(z) any(do.call (sign, list(z, trsh))))
    report <- x.data[ myrows , cln.print ]

    # Table export
    write.table(report, paste0(res.folder,"/",filename), row.names=FALSE, quote=FALSE, sep="\t")
}


#-------------------------------------------------------------#
# Detail and summary tables listing genes and variants in QTL #
#-------------------------------------------------------------#

# qtl.tb <- function(x, x.fltr, chr, pv.cln, bf.cor, ann.tb, expr.tb, expr.cln, expr.nm=NULL, cln.print) {
qtl.tb <- function(x, chr, ann.tb, expr.tb, expr.cln, expr.nm=NULL, cln.print) {

    # Usage
    ## x            a table with allele frequencies, genotypes, and p-values
    ## x.flt        filtered x table
    ## chr          chromosome to filter on
    ## pv.cln       column(s) that contains p-values to filter on       
    ## bf.cor       Bonferroni correction
    ## ann.tb       annotation table
    ## expr.tb      gene expression table
    ## expr.cln     column of the gene expression table to include
    ## expr.nm      name of expression stage
    ## cln.print    column to print in the final table

    cat("QTL tables: processing ", chr, "...\n", sep="")

    # Select rows in the QTL region
    x <- x[ x[,1] == chr, ]

    if(is.null(expr.nm)) { expr.nm <- colnames(expr)[expr.cln] }

    k <- length(expr.cln)

    myinfo <- strsplit(x[,8], ";", fixed = TRUE)
    
    # Prepare parallelization
    myinfo.lg <- 1:length(myinfo)
    n <- detectCores()-1
    myinfo.lg <- split(myinfo.lg, sort(myinfo.lg%%n))
    
    old.workers <- getDoParWorkers()
    registerDoParallel()
    

    myinfo.ann <- foreach(i=1:length(myinfo.lg), .combine='c') %dopar% {
        lapply(myinfo[myinfo.lg[[i]]], function(x) unlist(grep("ANN",x, value=T) %>% strsplit(., "|", fixed=T)))
    }

    mygenes.ls <- vector("list", length(myinfo.ann))


    mygenes.ls <- foreach(j=1:length(myinfo.lg), .combine='c') %dopar% {
        myidx <- myinfo.lg[[j]]

        mygenes.ls.tmp <- mygenes.ls[myidx]
        
        for (i in 1:length(myidx)) {
            ann <- myinfo.ann[[myidx[i]]]
            
            mycln.tmp <- seq(0, by=15, length.out=round(length(ann)/15))
            
            if (length(mycln.tmp) > 0) {
                mymtx <- sapply(mycln.tmp, function(y) ann[y+c(2,3,4,7,8,10,11)])

                myRNA <- mymtx[3,] %>% gsub("transcript:", "", .) 
                myRNA.tb <- matrix(NA, nrow=2+length(expr.cln), ncol=ncol(mymtx))
                
                mymain.gene <- grep("-", myRNA, invert=TRUE, value=TRUE) %>% strsplit(., "[.]") %>% unlist() %>% .[1]

                for (r in 1:length(myRNA)) {
                    
                    # Corresponding gene
                    mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
                    if (is.na(mygene)) next 
                    
                    # Transcript data
                    ann.tmp <- ann.tb[ ann.tb[,1] == myRNA[r] , ]
                    # If no transcript detected
                    if (length(ann.tmp) == 0) { ann.tmp <- ann[ ann[,1] == mygene , ] }
                    # If really nothing detected
                    if (dim(ann.tmp)[1] == 0) { ann.tmp <- matrix(NA, ncol=3, nrow=1) }
                    # Add data to the matrix
                    myRNA.tb[1:2, r] <- t(as.matrix(ann.tmp[,2:3]))

                    # Expression data
                    expr.tmp <- expr.tb[ expr.tb[,1] == myRNA[r] , expr.cln ]
                    # If no expression found
                    if (k == 1 && length(expr.tmp) == 0 ) { expr.tmp <- rep(NA, length(expr.cln)) }
                    if (k > 1  && (length(expr.tmp) == 0 | dim(expr.tmp)[1] == 0)) { expr.tmp <- rep(NA, length(expr.cln)) }
                    # Add data to matrix
                    myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(expr.tmp)

#                    # Add info for main gene table
#                    if (mygene == mymain.gene & isTRUE(todo)) {
#                        mymain.genes.tb[i,] <- cbind(mygene, t(myRNA.tb[,r]))
#
#                        # If expression is NA, is there any other expression data for this given gene?
#                        if (all(is.na(expr.tmp)) | sum(expr.tmp) == 0) {
#                            expr.tmp <- expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ]
#                            if (k == 1 && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- max(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE) }
#                            if (k > 1  && ! all(is.na(expr.tmp))) { mymain.genes.tb[i,4:ncol(mymain.genes.tb)] <- apply(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
#                        }
#
#                        # Don't do anything else for this gene in the future
#                        todo <- FALSE
#                    }
                }

                # Combine matrix
                mymtx <- rbind(mymtx, myRNA.tb)
                
                # Reorder rows
                mymtx <- mymtx[c(3, 8:(9+k), 4, 1:2, 5:7), ]
                
#                mygenes.ls[[i]] <- mymtx
                mygenes.ls.tmp[[i]] <- mymtx
            }
        }

        return(mygenes.ls.tmp)

#        rm(mygenes.ls.tmp)
#        gc(reset=TRUE)

    }
    
    cln.nm <- c("Gene", "GFF_annotation", "HHPred_annotation", expr.nm)
    mymain.genes.tb <- foreach(j=1:length(myinfo.lg), .combine='rbind') %dopar% {
    # for (j in (1:length(myinfo.lg))) {
        myidx <- myinfo.lg[[j]]

        mymain.genes.tb.tmp <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(myidx)))
        
        for (i in 1:length(myidx)) {
            ann <- myinfo.ann[[myidx[i]]]

            mycln.tmp <- seq(0, by=15, length.out=round(length(ann)/15))
            
            if (length(mycln.tmp) > 0) {
                mymtx <- sapply(mycln.tmp, function(y) ann[y+c(2,3,4,7,8,10,11)])

                myRNA <- mymtx[3,] %>% gsub("transcript:", "", .) # Change in the GFF file
                myRNA.tb <- matrix(NA, nrow=2+length(expr.cln), ncol=ncol(mymtx))
                
                mymain.gene <- grep("-", myRNA, invert=TRUE, value=TRUE) %>% strsplit(., "[.]") %>% unlist() %>% .[1]

                for (r in 1:length(myRNA)) {
                    
                    # Corresponding gene
                    mygene <- unlist(strsplit(myRNA[r], ".", fixed=TRUE))[1]
                    if (is.na(mygene)) next 
                    
                    # Transcript data
                    ann.tmp <- ann.tb[ ann.tb[,1] == myRNA[r] , ]
                    # If no transcript detected
                    if (length(ann.tmp) == 0) { ann.tmp <- ann[ ann[,1] == mygene , ] }
                    # If really nothing detected
                    if (dim(ann.tmp)[1] == 0) { ann.tmp <- matrix(NA, ncol=3, nrow=1) }
                    # Add data to the matrix
                    myRNA.tb[1:2, r] <- t(as.matrix(ann.tmp[,2:3]))

                    # Expression data
                    expr.tmp <- expr.tb[ expr.tb[,1] == myRNA[r] , expr.cln ]
                    # If no expression found
                    if (k == 1 && length(expr.tmp) == 0 ) { expr.tmp <- rep(NA, length(expr.cln)) }
                    if (k > 1  && (length(expr.tmp) == 0 | dim(expr.tmp)[1] == 0)) { expr.tmp <- rep(NA, length(expr.cln)) }
                    # Add data to matrix
                    myRNA.tb[3:nrow(myRNA.tb), r] <- as.matrix(expr.tmp)

                    # Add info for main gene table
                    # mymain.genes.tb.tmp[i,] <- cbind(mygene, t(myRNA.tb[,r]))
                    if (length(mymain.gene) > 0 && mygene == mymain.gene) {
                    # if (! is.null(mymain.gene) & mygene == mymain.gene) {
                        mymain.genes.tb.tmp[i,] <- cbind(mygene, t(myRNA.tb[,r]))

                        # If expression is NA, is there any other expression data for this given gene?
                        if (all(is.na(expr.tmp)) | sum(expr.tmp) == 0) {
                            expr.tmp <- expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ]
                            if (k == 1 && ! all(is.na(expr.tmp))) { mymain.genes.tb.tmp[i,4:ncol(mymain.genes.tb.tmp)] <- max(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE) }
                            if (k > 1  && ! all(is.na(expr.tmp))) { mymain.genes.tb.tmp[i,4:ncol(mymain.genes.tb.tmp)] <- apply(expr.tb[ grep(paste0(mygene,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE)) }
                        }

                        break
                    }
                }

            }
        }

        return(mymain.genes.tb.tmp)

    }
    
    colnames(mymain.genes.tb) <- cln.nm

    # Reset initial workers
    registerDoParallel(old.workers)

    mygenes.ls[ sapply(mygenes.ls, is.null) ] <- NA
    mymax      <- max(sapply(mygenes.ls, length))
    mygenes.tb <- lapply(mygenes.ls, function(x) x[1:mymax])
    mygenes.tb <- t(do.call(cbind, mygenes.tb))
    mygenes.tb[ mygenes.tb == "" ] <- NA
    mygenes.tb <- as.data.frame(mygenes.tb)

    cln.nm.b <- c("Gene (hit ", "GFF annotation (hit ", "HHPred annotation (hit ", paste(expr.nm, "(hit "), "Transcript (hit ", "Region (hit ", "Impact (hit ", "Type (hit ", "DNA mutation (hit ", "Protein mutation (hit ")
    cln.nm <- NULL
    for (i in 1:(mymax/length(cln.nm.b))) { cln.nm <- c(cln.nm, paste0(cln.nm.b, i, ")")) }
    colnames(mygenes.tb) <- cln.nm

    # Final table
    myfinal.tb <- cbind(x[, c(1,2,4,5)], mymain.genes.tb, x[, cln.print], mygenes.tb)
    myfinal.tb.qtl <- myfinal.tb 


    ## Select only informative columns then
    myclns <- ! apply(myfinal.tb.qtl, 2, function(x) all(is.na(x)))
    myfinal.tb.qtl <- myfinal.tb.qtl[ , myclns ]

    return(invisible(myfinal.tb.qtl))

}


#---------------#
# Summary table #
#---------------#

qtl.tb.sum <- function(x, ann.tb, expr.tb, expr.cln, expr.nm=NULL, gff, baits.bed=NULL) {

    # Usage
    ## x            QTL table from qtl.tb
    ## gff          GFF file
    ## baits.bed    BED coordinated of the baits
    
    chr <- unique(x[,1])
    bf.lim <- c(as.numeric(x[1,2]), as.numeric(x[nrow(x),2]))

    cat("Summary QTL table: processing ", chr, "...\n", sep="")
    
    if(is.null(expr.nm)) { expr.nm <- colnames(expr)[expr.cln] }

    k <- length(expr.cln)

    # GFF
    gff.genes <- mygff[ mygff[,3] == "gene", ]
    gff.genes.qtl <- gff.genes[gff.genes[,1] == chr & gff.genes[,5] >= bf.lim[1] & gff.genes[,4] <= bf.lim[2],]

    # mygenes.occ <- unique(gff.genes.qtl[,9])
    mygenes.occ <- unique(gff.genes.qtl[,10])
    mygenes.occ <- na.omit(unique(c(mygenes.occ, x[,5])))

    myimpct.cln <- grep("Impact", colnames(x))
    myimpct.tb  <- as.data.frame(matrix(c("LOW", "MODERATE", "HIGH", "MODIFIER", 0.3, 0.6, 0.9, 0.3), ncol=2))

    # myfreq.cln  <- grep("pol.freq",colnames(x))
    mylod.cln <- grep("lod", colnames(x))

    cln.nm   <- c("Chr.", "Gene ID", "v5 gene ID", "Start", "End", "GFF annotation", "HHPred annotation", expr.nm, "Captured", "Nb of variable sites", "Impact score", "Weighted impact score", "Mean LOD", "Median LOD", "Max LOD", "Global score gene", "Global score CDS")
    mysum.tb <- as.data.frame(matrix(NA, ncol=length(cln.nm), nrow=length(mygenes.occ)))
    for (g in mygenes.occ) {
    # for (g in "Smp_057210") {
        myidx <- match(g, mygenes.occ)

        if (! is.null(baits.bed)) {
            # Is the gene in the capture array?
            # g.crdt <- as.data.frame(mygff.genes[ mygff.genes[,9] == g, c(1,4:5) ])
            g.crdt <- as.data.frame(mygff.genes[ mygff.genes[,10] == g, c(1,4:5) ])
            if (dim(g.crdt)[1] > 0) {
                capture <- any( findInterval(baits.bed[ baits.bed[,1] == g.crdt[,1], 3], g.crdt[,2:3]) == 1 )
            } else {
                capture <- FALSE
            }
            if (capture) { capture <- "yes" } else { capture <- "no" }
        } else {
            capture <- NA
        }
        
        # Is there any GFF v5 ID(s) associated to the current gene?
        # old.id <- as.data.frame(gff.genes[ gff.genes[,9] == g, 13 ])[,3]
        old.id <- as.data.frame(gff.genes[ gff.genes[,10] == g, 13 ])[,3]
        if (length(old.id) == 0) { 
            old.id <- NA 
        } else {
            old.id <- paste(old.id, collapse=", ")
        }
        
        # Pull out information if the gene is in the variant data table
        if (any(x[,5] == g, na.rm=TRUE)) {


            #------#
            # Gene #
            #------#

            mytb.tmp <- x[ ! is.na(x[,5]) & x[,5] == g, ]
        
            mypos <- as.data.frame(gff.genes[ grepl(g, gff.genes[,9]), c(4,5) ])
            if (dim(mypos)[1] == 0) { mypos <- c(NA,NA) }

            nb.var <- nrow(mytb.tmp)

            # Impact related stat
            mytb.impct <- mytb.tmp[,myimpct.cln]
            nb.iso     <- apply(mytb.impct, 2, function(x) ! all(is.na(x))) %>% sum()
            for (i in 1:nrow(myimpct.tb)) { mytb.impct[ mytb.impct == myimpct.tb[i,1] ] <- myimpct.tb[i,2] }
            myimpct.vec <- na.omit(as.numeric(unlist((mytb.impct))))
            myimpct.sc  <- sum(myimpct.vec)
            myimpct.sc2 <- round(sum(myimpct.vec) / ((length(myimpct.vec) / nb.var) * nb.iso), digits=1)

            # LOD related stat
            if (all(is.na(mytb.tmp[, mylod.cln]))) {
                mylod.mean <- mylod.median <- mylod.max <- NA
            } else {
                mylod.mean   <- mean(mytb.tmp[, mylod.cln], na.rm = TRUE)
                mylod.median <- median(mytb.tmp[, mylod.cln], na.rm = TRUE)
                mylod.max    <- max(mytb.tmp[, mylod.cln], na.rm = TRUE)
            }

            # Global score
            mytb.sc <- mytb.tmp[, 11+k+(1:2)]
            mytb.sc[mytb.sc == 0] <- 1
            mytb.sc <- mytb.sc[, 1] * mytb.sc[, 2]
            mytb.sc[mytb.sc == 1] <- 0
            myimpct.vec <- apply(mytb.impct, 1, function(x) max(x, na.rm=TRUE)) %>% as.numeric()
            mytb.sc <- mytb.sc * myimpct.vec
            mygb.sc <- sum(mylod.max, sum(mytb.sc) / diff(unlist(mypos)), na.rm = TRUE)
            # mygb.sc <- sum(mylod.max, sum(mytb.tmp[, 11+k+(1:2)]), na.rm = TRUE)
            if (sum(as.numeric(mytb.tmp[1, 8:(7+k)])) == 0) { mygb.sc <- -mygb.sc }


            #-----#
            # CDS #
            #-----#

            # Exon related information
            cds.pos <- gff[grepl(paste0("cds:", g), gff[,"ID"]), ][,4:5] %>% as.matrix() %>% unique() %>% .[order(.[,1]), , drop = FALSE]

            # Check if some CDS are contained in others
            if (length(unique(cds.pos[,1])) != length(cds.pos[,1]) | length(unique(cds.pos[,2])) != length(cds.pos[,2])) {
                cds.db <- matrix(NA, ncol = nrow(cds.pos), nrow = nrow(cds.pos))
                for (i in 1:nrow(cds.pos)) {
                    cds.db[, i] <- apply(cds.pos, 1, function(x) findInterval(x, cds.pos[i,])) %>% colSums(.) == 3
                }

                # Update CDS position to remove nested CDS
                cds.pos <- cds.pos[! colSums(cds.db) > 1, ]
            }

            if (dim(cds.pos)[1] == 0) {
                myimpct.sc  <- NA
                myimpct.sc2 <- NA
                mygb.cds.sc <- NA
            } else {
                cds.vec      <- apply(cds.pos, 1, function(x) findInterval(mytb.tmp[, 2], x) == 1) %>% as.matrix() %>% apply(., 1, any)
                mytb.cds.tmp <- mytb.tmp[cds.vec, ]

                cds.length <- apply(cds.pos, 1, diff) %>% sum()

                cds.nb.var <- nrow(mytb.cds.tmp)

                # Impact related stat
                mytb.impct <- mytb.cds.tmp[,myimpct.cln]
                if (dim(mytb.impct)[1] > 0) {
                    for (i in 1:nrow(myimpct.tb)) { mytb.impct[ mytb.impct == myimpct.tb[i,1] ] <- myimpct.tb[i,2] }
                    myimpct.vec <- na.omit(as.numeric(unlist((mytb.impct))))
                }

                # Global score
                mytb.sc <- mytb.cds.tmp[, 11+k+(1:2)]
                if (dim(mytb.sc)[1] == 0) {
                    mygb.cds.sc <- 0
                } else {
                    mytb.sc[mytb.sc == 0] <- 1
                    mytb.sc <- mytb.sc[, 1] * mytb.sc[, 2]
                    mytb.sc[mytb.sc == 1] <- 0
                    myimpct.vec <- apply(mytb.impct, 1, function(x) max(x, na.rm=TRUE)) %>% as.numeric()
                    mytb.sc <- mytb.sc * myimpct.vec
                    mygb.cds.sc <- sum(mylod.max, sum(mytb.sc) / cds.length, na.rm = TRUE)
                    if (sum(as.numeric(mytb.tmp[1, 8:(7+k)])) == 0) { mygb.cds.sc <- -mygb.cds.sc }
                }
            }


            #----------------#
            # Building table #
            #----------------#

            mysum.tb[myidx,1:2]  <- mytb.tmp[1,c(1,5)] 
            mysum.tb[myidx,3]    <- old.id
            mysum.tb[myidx,4:5]  <- mypos 
            mysum.tb[myidx,6:7]  <- mytb.tmp[1,6:7]
            mysum.tb[myidx,8:(7+k)]  <- as.matrix(mytb.tmp[1,8:(7+k)])
            mysum.tb[myidx,8+k]  <- capture
            mysum.tb[myidx,9+k]  <- nb.var
            mysum.tb[myidx,10+k] <- myimpct.sc
            mysum.tb[myidx,11+k] <- myimpct.sc2
            mysum.tb[myidx,12+k] <- mylod.mean
            mysum.tb[myidx,13+k] <- mylod.median
            mysum.tb[myidx,14+k] <- mylod.max
            mysum.tb[myidx,15+k] <- mygb.sc
            mysum.tb[myidx,16+k] <- mygb.cds.sc

        } else {
            myexpr.tmp <- expr.tb[ grep(g, expr.tb[,1]) , expr.cln ]
            if (k == 1 && ! all(is.na(myexpr.tmp))) {
                myexpr.tmp <- max(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], na.rm=TRUE)
            } else if (k > 1  && ! all(is.na(myexpr.tmp))) {
                myexpr.tmp <- apply(expr.tb[ grep(paste0(g,"."), expr.tb[,1], fixed=TRUE) , expr.cln ], 2, function(x) max(x, na.rm=TRUE))
            } else {
                myexpr.tmp <- rep(NA, k)
            }

            mysum.tb[myidx,1]        <- as.character(gff.genes.qtl[ grep(g, gff.genes.qtl[,9]), 1 ])
            mysum.tb[myidx,c(2,4,5)] <- as.data.frame(gff.genes.qtl[ grep(g, gff.genes.qtl[,10]), c(10,4:5) ])
            mysum.tb[myidx,3]        <- old.id
            mysum.tb[myidx,6:7]      <- ann.tb[ grep(g,ann.tb[,1]), 2:3 ][1,]
            mysum.tb[myidx,8:(7+k)]  <- myexpr.tmp
            mysum.tb[myidx,8+k]      <- capture

        }

    }

    # Remove empty rows
    mysum.tb <- as.data.frame(mysum.tb[ ! rowSums(is.na(mysum.tb)) == ncol(mysum.tb), ])

    # Name column
    colnames(mysum.tb) <- cln.nm 

    myrows <- na.omit(rownames(mysum.tb)[mysum.tb[,1] == chr & mysum.tb[,5] >= bf.lim[1] & mysum.tb[,4] <= bf.lim[2]])
    mysum.tb.qtl <- mysum.tb[ myrows, ]

    return(invisible(mysum.tb.qtl))
}

