
#------------------------#
# Function to save plots #
#------------------------#

save.myplot <- function (filename, width=8, height=8) {
	# Folder creation for graphic output (if does not exist)
	if (file.exists("Graphs/") == FALSE) {dir.create("Graphs")}

	# Graphic saving (png and eps formats)
	dev.copy(png,paste("Graphs/",filename,".png",sep=""), width=width*72, height=height*72, bg="white")
	dev.off()
	dev.copy2pdf(file=paste("Graphs/",filename,".pdf",sep=""), width=width, height=height)
}


#-------------------------#
# Function to save tables #
#-------------------------#

save.mytable <- function (x, filename, ext=".tsv", sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE){
	# Folder creation for table output (if does not exist)
	if (file.exists("Results/") == FALSE) {dir.create("Results")}
	
	# Data export
	write.table(x, file=paste("Results/",filename, ext, sep=""), col.names=col.names, row.names=row.names, append=append, sep=sep)
}


#-------------------#
# Dynamic cross tab #
#-------------------#

mydata.fltr <- function(data.x, data.slct, pattern.data, sign, trsh, missing.data=0.2, na.string=NA) {

    # Usage
    ## data.x: dataframe to filter
    ## data.slct: dataframe used for row selection. Can't be the same as data.
    ## pattern.data: pattern to select column(s) in the dataframe to filter and used for selection
    ## sign: sign for the comparison
    ## trsh: numerical threshold value used for selection
    ## missing.data: proportion of missing data allowed
    ## na.rm: remove any row that contain NA in any selected columns

    # Check steps
    if (missing(data.x) || ! is.data.frame(data.x)) {stop("data.x must be present and must be a dataframe.")}
    if (missing(data.slct) || ! is.data.frame(data.slct)) {stop("data.slct must be present and must be a dataframe.")}
    if (missing(sign) || ! is.character(sign)) {stop("sign must be present and must be a character.")}
    if (missing(trsh) || ! is.numeric(trsh)) {stop("trsh must be present and must be numeric.")}

    # check if same number of row
    # check if same number of column selected
    if (! is.numeric(missing.data) | missing.data < 0 | missing.data >= 1) {stop("missing.data must be numeric and in [0;1[")}    # check if missing data is between 0 and 1

    # Tmp objects
    data.cln <- grep(pattern.data, colnames(data.x))
    data.x.tmp <- data.x[,data.cln]
    data.slct.tmp <- data.slct[,data.cln]

    # Put NA in tmp objects
    data.x.tmp[ ! eval(parse(text=paste("data.slct.tmp",sign,trsh))) ] <- na.string     # Warning: we want the opposite value (!) of the selection to be NA because the expression select the cells that pass the selection criteria

    # Replace values in the original data
    data.x[,data.cln] <- data.x.tmp
    
    # Remove rows with missing data over the threshold
    if (is.na(na.string)) {
        na.mt <- is.na(data.x[,data.cln])
    } else {
        na.mt <- data.x[,data.cln] == na.string
    }
    data.rw.md <- (rowSums(na.mt)/length(data.cln)) <= missing.data
    #data.x <- data.x[data.rw.md,]  # Keep only rows that pass the missing data test
    data.x[! data.rw.md,data.cln] <- na.string
    return(data.x)
}


#-------------------------------#
# Genotyping data quality check #
#-------------------------------#

quality.check <- function(x, est.rf=FALSE) {
    cat("Segregation distortion.\n")
    gt <- geno.table(x)
    gti <- nrow(gt[gt$P.value < 1e-7,])
    cat(paste("Number of sites showing segregation distortion (chi-sq < 1e-7): ", gti, ". These sites should be removed.", sep=""))

    cat("Pairwise genotype comparisons.\n")
    cg <- comparegeno(x)
    cgi <- nrow(which(cg >0.9, arr.ind=T))/2
    cat(paste("Number of individuals showing 0.9 indentity in their genotypes: ", cgi, ". This may indicate a mix-up or a very very closely related individuals.", sep=""))

    # est.rf takes quite a lot of time so it is why it is off by default
    if(est.rf) {
        cat("Pairwise recombinaison fraction (it will take a while...).\n")
        x <- est.rf(x)
        checkAlleles(x)
    }
}


#--------------#
# Scan-one QTL #
#--------------#

myscanone <- function(in.qtl, pheno.cln, methods, model, n.perm, n.cluster=NULL, addcovar=NULL, intcovar=NULL) {

    out.ls <- vector("list", length(methods))
    names(out.ls) <- methods

    for (m in methods) {
        
        # Methods used for the QTL analysis
        if (m == "em") { methods.nm <- "EM" }
        if (m == "mr") { methods.nm <- "marker regression" }

        cat("\t -Performing ", methods.nm, " analysis and permutation computation (", n.perm, ")...\n", sep="")
        out <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=model, addcovar = addcovar, intcovar = intcovar)
        cat("\t")
        set.seed(myseed + p + match(i, mycomp.ls) + match(m, methods))
        out.perm <- scanone(in.qtl, pheno.col=pheno.cln, method=m, model=model, addcovar = addcovar, intcovar = intcovar, n.perm=n.perm, n.cluster = n.cluster, verbose=FALSE)
        out.trsh <- as.numeric(sort(out.perm)[round(n.perm - n.perm * mylod.trsh)])
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

    return(out.ls)
}
