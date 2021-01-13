rename_chr_SmV7 <- function(x, cln) {
    y <- as.character(x[,cln])

    y <- gsub("SM_V7_", "", y)
    y <- gsub("(?!^U)[A-V,X,Y]+", "\\.un.SC_", y, perl = TRUE)
    y <- gsub("^U", "SC_", y)
    y <- gsub("^ZW", "W", y)
    y <- gsub("(?!^SC)^", "Chr_", y, perl = TRUE)

    # x <- cbind(original=x, renamed=y)
    x[,cln] <- y

    return(x)
}

