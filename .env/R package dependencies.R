# System
options("repos" = c(CRAN = "https://cran.revolutionanalytics.com"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref="7ea9440", upgrade="never")
