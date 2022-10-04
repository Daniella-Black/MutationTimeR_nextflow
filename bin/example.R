#!/usr/local/bin/Rscript
args = commandArgs(trailingOnly=TRUE)


sample <- args[1]
cnvpath <- args[2]
vcfpath <- args[3]

txt <- paste0(sample, cnvpath, vcfpath)
writeLines(txt, "outfile.txt")

