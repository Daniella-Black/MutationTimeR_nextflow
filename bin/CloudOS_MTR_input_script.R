#!/usr/local/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

sample = args[1]
txt <- "Hallo\nWorld"
writeLines(txt, paste0("out_", sample, ".txt")
