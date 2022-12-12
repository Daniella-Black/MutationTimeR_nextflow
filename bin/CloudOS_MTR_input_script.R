#!/usr/local/bin/Rscript
library(signature.tools.lib)
#library(MutationTimeR)
#library(mg14)
args = commandArgs(trailingOnly=TRUE)


sampleID <- args[1]
vcfpath <- args[2]

genome.v="hg38"
PR_threshold=8
# outdir
organ <- 'Breast'
smallvariants_VCF <- VariantAnnotation::readVcf(vcfpath,genome = genome.v)
e.vcf <- VariantAnnotation::expand(smallvariants_VCF)

# separate SNV and Indels
rd <- SummarizedExperiment::rowRanges(e.vcf)
e.snv <- e.vcf[nchar(as.character(rd$REF))==1 & nchar(as.character(rd$ALT))==1,]
selected_snv <- VariantAnnotation::fixed(e.snv)[,"FILTER"]=="PASS" & as.character(SummarizedExperiment::seqnames(e.snv)) %in% paste0("chr",c(1:22,"X","Y"))
e.snv <- e.snv[selected_snv,]


# collect SNV data into a table
rd <- SummarizedExperiment::rowRanges(e.snv)
rgs <- IRanges::ranges(e.snv)
vafs <- c(VariantAnnotation::info(e.snv)[,"VAF"])
pdf(file = paste0(sampleID, '_vaf_hist_all_muts.pdf'))
hist(vafs)
dev.off()
