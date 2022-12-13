#!/usr/local/bin/Rscript
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)


sampleID <- args[1]
vcfpath <- args[2]

genome.v="hg38"
PR_threshold=8
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
snvtab <- data.frame(chr=as.character(SummarizedExperiment::seqnames(e.snv)),
                     POS=BiocGenerics::start(rgs),
                     REF=as.character(rd$REF),
                     ALT=as.character(rd$ALT),
                     FILTER = VariantAnnotation::fixed(e.snv)[,"FILTER"],
                     DP=VariantAnnotation::geno(e.snv)$DP,
                     AU=VariantAnnotation::geno(e.snv)$AU,
                     TU=VariantAnnotation::geno(e.snv)$TU,
                     CU=VariantAnnotation::geno(e.snv)$CU,
                     GU=VariantAnnotation::geno(e.snv)$GU)

##makes the index column an ID column
snvtab <- cbind(ID = rownames(snvtab), snvtab)
rownames(snvtab) <- 1:nrow(snvtab)

snvtab$chr = gsub("chr", "", snvtab$chr)
sampleID_cols <- gsub("-", ".", sampleID)
snvtab = subset(snvtab, select = -c(9,11,13,15) )
names(snvtab)[names(snvtab) == sampleID_cols]<- "DP"
names(snvtab)[names(snvtab) == paste("AU.", sampleID_cols, ".1",sep='')] <- "A"
names(snvtab)[names(snvtab) == paste("TU.", sampleID_cols, ".1", sep='')] <- "T"
names(snvtab)[names(snvtab) == paste("CU.", sampleID_cols, ".1", sep='')] <- "C"
names(snvtab)[names(snvtab) == paste("GU.", sampleID_cols,  ".1", sep='')] <- "G"
names(snvtab)[names(snvtab) == 'chr'] <- '#CHROM'

snvtab['AD_REF'] <- rep(0, nrow(snvtab))
snvtab['AD_ALT'] <- rep(0, nrow(snvtab))

tumour_list<- c()

for (row in 1:nrow(snvtab)){
  bases <- c('A', 'T', 'C', 'G')
  for(nuc in 1:4){
    nuc <- bases[nuc]
    if(snvtab$REF[row] == nuc){
      snvtab$AD_REF[row] <-  snvtab$AD_REF[row] + as.integer(snvtab[row, nuc])
    }
    else if(snvtab$ALT[row] == nuc){
      snvtab$AD_ALT[row] <- snvtab$AD_ALT[row] +  as.integer(snvtab[row, nuc])
    }
  }
    }

snvtab['VAF'] <- snvtab['AD_ALT']/(snvtab['AD_ALT'] + snvtab['AD_REF'])

pdf(file = paste0(sampleID, '_vaf_hist_all_muts.pdf'))
hist(as.numeric(snvtab['VAF']))
dev.off()
