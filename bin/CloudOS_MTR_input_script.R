#!/usr/local/bin/Rscript
library(signature.tools.lib)
library(stringr)
#library(MutationTimeR)
#library(mg14)
args = commandArgs(trailingOnly=TRUE)


sampleID <- args[1]
vcfpath <- args[2]
cnvpath <- args[3]
tp <- args[4]
header <- args[5]


genome.v="hg38"
PR_threshold=8
# outdir
organ <- 'Breast'
header <- read.csv(header, sep='\n', header=FALSE)

#read in the vcf as a table to reduce memory usage compraed to VariantAnnotation::readVCf
snvtab <- read.table(vcfpath, sep='\t')

col_order <- c("chr", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO","FORMAT" , 'TUMOR')  
##change the colnames 
#if(length(colnames(snvtab)) == 11){
#  col_order <- c("chr", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO","FORMAT" , 'TUMOR' , 'NORMAL')
#} 
names(snvtab) <- col_order

##select only snvs
snvtab <- snvtab[nchar(as.character(snvtab$REF))==1 & nchar(as.character(snvtab$ALT))==1,]

##select variants which passed the filter
snvtab <- subset(snvtab, FILTER=='PASS')
##split TUMOR to get required values
snvtab[c('DP', 'FDP', 'SDP', 'SUBDP', 'AU', 'CU', 'GU', 'TU')] <- str_split_fixed(snvtab$TUMOR, ':', 8)
snvtab[c('A', 'AU_tier2')] <- str_split_fixed(snvtab$AU, ',',2)
snvtab[c('T', 'TU_tier2')] <- str_split_fixed(snvtab$TU, ',',2)
snvtab[c('C', 'CU_tier2')] <- str_split_fixed(snvtab$CU, ',',2)
snvtab[c('G', 'GU_tier2')] <- str_split_fixed(snvtab$GU, ',',2)

##getting the AD and DP values
snvtab['TUMOR'] <- NULL
tumour_list <- list()

snvtab$chr = gsub("chr", "", snvtab$chr)

snvtab['AD_REF'] <- rep(0, nrow(snvtab))
snvtab['AD_ALT'] <- rep(0, nrow(snvtab))
snvtab['FORMAT'] <- rep('DP:AD', nrow(snvtab))
for(i in 1:2){
  cols = c('QUAL', 'INFO')
  i = cols[i]
  snvtab[i] <- rep('.', nrow(snvtab))
}

write.table(snvtab,file = paste0(sampleID,"_SNVs.txt"),sep = "\t",quote = F,col.names = T,row.names = F)

#for (row in 1:nrow(snvtab)){
#  bases <- c('A', 'T', 'C', 'G')
#  for(nuc in 1:4){
#     nuc <- bases[nuc]
#     if(snvtab$REF[row] == nuc){
#        snvtab$AD_REF[row] <-  snvtab$AD_REF[row] + as.integer(snvtab[row, nuc])
#      }
#      else if(snvtab$ALT[row] == nuc){
#        snvtab$AD_ALT[row] <- snvtab$AD_ALT[row] +  as.integer(snvtab[row, nuc])
#     }
#    }
#  tumour_list <- append(tumour_list, paste(snvtab$DP[row], ':', snvtab$AD_REF[row], ',', snvtab$AD_ALT[row], sep=''))
#}

#snvtab['TUMOR'] <- tumour_list

#names(snvtab)[names(snvtab) =='chr'] <- '#CHROM'
#col_order <- c("#CHROM", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO","FORMAT" , 'TUMOR' )
#snvtab <- snvtab[, col_order]

#write.table(header, file = paste0(sampleID,"_SNVs.txt"),row.names = F,quote = F,sep = "\t", col.names=F)
#write.table(snvtab,file = paste0(sampleID,"_SNVs.txt"),sep = "\t",quote = F,col.names = T,row.names = F, append=T)

########################################################################################
##process cnv file
########################################################################################

#read in  the file
sv_vcf <- VariantAnnotation::readVcf(file = cnvpath, genome = genome.v)
# filter PASS and chromosomes
selected_sv <- VariantAnnotation::fixed(sv_vcf)[,"FILTER"]=="PASS" | VariantAnnotation::fixed(sv_vcf)[,"FILTER"]=="MGE10kb"
sv_vcf <- sv_vcf[selected_sv,]
select_chrom <- as.character(GenomeInfoDb::seqnames(sv_vcf)) %in% paste0("chr",c(1:22,"X","Y"))
sv_vcf <- sv_vcf[select_chrom,]

rr <- SummarizedExperiment::rowRanges(sv_vcf)

#select Canvas rows only
selection <- grepl(names(rr),pattern = "^Canvas")
cn_vcf <- sv_vcf[selection,]

rr <- SummarizedExperiment::rowRanges(cn_vcf)

#select only the copy number
sv_size <- VariantAnnotation::info(cn_vcf)$END - BiocGenerics::start(rr)
selection <- sv_size >= 10000
cn_vcf <- cn_vcf[selection,]


#construct cn df
chr <- sapply(as.character(GenomeInfoDb::seqnames(cn_vcf)), function(x) ifelse(grepl(x,pattern = "^chr"),substr(x,start = 4,stop = 5),x))
rr <- SummarizedExperiment::rowRanges(cn_vcf)
#sample_MCC <- VariantAnnotation::geno(cn_vcf)$MCC[,1]
#sample_MCC[is.na(sample_MCC)] <- VariantAnnotation::geno(cn_vcf)$CN[is.na(sample_MCC),1]
sample_MCC <- VariantAnnotation::geno(cn_vcf)$MCC[,2]
sample_MCC[is.na(sample_MCC)] <- VariantAnnotation::geno(cn_vcf)$CN[is.na(sample_MCC),2]
sv_df <- data.frame(
             seqnames = chr,
             start = BiocGenerics::start(rr),
             end = VariantAnnotation::info(cn_vcf)$END,
             #total.copy.number.inTumour = VariantAnnotation::geno(cn_vcf)$CN[,1],
             #minor_cn = VariantAnnotation::geno(cn_vcf)$CN[,1] - sample_MCC,stringsAsFactors = FALSE)
             total.copy.number.inTumour = VariantAnnotation::geno(cn_vcf)$CN[,2],
             minor_cn = VariantAnnotation::geno(cn_vcf)$CN[,2] - sample_MCC,
             stringsAsFactors = FALSE)
            

sv_df$start = sv_df$start +1
sv_df$width = sv_df$end - sv_df$start
sv_df$major_cn = sv_df$total.copy.number.inTumour - sv_df$minor_cn
sv_df['strand'] <- rep('*', nrow(sv_df))
sv_df['clonal_frequency'] <- rep(as.numeric(tp)/100, nrow(sv_df))
sv_df <- sv_df[ , -which(names(sv_df) %in% c("total.copy.number.inTumour"))]
col_order <- c('seqnames', 'start', 'end', 'width', 'strand', 'major_cn', 'minor_cn', 'clonal_frequency')
sv_df <- sv_df[, col_order]

write.table(sv_df,file = paste0(sampleID,"_CNVs.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)

sv_df <- NULL
sample_MCC <-NULL
rr <- NULL
sv_size <- NULL
selection <- NULL
cn_vcf <- NULL
snvtab <- NULL
rgs <-NULL
rd <- NULL
smallvariants_VCF <- NULL
e.vcf <- NULL
e.smv <- NULL
