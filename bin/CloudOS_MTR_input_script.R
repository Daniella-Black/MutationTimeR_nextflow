#!/usr/local/bin/Rscript
library(signature.tools.lib)
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
                     #VAF=VariantAnnotation::info(e.snv)[,"VAF"],
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
snvtab['FORMAT'] <- rep('DP:AD', nrow(snvtab))
for(i in 1:2){
  cols = c('QUAL', 'INFO')
  i = cols[i]
  snvtab[i] <- rep('.', nrow(snvtab))
}


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
  tumour_list <- append(tumour_list, paste(snvtab$DP[row], ':', snvtab$AD_REF[row], ',', snvtab$AD_ALT[row], sep=''))
    }

snvtab['TUMOR'] <- tumour_list

snvtab <- snvtab[ , -which(names(snvtab) %in% c("DP","A", "T", "G", "C", "AD_REF", "AD_ALT"))]

col_order <- c("#CHROM", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO","FORMAT" , 'TUMOR' )
snvtab <- snvtab[, col_order]
write.table(header, file = paste0(sampleID,"_SNVs.txt"),row.names = F,quote = F,sep = "\t", col.names=F)
write.table(snvtab,file = paste0(sampleID,"_SNVs.txt"),sep = "\t",quote = F,col.names = T,row.names = F, append=T)

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
sample_MCC <- VariantAnnotation::geno(cn_vcf)$MCC[,1]
sample_MCC[is.na(sample_MCC)] <- VariantAnnotation::geno(cn_vcf)$CN[is.na(sample_MCC),1]
#sample_MCC <- VariantAnnotation::geno(cn_vcf)$MCC[,2]
#sample_MCC[is.na(sample_MCC)] <- VariantAnnotation::geno(cn_vcf)$CN[is.na(sample_MCC),2]
sv_df <- data.frame(
             seqnames = chr,
             start = BiocGenerics::start(rr),
             end = VariantAnnotation::info(cn_vcf)$END,
             total.copy.number.inTumour = VariantAnnotation::geno(cn_vcf)$CN[,1],
             minor_cn = VariantAnnotation::geno(cn_vcf)$CN[,1] - sample_MCC,stringsAsFactors = FALSE)
             #total.copy.number.inTumour = VariantAnnotation::geno(cn_vcf)$CN[,2],
             #minor_cn = VariantAnnotation::geno(cn_vcf)$CN[,2] - sample_MCC,
             #stringsAsFactors = FALSE)
            

#sv_df$start = sv_df$start +1
#sv_df$width = sv_df$end - sv_df$start
#sv_df$major_cn = sv_df$total.copy.number.inTumour - sv_df$minor_cn
#sv_df['strand'] <- rep('*', nrow(sv_df))
#sv_df['clonal_frequency'] <- rep(as.numeric(tp)/100, nrow(sv_df))
#sv_df <- sv_df[ , -which(names(sv_df) %in% c("total.copy.number.inTumour"))]
#col_order <- c('seqnames', 'start', 'end', 'width', 'strand', 'major_cn', 'minor_cn', 'clonal_frequency')
#sv_df <- sv_df[, col_order]

#write.table(sv_df,file = paste0(sampleID,"_CNVs.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
