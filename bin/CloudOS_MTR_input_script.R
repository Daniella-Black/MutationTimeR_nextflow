#!/usr/local/bin/Rscript
library(signature.tools.lib)
args = commandArgs(trailingOnly=TRUE)


sampleID <- args[1]
vcfpath <- args[2]
cnvpath <- args[3]
tp <- args[3]

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

snvtab['VAF'] <- snvtab['AD_ALT']/(snvtab['AD_ALT'] + snvtab['AD_REF'])

snvtab <- snvtab[ , -which(names(snvtab) %in% c("DP","A", "T", "G", "C", "AD_REF", "AD_ALT"))]

col_order <- c("#CHROM", "POS", "ID", "REF","ALT", "QUAL", "FILTER", "INFO","FORMAT" , 'TUMOR' )
snvtab <- snvtab[, col_order]

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
            

sv_df$start = sv_df$start +1
sv_df$width = sv_df$end - sv_df$start
sv_df$major_cn = sv_df$total.copy.number.inTumour - sv_df$minor_cn
sv_df['strand'] <- rep('*', nrow(sv_df))
sv_df['clonal_frequency'] <- rep(as.numeric(tp)/100, nrow(sv_df))
sv_df <- sv_df[ , -which(names(sv_df) %in% c("total.copy.number.inTumour"))]
col_order <- c('seqnames', 'start', 'end', 'width', 'strand', 'major_cn', 'minor_cn', 'clonal_frequency')
sv_df <- sv_df[, col_order]

################################################################################################              
#####select mutations in regions of normal cn only
################################################################################################
              
#select only the copy number normal regions
cn_normal <- subset(sv_df, major_cn ==1 & minor_cn == 1)
#set up some lists to fill or iterate through
chr <- c('1', '2', '3', '4','5', '6', '7', '8','9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y')
muts_normal <- c()
for(chrom in chr){
  #to avoid iterating over all chromosomes do just one at a time - subset the vcf by mutations on the chrom of interest
  vcf_per_chrom = subset(snvtab, snvtab['#CHROM'] == chrom)
  #do the same with the cn
  cn_norm_per_chrom = subset(cn_normal, seqnames == chrom)
  if(nrow(vcf_per_chrom > 0)){
    for(mut in 1:nrow(vcf_per_chrom)){
      for(contig in 1:nrow(cn_norm_per_chrom)){
            if(vcf_per_chrom[mut, 'POS'] > cn_norm_per_chrom[contig, 'start'] & vcf_per_chrom[mut, 'POS'] <cn_norm_per_chrom[contig, 'end']){
              muts_normal <- append(muts_normal,vcf_per_chrom[mut, 'ID']) ##append the ID
          }
        }    
      }
   }
}

snvtab_normal <- snvtab[snvtab$V3 %in% muts_normal,]    

#######################################################################################################
##make the VAF histograms##
#######################################################################################################
#fileConn<-file("output.txt")
#writeLines(c(as.character(nrow(snvtab))), fileConn)
#close(fileConn)
              

write.table(snvtab,file = paste0(sampleID,"_SNVs.txt"),sep = "\t",quote = F,col.names = T,row.names = F)

write.table(snvtab_normal,file = paste0(sampleID,"_CNVs.tsv"),sep = "\t",quote = F,col.names = T,row.names = F)
             
bins=seq(0,1.0,by=0.01)

#pdf(file = paste0(sampleID, '_vaf_hist_all_muts.pdf'))
#hist(as.numeric(unlist(snvtab['VAF'])), breaks=bins, main = paste0(sampleID, ' (tumour purity = ', tp, ')'), xlab='VAF', col='#fadadd')
#dev.off()

#pdf(file = paste0(sampleID, '_vaf_hist_normal_cn.pdf'))
#hist(as.numeric(unlist(snvtab_normal['VAF'])), breaks=bins, main = paste0(sampleID, ' (tumour purity = ', tp, '), ',  '\n mutations in diploid regions = ', nrow(snvtab_normal), '/', nrow(snvtab)), 
#     xlab='VAF', col='#fadadd')
#dev.off()
