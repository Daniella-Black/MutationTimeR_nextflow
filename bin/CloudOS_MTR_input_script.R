#!/usr/local/bin/Rscript
args = commandArgs(trailingOnly=TRUE)


sample <- args[1]
cnvpath <- args[2]
vcfpath <- args[3]

callsBattenberg <- read.csv(cnvpath, sep='\t')
vcf <- read.Vcf(vcfpath)

bb <- makeGRangesFromDataFrame(callsBattenberg, keep.extra.columns=TRUE)
bb
mt <- mutationTime(vcf, bb, n.boot=100)

#######
head(mt$V)
table(mt$V$CLS)
vcf <- addMutTime(vcf, mt$V)
head(mt$T)
mcols(bb) <- cbind(mcols(bb),mt$T)

pdf(file=paste0(sample, '_vaf_expected_vaf.pdf'), height=7.5, width=7.8)
plotSample(vcf,bb)
dev.off()

mt$T <- mt$T[-c(1)]

write.csv(mt$T, paste0(sample,'_mT.csv'))
write.csv(mt$V, paste0(sample,'_mV.csv'))
write.csv(mt$V$CLS, paste0(sample,'_CLS.csv'))
