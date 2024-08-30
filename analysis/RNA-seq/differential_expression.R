#setwd("D:/diffexpress/")
library(edgeR)
library(tibble)


# countFile <- 'Ac_test_count.csv'
# controlSamples <- 'AC_NC_1,AC_NC_2,AC_NC_3'
# treatSamples <- 'AC_NS_1,AC_NS_2,AC_NS_3'
# lo2FCThreshold <- 1
# padjThreshold <- 0.05
# outputFile <- 'test.csv'



args <- commandArgs(trailingOnly = TRUE)
countFile <- args[1]
lo2FCThreshold <- as.numeric(args[2])
padjThreshold <- as.numeric(args[3])
controlSamples <- args[4]
treatSamples <- args[5]
outputFile <- args[6]


controlSamples_vector <- unlist(strsplit(controlSamples, split = ","))
treatSamples_vector <- unlist(strsplit(treatSamples, split = ","))
lo2FCThreshold <- as.numeric(lo2FCThreshold)
padjThreshold <- as.numeric(padjThreshold)

targets <- read.csv(countFile, row.names = 1)


# group <- rep(c('control', 'treat'), each = 3)
group <- c(rep('control',length(controlSamples_vector)),rep('treat',length(treatSamples_vector)))



dgelist <- DGEList(counts = targets, group = group)

keep <- rowSums(cpm(dgelist) > 1 ) >= 2
dgelist <- dgelist[keep, , keep.lib.sizes = FALSE]

dgelist_norm <- calcNormFactors(dgelist, method = 'TMM')


design <- model.matrix(~group)

dge <- estimateDisp(dgelist_norm, design, robust = TRUE)



if( length(group) > 2 ){

  #quasi-likelihood(QL) F-test tests:  bulk RNA-seq
  fit <- glmQLFit(dge, design, robust = TRUE)
  lt <- glmQLFTest(fit)
  lrt <- topTags(lt, n = nrow(dgelist$counts))
}else{
  #likelihood ratio test：scRNA-seq and no replicates data
  fit <- glmFit(dge, design, robust = TRUE)
  lt <- glmLRT(fit)
  lrt <- topTags(lt, n = nrow(dgelist$counts))
}



resdata <- lrt$table[,c(1,2,4,5)]  
# colnames(resdata) <- c('log2FoldChange','logCPM','F','pvalue','padj')
colnames(resdata) <- c('log2FoldChange','logCPM','Pvalue','Padj')

resdata[which(resdata$log2FoldChange >= lo2FCThreshold & resdata$Padj < padjThreshold),'significance'] <- 'up'
resdata[which(resdata$log2FoldChange <= -lo2FCThreshold & resdata$Padj < padjThreshold),'significance'] <- 'down'
resdata[which(abs(resdata$log2FoldChange) <= lo2FCThreshold | resdata$Padj >= padjThreshold),'significance'] <- 'none'
# filtered_data <- resdata[resdata$sig %in% c("up", "down"), ]

resdata <- rownames_to_column(resdata, var = "Gene")




write.csv(resdata,file = outputFile, row.names = F,quote = F)












