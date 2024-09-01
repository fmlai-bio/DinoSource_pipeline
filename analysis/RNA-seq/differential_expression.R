# setwd("D:/diffexpress/")

# 0、test
# countFile <- 'Ac_test_count.csv'
# controlSamples <- 'AC_NC_1'
# treatSamples <- 'AC_NS_3'
# lo2FCThreshold <- 1
# padjThreshold <- 0.05
# outputFile <- 'test2.csv'


# 
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




# 1、function

deseq2_func <- function(targets,lo2FCThreshold,padjThreshold,outputFile,a,b){
  library(DESeq2)
  input_data <- round(targets,digits = 0)

  input_data <- as.matrix(input_data)
  condition <- factor(c(rep("control",a),rep("treat",b)))
  coldata <- data.frame(row.names = colnames(input_data),condition)

  dds <- DESeqDataSetFromMatrix(countData=input_data,colData=coldata,design=~condition)

  dds <- DESeq(dds)

  res <- results(dds,alpha=0.1)
  summary(res)

  res <- res[order(res$padj),]
  resdata <- merge(as.data.frame(res),as.data.frame(counts(dds,normalized=TRUE)),by="row.names",sort=FALSE)
  names(resdata)[1] <- "Gene"
  resdata[which(resdata$log2FoldChange >= lo2FCThreshold & resdata$padj < padjThreshold),'significance'] <- 'up'
  resdata[which(resdata$log2FoldChange <= -lo2FCThreshold & resdata$padj < padjThreshold),'significance'] <- 'down'
  resdata[which(abs(resdata$log2FoldChange) <= lo2FCThreshold | resdata$padj >= padjThreshold),'significance'] <- 'none'
  

  result <- resdata[, c(1:7, which(colnames(resdata) == "significance"))]
  colnames(result) <- c("Gene","baseMean","log2FoldChange","lfcSE","stat","Pvalue","Padj","significance")
  result <- na.omit(result) 
  
  result$baseMean <- round(result$baseMean, 2)
  result$log2FoldChange <- round(result$log2FoldChange, 4)
  result$lfcSE <- round(result$lfcSE, 4)
  result$stat <- round(result$stat, 4)

  result$Pvalue <- formatC(result$Pvalue, format = "e", digits = 4)
  result$Padj <- formatC(result$Padj, format = "e", digits = 4)
  
  
  
  write.csv(result,file = outputFile, row.names = F,quote = F)
}



des_func <- function(targets,lo2FCThreshold,padjThreshold,outputFile){
  library(edgeR)
  exprSet <- data.frame(GeneId=row.names(targets),targets[,1],targets[,2])
  colnames(exprSet)<-c("GeneId",colnames(targets))

  group <- 1:2
  y <- DGEList(counts = exprSet[,2:3],genes = exprSet[,1],group = group)
  

  keep <- rowSums(cpm(y)>1) >= 1
  y <- y[keep, , keep.lib.sizes=FALSE]

  y <- calcNormFactors(y)
  

  y_bcv <- y
  bcv <- 0.1
  et <- exactTest(y_bcv, dispersion = bcv ^ 2)

  gene1 <- decideTestsDGE(et, p.value = padjThreshold, lfc = lo2FCThreshold)
  summary(gene1)

  colnames(gene1) <- "significance"

  resdata <- cbind(y$genes,y$counts,et$table,gene1)
  resdata$significance <- ifelse(resdata$significance == 1, "up", ifelse(resdata$significance == -1, "down", "none"))
  

  result <- result <- resdata[, c("genes", "logFC", "logCPM",'PValue','significance')]
  colnames(result) <- c("Gene","log2FoldChange","logCPM","Pvalue","significance")
  result <- na.omit(result) 
  
  # 
  result$log2FoldChange <- round(result$log2FoldChange, 4)
  result$logCPM <- round(result$logCPM, 4)
  
  # 
  result$Pvalue <- formatC(result$Pvalue, format = "e", digits = 4)

  write.csv(result,file = outputFile, row.names = F,quote = F)

}











# 2、Run diffexpression analysis
if (length(controlSamples_vector) > 1 & length(treatSamples_vector) >1 ){
  a <- length(controlSamples_vector)
  b <- length(treatSamples_vector)
  deseq2_func(targets,lo2FCThreshold,padjThreshold,outputFile,a,b)
}else{
  # edger_func(targets,group,lo2FCThreshold,padjThreshold,outputFile)
  des_func(targets,lo2FCThreshold,padjThreshold,outputFile)
}








