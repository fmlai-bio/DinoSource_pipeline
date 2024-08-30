library(limma)
library(WGCNA)
library(FactoMineR)
library(factoextra)  
library(tidyverse) # ggplot2 stringer dplyr tidyr readr purrr  tibble forcats
# library(data.table) 
enableWGCNAThreads(nThreads = 0.75*parallel::detectCores()) 


# usage
# Rscript WGCNA.R


args <- commandArgs(trailingOnly = TRUE) 
inputFile=args[1]
outputFile=args[2]
batchFile=args[3]
keepMAD=args[4]
tomWeight = args[5]

keepMAD <- as.numeric(keepMAD)
tomWeight <- as.numeric(tomWeight)

# keepMAD
# tomWeight





#fpkm <- read.table(inputFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,row.names =1)
fpkm <- read.csv(inputFile, header = TRUE, stringsAsFactors = FALSE,row.names =1)


batchMar <- read.table(batchFile,header = T)
batch <- batchMar[,2]

# batch<- c(rep('AC',6),rep('GEO',6))
fpkm <- removeBatchEffect(fpkm,batch = batch)

fpkm <- fpkm[!apply(fpkm, 1, function(row) any(row < 0)), ]



data <- log2(fpkm+1)

if(keepMAD == "all"){
  keep_data <- fpkm
}else{
  keep_data <- data[order(apply(data,1,mad), decreasing = T)[1:keepMAD],]
}

datExpr <- as.data.frame(t(keep_data))


gsg <- goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes],
                                              collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
}
gsg <- goodSamplesGenes(datExpr,verbose = 3)
gsg$allOK




R.sq_cutoff = 0.8  
powers <- c(seq(1,20,by = 1), seq(22,30,by = 2))    
sft <- pickSoftThreshold(datExpr, 
                         networkType = "unsigned",
                         powerVector = powers, 
                         RsquaredCut = R.sq_cutoff,  
                         verbose = 5)
#SFT.R.sq > 0.8 , slope ≈ -1
# # Plot the results: 
# par(mfrow = c(1,2));
# cex1 = 0.9;
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n")
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red")
# abline(h=R.sq_cutoff ,col="red")
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n")
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")



sft$powerEstimate  
wgcnapower = sft$powerEstimate
# wgcnapower = 11


if(is.na(wgcnapower)){
  type = "unsigned"
  nSamples=nrow(datExpr)
  wgcnapower = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                      ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                             ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                    ifelse(type == "unsigned", 6, 12))      
                      )
  )
}

promot <- '选择的power阈值是：'
promot
wgcnapower


net <- blockwiseModules(
  datExpr,
  power = wgcnapower,
  maxBlockSize = ncol(datExpr),
  corType = "pearson", #默认为"pearson","bicor"则更能考虑离群点的影响
  networkType = "unsigned",
  TOMType = "unsigned", 
  minModuleSize = 5,    ##越大模块越少
  mergeCutHeight = 0, ##越大模块越少
  numericLabels = TRUE, 
  saveTOMs = F,
  verbose = 3
)




# IMConn = softConnectivity(datExpr) 


corMatrix <- abs(cor(datExpr, use = "p"))


TOM <- TOMsimilarityFromExpr(datExpr, power = wgcnapower)


# highWeightMatrix <- TOM >= 0.25
highWeightMatrix <- TOM >= tomWeight


edgeIndices <- which(highWeightMatrix & !diag(nrow(highWeightMatrix)), arr.ind = TRUE)


geneNames <- colnames(datExpr)


edgeList <- data.frame(
  sourceGene = geneNames[edgeIndices[, 2]],
  targetGene = geneNames[edgeIndices[, 1]],
  weight = TOM[edgeIndices],
  pearson = corMatrix[edgeIndices]
)


edgeList <- edgeList %>%
  arrange(sourceGene, desc(weight), desc(pearson))



moduleColors <- net$colors
names(moduleColors) <- colnames(datExpr)


edgeList$module1 <- moduleColors[edgeList$sourceGene]
edgeList$module2 <- moduleColors[edgeList$targetGene]


edgeList <- edgeList[edgeList$module1 == edgeList$module2, ]


edgeList$module <-  paste0('module',edgeList$module1)

edgeList <- edgeList[, c("sourceGene", "targetGene", "weight", "pearson",'module')]




edgeList <- edgeList %>% mutate(id = 1:n())


edgeList <- edgeList %>% relocate(id, .before = sourceGene)




write.csv(edgeList, file = outputFile, row.names = FALSE,quote = F)
