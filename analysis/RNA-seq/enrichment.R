
library(clusterProfiler)     
library(GO.db)
library(tidyr)


args <- commandArgs(trailingOnly = TRUE)  

# geneList <- unlist(strsplit(args[1], ","))
geneFile <- args[1]
backgroundFile=args[2]
type=args[3]
pvalue=as.numeric(args[4])
qvalue=as.numeric(args[5])
outputFile=args[6]
minCount=as.numeric(args[7])

geneList <- read.csv(geneFile,header = F)[,1] 

#setwd("/home/lfm/DinosourceWeb/enrichGOKEGG/")
# geneList <- c("Acart01G0064","Acart01G0382","Acart01G0764","Acart01G0821","Acart01G1595","Acart01G1597","Acart01G1646","Acart02G0188","Acart03G0433","Acart03G0819","Acart03G0886","Acart03G1149","Acart03G1181","Acart03G1184","Acart04G0619","Acart04G0684","Acart05G0327","Acart05G1316","Acart06G0676","Acart06G0700","Acart06G0929","Acart06G0932","Acart07G0367","Acart07G0369","Acart07G0647","Acart08G0567","Acart09G0087","Acart09G0197","Acart09G0847")
#geneList <- c("Acart03G1149","Acart03G1181","Acart03G1184","Acart06G0929","Acart16G0031","Acart16G0929","Acart20G0387","Acart20G0388","Acart24G0088","Acart30G0642","Acart39G0500","Acart40G0267","Acart40G0269","Acart41G0301")
#backgroundFile <- '/home/lfm/DinosourceWeb/enrichGOKEGG/public/annotation.GO.KEGG.clean.txt'
#type <- "kegg"
#pvalue <- as.numeric('0.05')
#qvalue <- as.numeric('0.05')

#outputFile <- '/home/lfm/DinosourceWeb/enrichGOKEGG/results/test.20240719.kegg.csv'
# minCount <- 1


# pubilc file
koidKidFile <- 'D:/enrichGOKEGG/public/koid-kid.txt'
koInfoFile <- 'D:/enrichGOKEGG/public/ko_info.txt'





goterm <- Term(GOTERM)
golist <- as.data.frame(goterm)

background1 <- read.table(backgroundFile,header = T)
gene2go2kegg <- background1[,c(1,3,2)]
colnames(gene2go2kegg) <- c('query','GOs','KEGG_ko')
go2name <- data.frame(Go_id=rownames(golist),Description=golist$goterm)
go2gene <- data.frame(Go_id=gene2go2kegg$GOs,Gene=gene2go2kegg$query)

go2gene<-go2gene[go2gene$Go_id!="-",]

go2gene <- go2gene %>% separate_rows(Go_id, sep = ",")
go2gene <- unique(go2gene) 
#table(duplicated(go2gene))  
diffgo <- setdiff(go2gene$Go_id,go2name$Go_id)
length(diffgo)

go2gene <- subset(go2gene, !(Go_id %in% diffgo))

diffgo <- setdiff(go2gene$Go_id,go2name$Go_id)
length(diffgo)



kid2gene <- data.frame(Kegg_id=gene2go2kegg$KEGG_ko,Gene=gene2go2kegg$query)
kid2gene<-kid2gene[kid2gene$Kegg_id!="-",]

kid2gene <- kid2gene %>% separate_rows(Kegg_id, sep = ",")
kid2gene <- unique(kid2gene)
#table(duplicated(kid2gene)) 
koid2kid <- read.table(koidKidFile,header = T)
merged_data <- merge(kid2gene, koid2kid, by.x = "Kegg_id", by.y = "Kid", all.x = F)
kegg2gene <- merged_data[,c(3,2)]
kegg2name <- read.table(koInfoFile,header = T,sep = "\t")

# save(go2gene,go2name,kegg2gene,kegg2name,file="./results/GO_KEGG.RData")




#load("./results/GO_KEGG.RData")
enrich_func_go <- function(gene_set,outputFile,pvalue,qvalue,minCount){
  GO <- enricher(gene_set,
                 TERM2GENE=go2gene,
                 TERM2NAME=go2name,
                 pvalueCutoff = pvalue, 
                 qvalueCutoff = qvalue,
                 minGSSize = minCount)
  if ( is.null(GO)||is.null(GO@result)){
    x <- "there is no gene is enriched to GO."
    print(x)

    #write.csv(x,file = outputFile,row.names = F)
  }else{
    go <- GO@result   
    go$Count <- as.numeric(go$Count)
    go <-  go[go$Count >= minCount, ]
    write.csv(go,file = outputFile,row.names = F)
  }
}

enrich_func_kegg <- function(gene_set,outputFile,pvalue,qvalue,minCount){
  KEGG <- enricher(gene_set,
                   TERM2GENE=kegg2gene,
                   TERM2NAME=kegg2name,
                   pvalueCutoff = pvalue, 
                   qvalueCutoff = qvalue,
                   minGSSize = minCount)
  if (is.null(KEGG) || is.null(KEGG@result)){
    x <- "there is no gene is enriched to KEGG."
    print(x)
    #write.csv(x,file = outputFile,row.names = F)
  }else{
  
    kegg <- KEGG@result       
    kegg$Count <- as.numeric(kegg$Count)
    kegg <- kegg[kegg$Count >= minCount, ]
    write.csv(kegg,file = outputFile,row.names = F)
  }
}





if(type == 'go'){
  enrich_func_go(gene_set = geneList,outputFile = outputFile,pvalue = pvalue,qvalue = qvalue,minCount=minCount)  
}else if (type == 'kegg'){
  enrich_func_kegg(gene_set = geneList,outputFile = outputFile,pvalue = pvalue,qvalue = qvalue,minCount=minCount)  
}











