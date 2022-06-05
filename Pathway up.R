library(dplyr)
dat = read.csv("input.csv")
#将数据按要求进行筛选，根据具体组别设置决定这里选择的是logFC<0或者>0的数据
#p<0.05
subDatup_05<- filter(dat, adj.P.Val <= 0.05 & SYMBOL != "NA"& logFC > "2")
write.table(subDatup_05,"pathway up.txt",row.names =FALSE)
test1<-subDatup_05[,c(1,4,5)]
sampleNames <- c("SYMBOL","ENSEMBL", "ENTREZID")
names(test1)[1:3] <- sampleNames
head(test1)

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

#运行GO分析，单独组别
ego_MF <- enrichGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db,ont = "MF", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.5,
                   keyType = 'ENTREZID', readable = TRUE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
write.csv(summary(ego_MF),"GO MF-enrich up.csv",row.names =FALSE)
dotplot(ego_MF,showCategory=20,title="EnrichmentGO_MF_dot up")
barplot(ego_MF, showCategory=20,title="EnrichmentGO_MF_bar up")

ego_BP <- enrichGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db,ont = "BP", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.5,
                   keyType = 'ENTREZID', readable = TRUE)
ego_BP1 <- setReadable(ego_BP, OrgDb = org.Hs.eg.db)
write.csv(summary(ego_BP),"GO BP-enrich up.csv",row.names =FALSE)
dotplot(ego_BP,showCategory=20,title="EnrichmentGO_BP_dot up")
barplot(ego_BP, showCategory=20,title="EnrichmentGO_BP_bar up")

ego_CC <- enrichGO(gene = test1$ENTREZID, OrgDb = org.Hs.eg.db,ont = "CC", 
                   pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.5,
                   keyType = 'ENTREZID', readable = TRUE)
ego_CC1 <- setReadable(ego_CC, OrgDb = org.Hs.eg.db)
write.csv(summary(ego_CC),"GO CC-enrich up.csv",row.names =FALSE)
dotplot(ego_CC,showCategory=20,title="EnrichmentGO_CC_dot up")
barplot(ego_CC, showCategory=20,title="EnrichmentGO_CC_bar up")


#运行KEGG分析
library(stringr)
kk <- enrichKEGG(gene = test1$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05,
                 pAdjustMethod = 'BH',
                 keyType = 'kegg',
                 minGSSize = 3,
                 maxGSSize = 500,
                 qvalueCutoff = 0.05,
                 use_internal_data = FALSE)
head(kk,2)
write.csv(summary(kk),"KEGG-enrich up.csv",row.names =FALSE)
dotplot(kk,showCategory=20,title="Enrichment KEGG_dot up")
barplot(kk,showCategory=20,title="Enrichment KEGG_bar up")
