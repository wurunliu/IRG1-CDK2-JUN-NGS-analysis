library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library(enrichplot)

sig.gene<-read.csv(file="subDat.csv")
genelist<-sig.gene$log2FoldChange
names(genelist)<-as.character(sig.gene$ENTREZID)
genelist<-sort(genelist,decreasing=TRUE)
kk2<- gseKEGG(geneList     = genelist,
              organism     = 'hsa',
              nPerm        = 1000,
              minGSSize    = 10,
              pvalueCutoff = 0.5,
              verbose      = FALSE)
head(kk2)
write.csv(summary(kk2),"KEGG-GSEA.csv",row.names =FALSE)
gseaplot(kk2,geneSetID = "hsa04668")
gseaplot2(kk2, geneSetID = "hsa04668",color = "firebrick",pvalue_table=T)

hsa04510<- pathview(gene.data = genelist,
                    pathway.id = "hsa04510",
                    species = "hsa",
                    limit = list(gene=max(abs(genelist)), cpd=1))

#GO GSEA
geneList2<-sig.gene$log2FoldChange
names(geneList2)<-as.character(sig.gene$ENSEMBL)
geneList2<-sort(geneList2,decreasing=TRUE)
gsemf <- gseGO(geneList2, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont="MF", minGSSize = 30)
head(gsemf)
gsemf1 <- setReadable(gsemf, OrgDb = org.Hs.eg.db)
write.csv(summary(gsemf1),"GO MF-GSEA.csv",row.names =FALSE)
gseaplot(gsemf, geneSetID="GO:0001228")
gseaplot2(gsemf, geneSetID="GO:0007249",pvalue_table=T)


gsebp <- gseGO(geneList2, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont="BP", minGSSize = 30)
head(gsebp)
gsebp1 <- setReadable(gsebp, OrgDb = org.Hs.eg.db)
write.csv(summary(gsebp1),"GO BP-GSEA.csv",row.names =FALSE)
gseaplot(gsebp, geneSetID="GO:0007254")
gseaplot2(gsebp, geneSetID="GO:0007254",pvalue_table=T)


gsecc <- gseGO(geneList2, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont="CC", minGSSize = 30)
head(gsecc)
gsecc1 <- setReadable(gsecc, OrgDb = org.Hs.eg.db)
write.csv(summary(gsecc1),"GO CC-GSEA.csv",row.names =FALSE)
gseaplot(gsecc, geneSetID="GO:0002253")
gseaplot2(gsecc, geneSetID="GO:0012507",pvalue_table=T)


