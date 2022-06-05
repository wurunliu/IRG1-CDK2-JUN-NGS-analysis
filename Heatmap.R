library(DESeq2)
mycounts<- read.table("mergedcounts M0 M3 C3 J3.txt", header = TRUE, quote = '\t',skip =0)
sampleNames <- c("M01","M02","M03","M31","M32","M33","C31","C32","C33","J31","J32","J33")
names(mycounts)[2:13] <- sampleNames
head(mycounts)

rownames(mycounts)<-mycounts[,1] 
mycounts<-mycounts[,-1] 
head(mycounts)

#变成矩阵
countMatrix <- as.matrix(mycounts[1:12])
table2 <- data.frame(name = c("M01","M02","M03","M31","M32","M33","C31","C32","C33","J31","J32","J33"),
                     condition = c("M0","M0","M0","M3","M3","M3","C3","C3","C3","J3","J3","J3"))
rownames(table2) <- sampleNames
head(countMatrix)

#把counts矩阵转化为DESeq2的数据格式
dds <- DESeqDataSetFromMatrix(countMatrix, colData=table2, design= ~ condition)
#过滤掉counts<0的数值
dds <- dds[rowSums(counts(dds)) >10,]

dds_norm <- DESeq(dds) #dds标准化
dds_norm #查看标准化后的数据

#获取标准化后的数据
normalized_counts <- counts(dds_norm, normalized=TRUE)
head(normalized_counts)

#根据基因在不同的样本中表达变化的差异程度mad值对数据排序，差异越大的基因排位越前
normalized_counts_mad <- apply(normalized_counts, 1, mad)
normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
head(normalized_counts)
# 标准化后的数据输出
write.table(normalized_counts, file="dds_normalized_counts.csv",quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'","dds_normalized_counts.csv"))
# 标准化后的数据取输出
rld<- rlog(dds_norm, blind=FALSE)
rlgMat<- assay(rld)
rlgMat<-rlgMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlgMat, file="dds_normalized_counts_rlog.xls",quote=F, sep="\t", row.names=T, col.names=T)
system(paste("sed -i '1 s/^/ID\t/'","dds_normalized_counts_rlog.xls"))

#火山图
library('pheatmap')
select <- order(rowMeans(normalized_counts),decreasing=T)
ntd <- normTransform(dds_norm) #default to log2FC+1
log2.norm.counts <- assay(ntd)[select,]
df <- as.data.frame(colData(dds_norm)[,c("name","condition")])
out<-pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F,
         cluster_cols=F, annotation_col=df,scale="row",
         treeheight_col=0,clustering_method = "average",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(491),
         clustering_distance_rows = "maximum",cutree_cols = 2,
         main="M0 M3 C3 J3")
#火山图top100
library('pheatmap')
select <- order(rowMeans(normalized_counts),decreasing=T)[1:200]
ntd <- normTransform(dds_norm) #default to log2FC+1
log2.norm.counts <- assay(ntd)[select,]
df <- as.data.frame(colData(dds_norm)[,c("name","condition")])
out<-pheatmap(log2.norm.counts, cluster_rows=T, show_rownames=F,
              cluster_cols=F, annotation_col=df,scale="row",
              treeheight_col=0,clustering_method = "average",
              color = colorRampPalette(c("navy", "white", "firebrick3"))(491),
              clustering_distance_rows = "maximum",cutree_cols = 2,
              main="M0 M3 C3 J3 top200")

#从这里开始为后面分析通路输出要用的表格
res <- results(dds_norm,contrast=c("condition","2T","2C"))
head(res)
#write.table(res,"result.txt", sep = ",", row.names = TRUE)
#res2 <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
#head(res2)
#sum(res2$padj < 0.05, na.rm=TRUE)
res3 <- cbind(ID=rownames(res), as.data.frame(res))
head(res3)
res4 <- res3[order(res3$log2FoldChange),]
head(res4)

#给ENSEMBL加SYMBOL
library('biomaRt')
library("curl")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
my_ensembl_gene_id<-row.names(res4 )
head(my_ensembl_gene_id)
hg_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name'),
                   filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(hg_symbols)
ensembl_gene_id<-rownames(res4 )
diff_gene_deseq2<-cbind(ensembl_gene_id,res4 )
colnames(diff_gene_deseq2)[1]<-c("ensembl_gene_id")
diff_name<-merge(diff_gene_deseq2,hg_symbols,by="ensembl_gene_id") #把注释文件和基因表达量文件合并起来
head(diff_name)
res4_order<- diff_name[order(diff_name$log2FoldChange),]
head(res4_order)
write.csv(res4_order,file="res4_log_order_2T_vs_2C.csv") #差异表达基因且有SYMBOL

library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(tidyverse)
library(enrichplot)

#给symbol注释GeneID
res4_order$external_gene_name <- as.character(res4_order$external_gene_name)
test2 = bitr(res4_order$external_gene_name, fromType="SYMBOL", toType=c("ENSEMBL","ENTREZID"), OrgDb="org.Hs.eg.db")
names(test2) = c("external_gene_name", "ENSEMBL","ENTREZID")
head(test2)
library(dplyr)
join<- left_join(res4_order, test2, by = "external_gene_name")

#输出一个Subdat和Input表，做后续的通路分析
subDat<- filter(join, padj <= 0.05 & external_gene_name != "NA")
write.csv(subDat,"subDat.csv",row.names =FALSE,quote = TRUE) #导出的Subdat表格用于做GSEA分析
head(subDat)
input <- subDat[,c(9,4,8,10,11)]
sampleNames <- c("SYMBOL","logFC","adj.P.Val","ENSEMBL", "ENTREZID")
names(input)[1:5] <- sampleNames
head(input)
write.table(input,"input.txt",sep="\t",row.names =F,col.name=F,quote = F)
write.csv(input,"input.csv",row.names =FALSE,quote = TRUE)
#一直要run到这里，写出一个input表格，用于做通路分析
