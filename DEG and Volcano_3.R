res = read.csv("res4_log_order_M3_vs_J3.csv")
res2 <-subset (res,pvalue >10e-300)
head(res)
head(res2)

#提取差异基因分析的结果
diff_gene_deseq2<-subset(res,padj<0.05&abs(log2FoldChange) > 1)
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file = "DEG_M3_vs_J3.csv") #输出为一个文件

#火山图
library(ggplot2)
library(EnhancedVolcano)
library(airway)
library(magrittr)
library(ggrepel)
EnhancedVolcano(res2,
                lab =res2$external_gene_name,
                selectLab = c('SHROOM3','CCL8','CXCL9','CXCL11','PTGS2',
                              'OLR1','TNIP3','CSF3','TNFAIP6','CXCL10',
                              'CXCL3','IL1B','CLSTN2','ACOD1','IL1A',
                              'TLR7','GBP1','ADORA2A-AS1','DNAAF1','LINC01055'),
                labCol = 'black',
                boxedLabels = F,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'grey30',
                x = 'log2FoldChange',
                y = 'pvalue',
                pCutoff = 10e-5,
                FCcutoff = 2.0,
                title = 'M3 vs J3',
                pointSize = 2.0,
                labSize = 4.0,
                colAlpha = 0.75,
                legendLabSize = 8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE)
#火山图run到这里就可以了


#备用
EnhancedVolcano(res2,
                lab = rownames(res2),
                x = 'log2FoldChange',
                y = 'pvalue',
                selectLab = c('ENSG00000102794') ,#IRG1
                drawConnectors = TRUE,
                widthConnectors = 0.4,
                colConnectors = 'black',
                xlim = c(-15, 15),
                title = 'M155 0h vs 6h volcano',
                pCutoff = 10e-2,
                FCcutoff = 2,
                col=c('black', 'blue', 'green', 'red1'),
                colAlpha = 1,
                legend=c('NS','Log2 FC','P value',
                         'P value & Log2 FC'),
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 5.0)