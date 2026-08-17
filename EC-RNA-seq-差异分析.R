library(limma)
library(DESeq2)
library("RColorBrewer")
library("ggrepel")
library(amap)
library(gplots)
library(ggplot2)
library(BiocParallel)
library(pheatmap)
library(dplyr)
#install.packages("tidyr")  # 重新安装
library(tidyr)  

data<-read.table("E:/数据分析/ZZ5-LHB-mix-20260613/RNA-seq/RNA-seq all_counts.txt")

YXD=data.frame(TN1=data[,6],TN2=data[,7],TN3=data[,8],
                 WT1=data[,9],WT2=data[,10],WT3=data[,11])
rownames(YXD)=data[,1]
group=c("TN","TN","TN","WT","WT","WT")
YXD_2 <- YXD[rowSums(YXD)>10,] 
conditions<-group
batch<-factor(c(1,1,1,1,1,1))
sample<-data.frame(conditions,batch)   
rownames(sample)<-colnames(YXD_2)
#产生DESeq数据集并计算标准化因子

ddsFullCountTable <- DESeqDataSetFromMatrix(countData =YXD_2,colData = sample,  design= ~ conditions)
dds <- DESeq(ddsFullCountTable)

normalized_counts <- counts(dds, normalized=TRUE)
#normalized_counts_mad <- apply(normalized_counts, 1, mad) 
#normalized_counts <- normalized_counts[order(normalized_counts_mad, decreasing=T), ]
#write.table(normalized_counts, file="E:/数据分析/ZZ5-LHB-mix-20260613/Data.normalized.xls",quote=F,sep="\t", row.names=T, col.names=T)
'''
## log转换后的结果并输出
rld <- rlog(dds, blind=FALSE)
rlogMat <- assay(rld)
rlogMat <- rlogMat[order(normalized_counts_mad, decreasing=T), ]
write.table(rlogMat, file="E:/DESeq2.normalized.rlog.xls",quote=F, sep="\t", row.names=T, col.names=T)
'''
######样本相关性热图绘制及PCA分析
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pearson_cor <- as.matrix(cor(normalized_counts, method="pearson"))
hc <- hcluster(t(normalized_counts), method="pearson")
heatmap.2(pearson_cor, Rowv=as.dendrogram(hc), symm=T, trace="none",
          col=hmcol, margins=c(11,11), main="Pearson correlation")
#pca_data <- plotPCA(dds, intgroup=c("conditions"), returnData=T, ntop=3000)

sampleA = "TN"
sampleB = "WT"
contrastV <- c("conditions", sampleA, sampleB)
res <- results(dds,contrast=contrastV)
#res$padj[is.na(res$padj)] <- 1  #校正后p-value为NA的赋值为1
res <- cbind(ID=rownames(res),normalized_counts ,as.data.frame(res))
res <- res[order(res$padj),]
res_de <- subset(res, res$padj<0.05, select=c('ID','WT1','WT2','WT3',"TN1","TN2","TN3",'log2FoldChange', 'padj'))
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1))



write.table(res,"E:/数据分析/ZZ5-LHB-mix-20260613/ALL-data-nor-log.xls",sep="\t", quote=F, row.names=T,col.names = T)
write.table(res_de_up,"E:/数据分析/ZZ5-LHB-mix-20260613/up.xls",sep="\t", quote=F, row.names=T)
write.table(res_de_dw,"E:/数据分析/ZZ5-LHB-mix-20260613/down.xls",sep="\t", quote=F, row.names=T)

cut_off_qvalue = 0.05
cut_off_logFC = 1
res$Sig <- ifelse(res$padj < cut_off_qvalue & 
                    abs(res$log2FoldChange) > cut_off_logFC,
                  ifelse(res$log2FoldChange > cut_off_logFC ,'Enriched','Depleted'),'Similar')
res <- data.frame(res)
tmp <- res %>% drop_na(Sig)
table(res$Sig)
###add name delete #
#gene_tmp <- c("nuoJ","mukF","nlpE","glpK","ssrS","rrsD","rrsG","rrsH","rrsE","rrsB","ssrA","rrsA","rrsC")
#gene_tmp <- data.frame(gene_tmp)
#gene_tmp$geneList <- gene_tmp$gene_tmp
#ID <- res$ID
#tmp <- res %>% left_join(gene_tmp,by = c("ID" = "gene_tmp"))

ggplot(res, aes(x = res$log2FoldChange, y = -log10(res$padj), colour=Sig)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue","red","black")) + 
  xlim(c(-10, 10)) + 
  ylim(c(0,150))+ 
  geom_vline(xintercept=c(-cut_off_logFC,cut_off_logFC),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_qvalue),
             lty=4,col="black",lwd=0.8) +
  labs(x="Log2(FC)",
       y="-log10(Padj)")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank() 
  )# +  
#  geom_text_repel(aes(label=geneList), fontface="bold",
#                  color="black", box.padding=unit(0.35, "lines"),
#                  point.padding=unit(5, "lines"),segment.color ="black",max.overlaps=Inf)
#dev.off()

#logCounts <- log2(res$baseMean+1)
#logFC <- res$log2FoldChange
#FDR <- res$padj
#plot(logFC, -1*log10(FDR), col=ifelse(FDR<=0.01, "red", "black"),
#     xlab="logFC", ylab="-log10 Padj", main="Differentially expressed gene",
#     pch=1)

###choose 1
res_de_up_sorted <- res_de_up %>% arrange(padj)
res_de_dw_sorted <- res_de_dw %>% arrange(padj)
res_de_up_top20_id <- as.vector(head(res_de_up_sorted$ID,20))
res_de_dw_top20_id <- as.vector(head(res_de_dw_sorted$ID,20))
#choose 2
#res_de_up_top20_id <- as.vector(head(res_de_up$ID,20))
#res_de_dw_top20_id <- as.vector(head(res_de_dw$ID,20))

res_de_top20 <- c(res_de_up_top20_id, res_de_dw_top20_id)   
res_de_top20_expr <- normalized_counts[rownames(normalized_counts) %in% res_de_top20,]
res_de_top20_expr <-YXD_2[rownames(YXD_2) %in% res_de_top20,]
pheatmap(res_de_top20_expr, cluster_row=T, scale="row", annotation_col=sample)
library(pheatmap)


library(clusterProfiler)
gene2ko <- download_KEGG(organism = "eco", keggType = "KEGG")$KEGG2Gene
colnames(gene2ko) <- c("KO", "gene_id")
# 2. 获取KO-eco通路映射
ko2path <- download_KEGG(organism = "eco", keggType = "KEGG")$Pathway2Gene
colnames(ko2path) <- c("pathway_id", "KO")
# 3. 三表合并：gene-KO-eco通路
gene_ko_path <- merge(gene2ko, ko2path, by = "KO")
# 4. 导出csv本地保存
write.csv(gene_ko_path, "Ecoli_gene_KO_ecoPath.csv", row.names = F)








up_genes <- rownames(res_de_up)


Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/up-go.txt")
#symbols <-bitr(Tn5, fromType = "ALIAS", toType = "ENTREZID",OrgDb = 'org.EcK12.eg.db')

Tn5<-Tn2[,1]

up_genes_entrez <- bitr(up_genes, fromType = "SYMBOL",toType = "ENTREZID", OrgDb = org.EcK12.eg.db)
ego <- enrichGO(
  gene          = Tn5,
  keyType = "ENTREZID",
  OrgDb         = org.EcK12.eg.db,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.5,
  readable      = TRUE)
dotplot(ego)#,showCategory = 50)
#barplot(go_up, showCategory =8, split = "ONTOLOGY") +
#  facet_grid(ONTOLOGY ~ ., scales ="free")
write.csv(ego,"E:/数据分析/ZZ7-LHB-mix-20260701/KEGG GO/up-GO-gc.csv")#GO-result-ALL.csv")

Tn2<-read.table("E:/数据分析/ZZ7-LHB-mix-20260701/down-kegg.txt")
#symbols <-bitr(Tn5, fromType = "ALIAS", toType = "ENTREZID",OrgDb = 'org.EcK12.eg.db')

Tn5<-Tn2[,1]
kegg <- enrichKEGG(
  gene = Tn5,# up_genes
  organism = 'eco', 
  pAdjustMethod = 'fdr', 
  pvalueCutoff = 0.05, #0.2  #0.4  #0.3  #0.6
  qvalueCutoff = 0.5, 
)
dotplot(kegg)
write.csv(kegg,"E:/数据分析/ZZ7-LHB-mix-20260701/KEGG GO/dw-KEGG-GC.csv")

