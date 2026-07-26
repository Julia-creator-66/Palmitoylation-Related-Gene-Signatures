#install.packages("colorspace")
#install.packages("stringi")
#install.packages("ggplot2")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("DOSE")
#BiocManager::install("clusterProfiler")
#BiocManager::install("enrichplot")
#BiocManager::install("pathview")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")

pvalueFilter=0.05      #p值过滤条件
adjPvalFilter=0.05     #矫正后的p值过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

setwd("C:\\Users\\Administrator\\Desktop\\Exosome\\12.KEGG")      #设置工作目录
rt=read.table("Exosome.diffGenes.txt", header=T, sep="\t", check.names=F)     #读取基因列表文件

#提取基因的名称, 将基因名称转换为基因id
genes=unique(as.vector(rt[,1]))
entrezIDs=mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
entrezIDs=as.character(entrezIDs)
rt=cbind(rt, entrezIDs)
colnames(rt)[1]="gene"
rt=rt[rt[,"entrezIDs"]!="NA",]      #去除基因id为NA的基因
gene=rt$entrezID
#gene=gsub("c\\(\"(\\d+)\".*", "\\1", gene)
logFC=rt$logFC
names(logFC)=gene

#通路富集分析
kk <- enrichKEGG(gene=gene, organism="hsa", pvalueCutoff=1, qvalueCutoff=1)
kk@result$Description=gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$gene[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$p.adjust<adjPvalFilter),]
#输出显著富集的结果
write.table(KEGG, file="KEGG.txt", sep="\t", quote=F, row.names = F)

#设置展示通路的数目
showNum=30
if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

#柱状图
pdf(file="barplot.pdf", width=8, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=8, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#绘制基因和通路的网络关系图
pdf(file="cnetplot.pdf", width=8, height=6.5)
kk2=setReadable(kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
cnet=cnetplot(kk2, layout=igraph::layout_in_circle, showCategory=10, foldChange=2^logFC)
print(cnet)
dev.off()

#绘制通路图
keggId="hsa05010"      #输入通路ID，进行通路图的绘制(需要修改)
pv.out=pathview(gene.data=logFC, pathway.id=keggId, species="hsa", out.suffix="pathview")




