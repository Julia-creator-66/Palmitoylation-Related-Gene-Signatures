install.packages("colorspace")
install.packages("stringi")
install.packages("tidytable")
install.packages("ggplot2")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("DOSE")
BiocManager::install("clusterProfiler")
BiocManager::install("enrichplot")


#引用包
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("tidytable")
library("ggplot2")

pvalueFilter=0.05      #p值过滤条件
adjPvalFilter=0.05     #矫正后p值的过滤条件

#定义图形的颜色
colorSel="p.adjust"
if(adjPvalFilter>0.05){
	colorSel="pvalue"
}

drugFile="DSigDB_All_detailed.txt"      #药物和基因的关系文件
hubFile="interGenes.txt"                  #交集特征基因的文件
setwd("C:\\Users\\ljz\\Desktop\\27.drugEnrich")      #设置工作目录

#读取交集特征基因的列表文件
rt=read.table(hubFile, header=F, sep="\t", check.names=F)
#提取基因的名称
genes=unique(as.vector(rt[,1]))

#读取药物和基因的关系文件
drugRT=read.table(drugFile, header=T, sep="\t", check.names=F, quote="", comment.char="")
drugRT=drugRT[,1:2]

#药物富集分析
kk=enricher(genes,
	pvalueCutoff=1, qvalueCutoff=1,
	minGSSize = 10, maxGSSize = 500,
	TERM2GENE=drugRT)

#输出显著富集的结果
DRUG=as.data.frame(kk)
DRUG=DRUG[(DRUG$pvalue<pvalueFilter & DRUG$p.adjust<adjPvalFilter),]
write.table(DRUG[,-c(3,4)], file="DRUG.enrich.xls", sep="\t", quote=F, row.names = F)

#设置展示药物的数目
showNum=30
if(nrow(DRUG)<showNum){
	showNum=nrow(DRUG)
}

#柱状图
pdf(file="barplot.pdf", width=11, height=7)
barplot(kk, drop=TRUE, showCategory=showNum, label_format=100, color=colorSel)
dev.off()

#气泡图
pdf(file="bubble.pdf", width=11, height=7)
dotplot(kk, showCategory=showNum, orderBy="GeneRatio", label_format=100, color=colorSel)
dev.off()

#输出网络关系文件
DRUG=DRUG[,c("ID", "geneID")]
DRUG2=separate_rows(DRUG, geneID, sep="/")
DRUG2=DRUG2 %>% group_by(geneID) %>% slice_head(n=20)     #设置展示药物的数目
networkTab=cbind(DRUG2, "Drug")
colnames(networkTab)=c("Node1", "Node2", "Type")
write.table(file="net.network.txt", networkTab, sep="\t", quote=F, row.names=F)
#输出节点属性文件
drugNode=cbind(unique(networkTab[,1]),"Drug")
colnames(drugNode)=c("Node", "Type")
geneNode=cbind(unique(networkTab[,2]),"Gene")
colnames(geneNode)=c("Node", "Type")
nodeTab=rbind(drugNode, geneNode)
write.table(file="net.node.txt", nodeTab, sep="\t", quote=F, row.names=F)





