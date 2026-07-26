#install.packages("ggplot2")
#install.packages("GseaVis")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler")
#BiocManager::install("fgsea")


#引用包
library(ggplot2)
library(fgsea)
library(clusterProfiler)
library(GseaVis)

pvalueFilter=0.05      #p值的过滤条件
adjPvalFilter=0.05     #矫正后p值的过滤条件
degfile="all.txt"      #所有基因差异分析的结果
setwd("C:\\Users\\Administrator\\Desktop\\Exosome\\13.GSEA")     #设置工作目录

#读取输入文件
data=read.table(degfile, header=T, sep="\t", check.names=F)
#提取基因名称和logFC
genelist=data$id
log2fc=data$logFC

#将基因名字转换成ID
ID <- clusterProfiler::bitr(genelist, fromType="SYMBOL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
mergedata=merge(ID,data.frame(SYMBOL=genelist,fc=log2fc),by="SYMBOL",all=F)
#根据logFC对基因进行排序
mergedata=mergedata[order(mergedata$fc, decreasing=T),]
geneinput=mergedata$fc
names(geneinput)=mergedata$ENTREZID

#GSEA富集分析
gseo=gseKEGG(geneinput, pvalueCutoff=1, by="fgsea")
out=as.data.frame(gseo@result)
out$core_enrichment = as.character(sapply(out$core_enrichment,function(x)paste(ID$SYMBOL[match(strsplit(x,"/")[[1]],as.character(ID$ENTREZID))],collapse="/")))
out=out[(out$pvalue<pvalueFilter & out$p.adjust<adjPvalFilter),]
write.table(out,"GSEA.result.xls",sep="\t",col.names=T,row.names=F,quote=F)

#对富集的结果进行可视化
for(i in 1:nrow(out)){
   p1 = gseaNb(gseo,geneSetID = out$ID[i],legend.position=right, lineSize=1.25,
            subPlot=3, curveCol=c("darkgreen","darkgreen"), addPval=T, pvalX=0.8, pvalY=0.8)
   ggsave(sprintf("GSEA-%s.pdf",out$ID[i]),p1,height=7,width=6.5)
}

#对富集最显著的前五个通路进行可视化
pdf(file="top5.pdf", width=9, height=6.5)
gseaNb(object = gseo,
       geneSetID = out$ID[1:5],
       subPlot=3,
       newGsea = F,
       curveCol=c("#76BA99", "#EB4747", "#996699", "#5C88DA", "#FFCD00"),
       rmHt = F,
       addPval = T,
       pvalX = 0.8,
       pvalY = 0.8)
dev.off()





