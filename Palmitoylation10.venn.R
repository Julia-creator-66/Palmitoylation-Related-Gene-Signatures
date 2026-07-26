install.packages("ggvenn")


#引用包
library(ggvenn)

setwd("C:\\Users\\Administrator\\Desktop\\肺癌棕榈酰化（外泌体）\\10.venn")    #设置工作目录
geneList=list()

#读取差异分析的结果文件
diffRT=read.table("diff.txt", header=T, sep="\t", check.names=F)
geneNames=as.vector(diffRT[,1])     #提取差异基因的名称
geneList[["DEG"]]=geneNames

#读取外泌体的基因列表文件
rt=read.table("gene.txt", header=F, sep="\t", check.names=F)
geneNames=as.vector(rt[,1])       #提取外泌体相关的基因
geneList[["Exosome"]]=geneNames

#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList, show_percentage = T,
	stroke_color = "white", stroke_size = 0.5,
	fill_color = c("#E41A1C","#1E90FF"),
	set_name_color =c("#E41A1C","#1E90FF"),
	set_name_size=6, text_size=4.5)
dev.off()

#输出交集基因
interGenes=Reduce(intersect, geneList)
outTab=diffRT[diffRT[,1] %in% interGenes,]
write.table(file="Exosome.diffGenes.txt", outTab, sep="\t", quote=F, row.names=F)





