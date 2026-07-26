#install.packages("ggvenn")


#引用包
library(ggvenn)

setwd("C:\\Users\\ljz\\Desktop\\18.venn")    #设置工作目录
files=list.files(pattern="*.txt$")     #获取目录下txt结尾的文件
geneList=list()

#读取所有txt文件中的基因信息，保存到geneList
for(inputFile in files){
	if(inputFile=="interGenes.txt"){next}
    rt=read.table(inputFile, header=F, sep="\t", check.names=F)      #读取输入文件
    geneNames=unlist(strsplit(as.vector(rt[,1]), " "))       #提取基因名称
    geneNames=gsub("^ | $","",geneNames)     #去掉基因首尾的空格
    uniqGene=unique(geneNames)                 #基因取unique
    header=unlist(strsplit(inputFile,"\\.|\\-"))
    geneList[[header[1]]]=uniqGene
}
#绘制venn图
pdf(file="venn.pdf", width=6, height=6)
ggvenn(geneList, show_percentage = T,
  stroke_color = "white", stroke_size = 0.5,
  fill_color = c("#E41A1C","#1E90FF","#FF8C00"),
  set_name_color =c("#E41A1C","#1E90FF","#FF8C00"),
  set_name_size=6, text_size=4.5)
dev.off()

#输出三种方法的交集特征基因
interGenes=Reduce(intersect, geneList)
write.table(file="interGenes.txt", interGenes, sep="\t", quote=F, col.names=F, row.names=F)




