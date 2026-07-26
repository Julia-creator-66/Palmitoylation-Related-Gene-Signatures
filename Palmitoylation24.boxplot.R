install.packages("reshape2")
install.packages("ggpubr")
install.packages("corrplot")


#引用包
library(reshape2)
library(ggpubr)
library(corrplot)

inputFile="immuneScore.txt"     #免疫细胞的打分文件
setwd("C:\\Users\\ljz\\Desktop\\24.boxplot")     #设置工作目录

#读取免疫细胞的打分文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)
rt=t(rt)

#对样品进行分组(对照组和实验组)
con=grepl("_Control", rownames(rt), ignore.case=T)
treat=grepl("_Treat", rownames(rt), ignore.case=T)
conData=rt[con,]
treatData=rt[treat,]
conNum=nrow(conData)
treatNum=nrow(treatData)
data=t(rbind(conData, treatData))

##################绘制箱线图##################
#把数据转换成ggplot2输入文件
Type=gsub("(.*)\\_(.*)", "\\2", colnames(data))
data=cbind(as.data.frame(t(data)), Type)
data=melt(data, id.vars=c("Type"))
colnames(data)=c("Type", "Immune", "Expression")
#绘制箱线图
group=levels(factor(data$Type))
bioCol=c("#008B45FF","#EE0000FF","#0066FF","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(group)]
boxplot=ggboxplot(data, x="Immune", y="Expression", fill="Type",
				  xlab="",
				  ylab="Fraction",
				  legend.title="Type",
				  #notch=T, add="point",
				  width=0.8,
				  palette=bioCol)+
				  rotate_x_text(60)+
	stat_compare_means(aes(group=Type),symnum.args=list(cutpoints=c(0, 0.001, 0.01, 0.05, 1), symbols=c("***", "**", "*", "")), label="p.signif")
#输出箱线图
pdf(file="immune.diff.pdf", width=8, height=6)
print(boxplot)
dev.off()




