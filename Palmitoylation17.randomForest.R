#install.packages("randomForest")
#install.packages("ggplot2")


#引用包
set.seed(12345)
library(randomForest)
library(ggplot2)

expFile="logistic.sigExp.txt"        #逻辑回归显著基因的表达文件
setwd("C:\\Users\\Administrator\\Desktop\\Exosome\\17.RF")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
data=rt[,-ncol(rt)]

#获取样品的分组信息(对照组和实验组)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
data=cbind(group, data)
data$group=factor(data$group, levels=c("Control","Treat"))

#随机森林树
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees

#获取基因重要性的评分
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)
importance=importance(x=rf2)

#筛选疾病的特征基因
rfGenes=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),]
rfGenes=names(rfGenes[rfGenes>4])        #挑选重要性评分大于4的基因
#rfGenes=names(rfGenes[1:10])           #挑选重要性评分最高的5个基因
write.table(rfGenes, file="RF.gene.txt", sep="\t", quote=F, col.names=F, row.names=F)

#绘制基因重要性的气泡图
importance=importance[order(importance[,"MeanDecreaseGini"], decreasing = TRUE),,drop=F]
rt=cbind(row.names(importance), as.data.frame(importance))
colnames(rt)=c("ID", "Importance")
rt=rt[1:30,]
rt$ID=factor(rt$ID, levels=rev(rt$ID))
rt$Importance=as.numeric(rt$Importance)
p=ggplot(rt, aes(Importance, ID)) + 
    geom_segment(aes(x=0, xend=Importance, y=ID, yend=ID), color="grey", cex=1.5) + 
    geom_point(aes(color=Importance), cex=3.5)+
    scale_colour_gradient(low="#5050FFFF", high="#BB0021FF") + 
    labs(x="Importance", y="")+
    theme(axis.text.x=element_text(color="black", size=10),
          axis.text.y=element_text(color="black", size=10))+
    theme_bw()

#输出图形
pdf(file="geneImportance.pdf", width=5, height=4.5)
print(p)
dev.off()




