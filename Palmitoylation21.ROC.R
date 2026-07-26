install.packages("glmnet")
install.packages("pROC")
install.packages("ggsci")


#引用包
library(glmnet)
library(pROC)
library(ggsci)

expFile="merge.normalize.txt"     #表达数据文件
geneFile="interGenes.txt"          #基因列表文件
setwd("C:\\Users\\ljz\\Desktop\\21.ROC")    #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#提取样品的分组信息(对照组和实验组)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", colnames(rt))
y=ifelse(y=="Control", 0, 1)

#读取基因列表文件
geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

#定义图形的颜色
bioCol=pal_aaas("default")(nrow(geneRT))

#对基因进行循环，绘制ROC曲线
aucText=c()
k=0
for(x in as.vector(geneRT[,1])){
	k=k+1
	#绘制ROC曲线
	roc1=roc(y, as.numeric(rt[x,]))     #得到ROC曲线的参数
	if(k==1){
		pdf(file="ROC.genes.pdf", width=5.5, height=5)
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}else{
		plot(roc1, print.auc=F, col=bioCol[k], legacy.axes=T, main="", add=TRUE, lwd=3)
		aucText=c(aucText, paste0(x,", AUC=",sprintf("%.3f",roc1$auc[1])))
	}
}
#绘制图例，得到ROC曲线下的面积
legend("bottomright", aucText, lwd=3, bty="n", cex=1, col=bioCol[1:(ncol(rt)-1)])
dev.off()

#构建逻辑模型的模型
rt=rt[as.vector(geneRT[,1]),]
rt=as.data.frame(t(rt))
logit=glm(y ~ ., family=binomial(link='logit'), data=rt)
pred=predict(logit, newx=rt)     #得到模型的打分

#绘制模型的ROC曲线
roc1=roc(y, as.numeric(pred))      #得到模型ROC曲线的参数
ci1=ci.auc(roc1, method="bootstrap")     #得到ROC曲线下面积的波动范围
ciVec=as.numeric(ci1)
pdf(file="ROC.model.pdf", width=5.5, height=5)
plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main="Model", lwd=3)
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
dev.off()




