install.packages("glmnet")


#引用包
library(glmnet)
expFile="merge.normalize.txt"         #表达数据文件
geneFile="Exosome.diffGenes.txt"     #基因列表文件
setwd("C:\\Users\\ljz\\Desktop\\14")    #设置工作目录

#读取输入文件
data=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)

#读取基因列表文件,提取交集基因的表达量
geneRT=read.table(geneFile, header=T, sep="\t", check.names=F)
data=data[as.vector(geneRT[,1]),]

#获取样品分组信息(对照组和实验组)
data=t(data)
group=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(data))
rt=as.data.frame(data)
rt$Type=ifelse(group=="Control", 0, 1)

#单因素逻辑回归分析
outTab=data.frame()
for(i in colnames(rt)[1:(ncol(rt)-1)]){
	rt2=rt[,c(i, "Type")]
	colnames(rt2)=c("Gene", "Type")
	fit=glm(Type ~ Gene, family="binomial", data=rt2)
	summ2=summary(fit)
	conf=confint(fit, level=0.95)
	geneTab=cbind(Gene=i,
	               OR=exp(summ2$coefficients[,"Estimate"])[2],
	               OR.95L=exp(conf[,1])[2],
	               OR.95H=exp(conf[,2])[2],
	               pvalue=summ2$coefficients[,"Pr(>|z|)"][2])
	outTab=rbind(outTab, geneTab)
}

#输出逻辑回归分析结果文件
outTab=outTab[as.numeric(outTab$pvalue)<0.05,]
outTab=outTab[order(as.numeric(outTab$pvalue)),]
write.table(outTab, file="logistic.result.txt", sep="\t", quote=F, row.names=F)
#输出逻辑回归显著基因的表达数据
sigExp=rt[,c(as.vector(outTab[,"Gene"]),"Type")]
sigExpOut=cbind(id=row.names(sigExp), sigExp)
write.table(sigExpOut, file="logistic.sigExp.txt", sep="\t", quote=F, row.names=F)





