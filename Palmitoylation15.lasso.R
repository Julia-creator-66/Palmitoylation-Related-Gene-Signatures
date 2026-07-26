#install.packages("glmnet")


#引用包
library(glmnet)
#set.seed(12345)

expFile="logistic.sigExp.txt"        #逻辑回归显著基因的表达文件
setwd("C:\\Users\\ljz\\Desktop\\15.lasso")      #设置工作目录

#读取表达数据文件
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
rt=rt[,-ncol(rt)]

#构建lasso回归模型
x=as.matrix(rt)
y=gsub("(.*)\\_(.*)\\_(.*)", "\\3", row.names(rt))
fit=glmnet(x, y, family = "binomial", alpha=1)
#绘制Lasso回归的图形
pdf(file="lasso.pdf", width=6, height=5.5)
plot(fit)
dev.off()
#绘制交叉验证的图形
cvfit=cv.glmnet(x, y, family="binomial", alpha=1, type.measure='deviance', nfolds=10)
pdf(file="cvfit.pdf", width=6, height=5.5)
plot(cvfit)
dev.off()

#输出lasso回归找到疾病的特征基因
coef=coef(fit, s=cvfit$lambda.min)
index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)





