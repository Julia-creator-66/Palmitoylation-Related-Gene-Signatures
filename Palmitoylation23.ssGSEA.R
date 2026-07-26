if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSVA")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GSEABase")


#引用包
library(GSVA)
library(GSEABase)

expFile <- "merge.normalize.txt"
gmtFile <- "immune.gmt"
setwd("C:\\Users\\ljz\\Desktop\\23.ssGSEA")

# 读取并处理表达矩阵
rt <- read.table(expFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
mat <- as.matrix(rt)
mat <- matrix(as.numeric(mat), 
             nrow = nrow(mat),
             dimnames = list(rownames(mat), colnames(mat)))
mat <- mat[!apply(mat, 1, function(x) any(is.na(x))), ]  # 移除NA值
mat <- avereps(mat)
mat <- mat[rowMeans(mat) > 0, ]

# 读取基因集
geneSet <- getGmt(gmtFile, geneIdType=SymbolIdentifier())

# 关键修正：使用新版参数调用方式
params <- ssgseaParam(exprData = mat,
                     geneSets = geneSet,
                     alpha = 0.25,      # ssGSEA特有参数，控制权重分布
                     normalize = TRUE)  # 是否进行内置标准化

# 执行分析
ssgseaScore <- gsva(params)

# 自定义0-1标准化（按列）
normalize <- function(x) { (x - min(x)) / (max(x) - min(x)) }
normalScore <- apply(ssgseaScore, 2, normalize)

# 输出结果
ssgseaOut <- rbind(colnames(normalScore), normalScore)
write.table(ssgseaOut, "immuneScore.txt", sep="\t", quote=FALSE, col.names=FALSE)



