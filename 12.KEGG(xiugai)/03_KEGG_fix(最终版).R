# ============================================================
# Reviewer Comment 8 fix: remove leftover template artifacts.
# Original bugs in Palmitoylation12.KEGG.R:
#   (1) keggId = "hsa05010"  -> Alzheimer disease (template placeholder)
#   (2) reads "Exosome.diffGenes.txt" instead of "Palmitoylation.diffGenes.txt"
#   (3) gene = rt$entrezID  -> column is actually "entrezIDs"
# ============================================================
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")

setwd("E:\\F盘\\5.桌面文件\\3.私人资料\\论文第六次修改\\导入github的数据R语言（R4.4.3）\\12.KEGG(xiugai)")   # <<<< CHANGE

rt <- read.table("Palmitoylation.diffGenes.txt", header = TRUE, sep = "\t",
                 check.names = FALSE)                     # (2) fixed
genes <- unique(as.vector(rt[, 1]))
entrezIDs <- as.character(mget(genes, org.Hs.egSYMBOL2EG, ifnotfound = NA))
rt <- cbind(rt, entrezIDs)
colnames(rt)[1] <- "gene"
rt <- rt[rt[, "entrezIDs"] != "NA", , drop = FALSE]
gene <- rt$entrezIDs                                     # (3) fixed
logFC <- rt$logFC
names(logFC) <- gene

options(timeout = 99999)
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 1, qvalueCutoff = 1)
kk@result$Description <- gsub(" - Homo sapiens \\(human\\)", "", kk@result$Description)
KEGG <- as.data.frame(kk)
KEGG$geneID <- as.character(sapply(KEGG$geneID, function(x)
  paste(rt$gene[match(strsplit(x, "/")[[1]], as.character(rt$entrezIDs))], collapse = "/")))
KEGG <- KEGG[(KEGG$pvalue < 0.05 & KEGG$p.adjust < 0.05), ]
write.table(KEGG, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# --- single-pathway figure: pick the TOP significant pathway, not a hard-coded ID ---
top <- KEGG$ID[1]                                          # (1) fixed
if (!is.na(top)) {
  kk2 <- setReadable(kk, OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
  pdf(file = "pathview_top.pdf", width = 8, height = 8)
  print(cnetplot(kk2, showCategory = 10, foldChange = 2^logFC))
  dev.off()
  cat("Top KEGG pathway:", top, "\n")
}
