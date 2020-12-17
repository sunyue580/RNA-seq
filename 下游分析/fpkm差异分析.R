setwd('C:\\Users\\think\\Desktop\\data3')

rt <- read.csv("gene_expression_TCGA.csv",header = T,row.names = NULL)
rt <- rt[,-c(1,3)]

library(tibble)
rt <- column_to_rownames(rt,var = "gene_name")

expMatrix <- rt
fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm+1) - log(sum(fpkm+1)) + log(1e6))
}
tpms <- apply(expMatrix,2,fpkmToTpm)
tpms[1:3,]
colSums(tpms)

x <- substring(colnames(tpms),14,15)

tpms_con <- tpms[,x=="11"]
tpms_treat <- tpms[,x=="01"]
tpms2 <- cbind(tpms_con,tpms_treat)

group_list=c(rep('Normal',32),rep('Tumor',375))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Normal","Tumor"),ordered = F)
#表达矩阵数据校正
exprSet <- tpms2
pdf("x.pdf")
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
dev.off()
library(limma) 
exprSet=normalizeBetweenArrays(exprSet)
boxplot(exprSet,outline=FALSE, notch=T,col=group_list, las=2)
#判断数据是否需要转换
exprSet <- log2(exprSet+1)
#差异分析：
dat <- exprSet
design=model.matrix(~factor( group_list ))
fit=lmFit(dat,design)
fit=eBayes(fit)
options(digits = 4)
topTable(fit,coef=2,adjust='BH')
bp=function(g){
  library(ggpubr)
  df=data.frame(gene=g,stage=group_list)
  p <- ggboxplot(df, x = "stage", y = "gene",
                 color = "stage", palette = "jco",
                 add = "jitter")
  #  Add p-value
  p + stat_compare_means()
}
deg=topTable(fit,coef=2,adjust='BH',number = Inf)
head(deg) 
write.csv(deg,"diff_all.csv")
