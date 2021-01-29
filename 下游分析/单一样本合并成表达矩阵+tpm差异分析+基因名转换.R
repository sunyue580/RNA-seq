setwd('C:\\Users\\think\\Desktop\\data2\\GSE130036_RAW')
##以下以GSE130036为例,下载了每个样本的tsv.gz文件
##每个样本文件格式如下
# target_id           	length	     eff_length	   est_counts	  tpm
# ENST00000621489.1_1	   2538	       2052.65	      68.5741	    1.43133
# ENST00000623554.3_1	   1137	       897.13	        0	          0
# ENST00000623919.1_1	   1422        1301.28	      2	          0.0658495
# ENST00000623933.1_1 	 668	       428.191	      0	          0
# MSTRG.5.1	             14263	     12398.4	      298.334	    1.03094
# ENST00000624444.1_1	   2477	       1765.97	      0	          0
# ENST00000623374.1_1	   6430	       5163.09	      136.129	    1.12963
# ENST00000621924.4_1	   2216	       1587.11      	0	          0
# ENST00000619488.1_1	   1985	       1356.77	      70.5598	    2.22816
# ENST00000615804.1_1	   2168	       1525.54	      0	          0
# MSTRG.5.5	             14645	     13178.3	      226.42	    0.736119

a = list.files()
n = length(a)     

rt <- read.table(file = a[1],header=TRUE)
gene <- rt$target_id
rt <- rt[,"tpm",drop=F]
colnames(rt)[1] <- sapply(strsplit(a[1],"_"),"[",1)

for (i in 2:n){
  rt0 = read.table(file = a[i], header=TRUE)
  rt0 <- rt0[,"tpm",drop=F]
  colnames(rt0)[1] <- sapply(strsplit(a[i],"_"),"[",1)
  
  rt = cbind(rt,rt0)
}
rt <- cbind(gene=gene,rt)
rt <- rt[grep("ENST",rt$gene),]
x <- strsplit(rt$gene,"\\.")
y <- sapply(x,"[",1) #提取列表第1个元素
rt <- cbind(gene=y,rt)
rt <- rt[,-2]


rt<-as.matrix(rt)   # 变成矩阵
rownames(rt)=rt[,1]     # 第一列命名为行名
exp<-rt[,2:ncol(rt)]   # 去掉第一列，剩下的全部为表达量数据
dimnames<-list(rownames(exp),colnames(exp))    # 命名行和列
exp<-matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)    # 全部设置为数值类型，避免有文本值

# library("impute")
library(limma)
# mat=impute.knn(exp)
# rt=mat$data      # 以上三个，补全缺失数据
# dim(rt)  #32855    44
# rt=avereps(rt)# 相同基因取平均值
rt=avereps(exp)
dim(rt)  #16780    44
rt <- tibble::rownames_to_column(data.frame(rt),var = "gene")

library(biomaRt)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes <- getBM(attributes=c('ensembl_transcript_id','hgnc_symbol','chromosome_name'), mart = ensembl)
genes <- genes[genes$chromosome_name %in% c(seq(1,22,1),"X","Y"),]
dim(genes) #231944      3
genes <- genes[!duplicated(genes$ensembl_transcript_id),]
dim(genes) #231920      3
genes <- genes[,-3]
genes <- genes[genes$hgnc_symbol != "",]
dim(genes) #199314      2

rt2 <- dplyr::inner_join(genes,rt,by=c("ensembl_transcript_id"="gene"))
dim(rt2) #152015     39
write.table(rt2,file = "../GSE130036.txt",row.names=F,quote = F,sep = "\t")

rt3 <- rt2[,-2]
rt3 <- tibble::column_to_rownames(rt3,var = "ensembl_transcript_id")
rt3 <- rt3[,c(29:37,1:28)]

group_list=c(rep('Normal',9),rep('Tumor',28))
## 强制限定顺序
group_list <- factor(group_list,levels = c("Normal","Tumor"),ordered = F)
#表达矩阵数据校正
exprSet <- rt3
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
dim(deg)
deg_diff <- deg[abs(deg$logFC)>1 & deg$adj.P.Val < 0.05,]
dim(deg_diff) #455   6

deg <- tibble::rownames_to_column(deg,var = "gene")
deg_diff <- tibble::rownames_to_column(deg_diff,var = "gene")

deg <- dplyr::inner_join(genes,deg,by=c("ensembl_transcript_id"="gene"))
deg_diff <- dplyr::inner_join(genes,deg_diff,by=c("ensembl_transcript_id"="gene"))

write.csv(deg,"../diff_all.csv",row.names = F,quote = F)
write.csv(deg_diff,"../diff.csv",row.names = F,quote = F)

