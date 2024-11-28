setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/09.Functional.data/02.GTEX.eQTL/")
library(data.table)

###################GTEX的chcrX的eQTL源文件:坐标hg38
eqtl=fread("chrX.eQTL.txt")
head(eqtl)
data=do.call(rbind,strsplit(eqtl$variant_id,"_"))
eqtl$SNP.BP=as.numeric(data[,2])


###########Female 的XWAS结果:hg38
SNP=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.2.Female.XWAS.analysis/Female.XWAS.10E-5.hg38.txt")

all=SNP
head(all)
head(eqtl)
all$position=all$BP


#############注释

sub=subset(eqtl,SNP.BP %in% all$position)
head(sub)
sub$position=as.numeric(sub$SNP.BP)

#######合并注释结果和原始的summary 结果
result=merge(sub,all,by="position")

###############对于重复的注释行，选择样本量最大的行 (这是嵌套循环：首先选择position重复的行，然后根gene_ID选择子集,然后再在该子集里面选择样本量最大的行)########该步骤可忽略
#pos=unique(result$position)
#frame=c()

#for ( i in 1:length(pos)) { file=subset(result, position %in% pos[i] ); final=c();for ( j in 1:length(unique(file$gene_id))) { data=file[which(file$gene_id==unique(file$gene_id)[j]),];data01=data[which(data$ma_samples==max(data$ma_samples)),]; final= as.data.frame(rbind(final, data01))}; frame= as.data.frame(rbind(frame, final))}
#df=unique(frame)

write.table(result,"Female.XWAS.eQTL.Gene.hg19.txt",row.names=F)

###################将基因版本转换为 gene symbol
library(org.Hs.eg.db)
library('clusterProfiler')


gene<- bitr(result$V2, fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)
head(gene)
 names(gene)[1]="geneId"
