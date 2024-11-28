setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/09.Functional.data/02.GTEX.eQTL/")
library(data.table)

###################GTEX的chcrX的eQTL源文件
eqtl=fread("chrX.eQTL.txt")
head(eqtl)
data=do.call(rbind,strsplit(eqtl$variant_id,"_"))
eqtl$SNP.BP=as.numeric(data[,2])

###########cross sex 的XWAS结果:hg38
SNP=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results/ESRD.SNP.10E-5.hg38.txt")

head(SNP)
vector=unlist(strsplit(SNP$SNP, "_"))
data <- as.data.frame(matrix(vector, ncol = 4, byrow = TRUE))
names(data)=c("chromosome","position","Alt","Ref")
all=as.data.frame(cbind(SNP,data))
head(all)
head(eqtl)


###########cross sex 的XWAS结果:hg19
SNP=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results/Meta.XWAS.10E-5.SNP.hg38.txt")

all=SNP



#############注释，过滤cross sex独立的两位点
pos=c(91895380,108694045)
sub=subset(eqtl,SNP.BP %in% all$position)
head(sub)
sub$position=as.numeric(sub$SNP.BP)

#######合并注释结果和原始的summary 结果

result=merge(sub,all,by="position")
result$SNP.BP.y=c()
result$SNP.BP.x=c()
result$Alt=c()
result$Ref=c()

###############对于重复的注释行，选择样本量最大的行 (这是嵌套循环：首先选择position重复的行，然后根gene_ID选择子集,然后再在该子集里面选择样本量最大的行)########该步骤可以忽略
###pos=unique(result$position)
###frame=c()

#for ( i in 1:length(pos)) { file=subset(result, position %in% pos[i] ); final=c();for ( j in 1:length(unique(file$gene_id))) { data=file[which(file$gene_id==unique(file$gene_id)[j]),];data01=data[which(data$ma_samples==max(data$ma_samples)),]; final= as.data.frame(rbind(final, data01))}; frame= as.data.frame(rbind(frame, final))}
#df=unique(frame)

write.table(result,"ESRD.XWAS.eQTL.Gene.hg19.txt",row.names=F)






write.table(result,"ESRD.XWAS.eQTL.Gene.hg19.txt",row.names=F)

###################将基因版本转换为 gene symbol
library(org.Hs.eg.db)
library('clusterProfiler')

V2=sub("\\..*", "", result$V2)
gene<- bitr(V2, fromType = "ENSEMBL", toType=c("SYMBOL"),OrgDb = org.Hs.eg.db)
head(gene)
 names(gene)[1]="geneId"
###两个cross sex独立的SNP位点为chrX_91895380_G_T 和 chrX_108694045_A_G: 无注释，能被注释到的都是非独立位点

