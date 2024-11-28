######################################该脚本是对lead SNP 附近1Mb的基因做共定位分析；

library(coloc)
library(data.table)
library(locuscomparer) #####作图用

################加载meta分析结果
meta=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results/Meta.Results.All.SNP.txt")
gwas=meta[,c(3,1,12,10,11)]
head(gwas)
names(gwas)=c("rs_id","SNP","pval_nominal","beta","varbeta")
str=do.call(rbind,strsplit(gwas$SNP,"_"))
gwas$position=as.numeric(str[,2])

##############加载Tissue的eQTL数据;注意：该Tissue文件不存在，用“Tissue"只是为了后续写代码的时候，批量替换为相应的组织。

eqtl=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/GTEX.eQTL.all.assoc/Tissue.chrX.eQTL.hg38.txt")

dbSNP=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/chrX.NPR.dbSNP.txt")

names(dbSNP)=c("chr","position","rs_id")
head(eqtl)
data=do.call(rbind,strsplit(eqtl$V2,"_"))
head(data)
eqtl$position=as.numeric(data[,2])
eqtl.rsid=merge(eqtl,dbSNP,by="position")
eqtl.rsid$tissue="Tissue"

eqtl.final=eqtl.rsid[,c(12,3,8,9,10,5,7,1,2,13)]
names(eqtl.final)=c("rs_id","variant_id","pval_nominal","beta","varbeta","N","MAF","position","gene","tissue")

##################################截取GWAS lead SNP 上下游1Mb的SNP位点的summary statistics: 分为两部分：1)Meta的lead SNP , 2)Male 的lead SNP

#########Meta lead SNP1: rs3138874 , position:108694045
gwas=gwas[order(gwas$position),]
pos <- 108694045
block=c((pos-1000000):(pos+1000000))
region=subset(gwas, position %in% block)

############################合并region Gwas和eQTL 数据
all=merge(eqtl.final, region, by="rs_id", all=FALSE,suffixes=c("_eqtl","_gwas"))
all=all[-which(duplicated(all[,c(3:9,12:15)])),]
if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}

gene=unique(all$gene)
frame=c()
for ( i in gene) { input=subset(all, gene == i);input$se.gwas=input$varbeta_gwas;input$varbeta_gwas=input$se.gwas*input$se.gwas;input$se.eqtl=input$varbeta_eqtl;input$varbeta_eqtl=input$se.eqtl*input$se.eqtl;
if (length(which(duplicated(input$SNP))) > 0) {input=input[-(which(duplicated(input$SNP))),]};
if (length(which(duplicated(input$position_eqtl))) > 0) {input=input[-(which(duplicated(input$position_eqtl))),]};
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas,snp=input$rs_id,beta=input$beta_gwas,varbeta=input$varbeta_gwas,type="cc",  N=2750),dataset2=list(pvalues=input$pval_nominal_eqtl, snp=input$rs_id, type="quant", N=input$N, MAF=input$MAF));
sum=as.data.frame(t(result$ summary));
sum$gene= i; sum$tissue=input$tissue[1]; sum$class="Meta"; sum$SNP="rs3138874"; frame=as.data.frame(rbind(frame,sum))
}

write.table(frame,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/Tissue.Meta.rs3138874.Coloc.Result.txt",row.names=F)

###########Meta lead SNP2:rs5941025 , position:91895380
gwas=gwas[order(gwas$position),]
pos=91895380
block=c((pos-1000000):(pos+1000000))
region=subset(gwas, position %in% block)


#############################Male lead SNP1: rs73250616 , position:100738288
Male=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/tmp/Male.XWAS.rsID.hg38.txt")
gwas=Male[,c(11,2,8,6,7,1)]
head(gwas)
names(gwas)=c("rs_id","SNP","pval_nominal","beta","varbeta","position")
gwas$varbeta=var(gwas$beta)

eqtl.rsid$tissue="Tissue"
eqtl.final=eqtl.rsid[,c(12,3,8,9,10,5,7,1,2,13)]
names(eqtl.final)=c("rs_id","variant_id","pval_nominal","beta","varbeta","N","MAF","position","gene","tissue")

gwas=gwas[order(gwas$position),]
pos=100738288
block=c((pos-1000000):(pos+1000000))
region=subset(gwas, position %in% block)

############################合并region Gwas和eQTL 数据
all=merge(eqtl.final, region, by="rs_id", all=FALSE,suffixes=c("_eqtl","_gwas"))
all=all[-which(duplicated(all[,c(3:9,12:15)])),]
if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}


gene=unique(all$gene)
frame=c()
for ( i in gene) { input=subset(all, gene == i);input$se.gwas=input$varbeta_gwas;input$varbeta_gwas=input$se.gwas*input$se.gwas;input$se.eqtl=input$varbeta_eqtl;input$varbeta_eqtl=input$se.eqtl*input$se.eqtl;
if (length(which(duplicated(input$SNP))) > 0) {input=input[-(which(duplicated(input$SNP))),]};
if (length(which(duplicated(input$position_eqtl))) > 0) {input=input[-(which(duplicated(input$position_eqtl))),]};
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas,snp=input$rs_id,beta=input$beta_gwas,varbeta=input$varbeta_gwas,type="cc",  N=1623),dataset2=list(pvalues=input$pval_nominal_eqtl, snp=input$rs_id, type="quant", N=input$N, MAF=input$MAF));
sum=as.data.frame(t(result$ summary));
sum$gene= i; sum$tissue=input$tissue[1]; sum$class="Male"; sum$SNP="rs73250616"; frame=as.data.frame(rbind(frame,sum))
}

write.table(frame,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/Tissue.Male.rs73250616.Coloc.Result.txt",row.names=F)

#############################Male lead SNP2: rs3138874, position:108694045
gwas=gwas[order(gwas$position),]
pos=108694045
block=c((pos-1000000):(pos+1000000))
region=subset(gwas, position %in% block)

############################合并region Gwas和eQTL 数据
all=merge(eqtl.final, region, by="rs_id", all=FALSE,suffixes=c("_eqtl","_gwas"))
all=all[-which(duplicated(all[,c(3:9,12:15)])),]
if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}

gene=unique(all$gene)
frame=c()
for ( i in gene) { input=subset(all, gene == i);input$se.gwas=input$varbeta_gwas;input$varbeta_gwas=input$se.gwas*input$se.gwas;input$se.eqtl=input$varbeta_eqtl;input$varbeta_eqtl=input$se.eqtl*input$se.eqtl;
if (length(which(duplicated(input$position_eqtl))) > 0) {input=input[-(which(duplicated(input$position_eqtl))),]};
if (length(which(duplicated(input$SNP))) > 0) {input=input[-(which(duplicated(input$SNP))),]};
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas,snp=input$rs_id,beta=input$beta_gwas,varbeta=input$varbeta_gwas,type="cc",  N=1623),dataset2=list(pvalues=input$pval_nominal_eqtl, snp=input$rs_id, type="quant", N=input$N, MAF=input$MAF));
sum=as.data.frame(t(result$ summary));
sum$gene= i; sum$tissue=input$tissue[1]; sum$class="Male"; sum$SNP="rs3138874"; frame=as.data.frame(rbind(frame,sum))
}

write.table(frame,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/Tissue.Male.rs3138874.Coloc.Result.txt",row.names=F)

#############################Male lead SNP3:rs142591759, position:153994161
gwas=gwas[order(gwas$position),]
pos=153994161
block=c((pos-1000000):(pos+1000000))
region=subset(gwas, position %in% block)

############################合并region Gwas和eQTL 数据
all=merge(eqtl.final, region, by="rs_id", all=FALSE,suffixes=c("_eqtl","_gwas"))
all=all[-which(duplicated(all[,c(3:9,12:15)])),]
if ( length(which(all$MAF == 0)) > 0 ) {all=all[-which(all$MAF == 0),]}

gene=unique(all$gene)
frame=c()
for ( i in gene) { input=subset(all, gene == i);input$se.gwas=input$varbeta_gwas;input$varbeta_gwas=input$se.gwas*input$se.gwas;input$se.eqtl=input$varbeta_eqtl;input$varbeta_eqtl=input$se.eqtl*input$se.eqtl;
if (length(which(duplicated(input$position_eqtl))) > 0) {input=input[-(which(duplicated(input$position_eqtl))),]};
if (length(which(duplicated(input$SNP))) > 0) {input=input[-(which(duplicated(input$SNP))),]};
result <- coloc.abf(dataset1=list(pvalues=input$pval_nominal_gwas,snp=input$rs_id,beta=input$beta_gwas,varbeta=input$varbeta_gwas,type="cc",  N=1623),dataset2=list(pvalues=input$pval_nominal_eqtl, snp=input$rs_id, type="quant", N=input$N, MAF=input$MAF));
sum=as.data.frame(t(result$ summary));
sum$gene= i; sum$tissue=input$tissue[1]; sum$class="Male"; sum$SNP="rs142591759"; frame=as.data.frame(rbind(frame,sum))
}

write.table(frame,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/Tissue.Male.rs142591759.Coloc.Result.txt",row.names=F)
