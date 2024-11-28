library(data.table)
#########加载性别合并的PLINK文件
combine=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/Combined.Sex.Imputed.frq")

#######加载meta分析的结果
meta=fread("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results/Meta.Results.All.SNP.txt")


############由于PLINK计算出来的Freq是最小等位基因频率MAF，而EAF是效应等位基因频率,根据PLINK logistic结果文件来看，A1就是效应等位基因，A1对应的频率和MAF相等

combine=combine[,2:5]
all=merge(meta,combine,by="SNP")
all$EAF=all$MAF

write.table(all,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results/Meta.Results.All.SNP.txt",row.names=F)
