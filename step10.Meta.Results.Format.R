##########该脚本是将最后meta的结果进行整理：1）去除I2>20的位点 
 setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/")
meta=read.table("Meta.Results.txt",header=T)

result=subset(meta,meta.I2 <= 20 & meta.Qpval > 0.1/nrow(meta))
 p_value=as.numeric(result$meta.pval)
median(qchisq(p_value, df=1, lower.tail=T)) / qchisq(0.5, 1)

p.adj=p.adjusted(result$meta.pval,method="fdr")
result$meta.pval.fdr=p.adj


setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.XWAS.Meta.Results")

write.table(result,"Meta.Results.txt",row.names=F)
