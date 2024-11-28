####################提取填补后男性和女性的PLINK文件的每个位点的MAF; 注意这里的MAF只是最小等位基因，还需要根据Male.output.txt文件和Male.Imputed.frq文件来推断哪个是效应等位基因，以及效应等位基因的频率如何计算。EAF并不等于MAF

######男性
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/"
plink --bfile ${wd}Male.Imputed --freq --out ${wd}Male.Imputed

#####女性
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/"
plink --bfile ${wd}Female.Imputed --freq --out ${wd}Female.Imputed


setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/")
library(data.table)
frq=fread("Female.Imputed.frq")
output=fread("Female.output.txt")
frq.new=frq[,c(2,3,4,5)]
all=merge(output,frq.new,key="SNP")
for ( i in 1: nrow(all)) { if  ( (all[i,7] > 1)  ) {all$EAF[i] <- 1-all[i,11]} else { all$EAF[i] <- all[i,11]} }
write.table(all,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Female.EAF.txt",row.names=F)


setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/")
library(data.table)
frq=fread("Male.Imputed.frq")
output=fread("Male.output.txt")
frq.new=frq[,c(2,3,4,5)]
all=merge(output,frq.new,key="SNP")
for ( i in 1: nrow(all)) { if  ( (all[i,7] > 1)  ) {all$EAF[i] <- 1-all[i,11]} else { all$EAF[i] <- all[i,11]} }
write.table(all,"/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Male.EAF.txt",row.names=F)


