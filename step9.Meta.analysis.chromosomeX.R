 setwd("/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/")
female=read.table("Female.output.txt",header=T)
 male=read.table("Male.output.txt",header=T)

 female$Estimate=log(female$OR)
 female$SE=female$Estimate/female$STAT
 head(female)
 male$Estimate=log(male$OR)
 male$SE=male$Estimate/male$STAT
 female01=female[,c(2,10,11,9)]
 male01=male[,c(2,10,11,9)]
 names(female01)=c("SNP","beta.female","SE.female","P.female")
 names(male01)=c("SNP","beta.male","SE.male","P.male")
 inter=intersect(female01$SNP,male01$SNP)

all=merge(female01,male01,by="SNP")
all=subset(all, (SE.female !=0) & (SE.male !=0) )
se=cbind(all$SE.female,all$SE.male)
es=cbind(all$beta.female, all$beta.male)
library(metafor)
result=c()

for ( i in 1:nrow(es)) { file=data.frame(ID=c("female","male"),Estimate=es[i,],Std=se[i,]);meta=rma(yi=Estimate,data=file,sei=Std,method="FE");result=rbind(result,c(summary(meta)$beta,summary(meta)$se,summary(meta)$pval,summary(meta)$I2,summary(meta)$ QE,summary(meta)$ QEp))}


result=as.data.frame(result)
str(result)

colnames(result)=c("meta.beta","meta.se","meta.pval","meta.I2","meta.Q","meta.Qpval")
final=as.data.frame(cbind(all,result))
final=final[order(final$meta.pval),]
write.table(final,"Meta.Results.txt",row.names=F)


