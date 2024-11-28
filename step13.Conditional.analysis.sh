##################该脚本是conditional & joint association analysis using GWAS summary statistics；具体细节参考网页:https://yanglab.westlake.edu.cn/software/gcta/#Overview

#############Step1.input文件整理为特定格式/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/All.SNP.txt：##有beta及P值，计算SE：A$SE <- sqrt(((A$beta)^2)/qchisq(A$P, 1, lower.tail = F)); names(data)=c("SNP", "A1", "A2", "freq", "b", "se", "p"); 这里的freq是EAF，它的值等于MAF；这是因为PLINK的logistic回归分析结果文件就是将最小等位基因作为效应基因的（可通过logistic结果文件：*output.txt和PLINK输出的*.freq文件，以及bim文件进行比对，便可确定，三个文件对应列是一致的，。注意：bim文件倒数第二列总是最小等位基因，倒数第二列是主要等位基因）

#SNP A1 A2 freq b se p N
#chrX_2781514_C_A C A 0.57 0.02 0.043 0.7
#chrX_2781927_A_G G A 0.64 0.03 0.046 0.43


##################Step2.将男性和女性的基因型整合到一个大文件PLINK当中： /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/Combined.Sex.Imputed,并更新性别信息:/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/Updata.Sex.txt; 男性为1 女性为2！！！！！



################step3.1 对cross sex的结果进行conditional & joint association analysis

wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/"
gcta --bfile ${wd}Combined.Sex.Imputed --maf 0.01 --cojo-file ${wd}Combine.COJO.input.files.txt --cojo-p 0.0001 --update-sex ${wd}Updata.Sex.txt --cojo-wind 10000  --cojo-slct --out ${wd}Gender.Shared.chr23

################step3.2 对Female的结果进行conditional & joint association analysiswd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/"
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/08.1.XWAS.Meta.Conditional.analysis/"
gcta --bfile ${wd}Female.Imputed --maf 0.01 --cojo-file ${wd}Female.COJO.input.files.txt  --cojo-p 0.0001 --update-sex ${wd}Updata.Sex.txt --cojo-wind 10000  --cojo-slct --out ${wd}Female.Stratified.chr23


