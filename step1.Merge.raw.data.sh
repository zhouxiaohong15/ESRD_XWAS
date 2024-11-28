############Step1 首先要翻转省医的1160个样本的数据"ESRD_forward.bed"


wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data/Flip.Raw.data/"
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data/"

#########Step2. 合并省医的ESRD和所有对照组; 数据组成为：省医1600多个ESRD+CKD，中山医肿瘤547个肿瘤对照，山东500个健康对照，省医后面再匹配的900个对照;省医185个糖尿病样本; 总共3806个样本

#1）和省医900个对照合并
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}ESRD_forward --bmerge ${wd1}GDPH.Control  --make-bed  --allow-no-sex --out ${wd1}GDPH.Raw.ASA;

wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Merge.Raw.data"
##2）和省医174个DM对照合并
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}GDPH.Raw.ASA --bmerge ${wd1}GDPH.DM.Control --make-bed  --allow-no-sex --out ${wd1}All.GDPH.ESRD.CKD.DM.Raw.ASA;

##3）和山东500个对照合并
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile All.GDPH.ESRD.CKD.DM.Raw.ASA --bmerge Shandong.Control --make-bed  --allow-no-sex --out GDPH.Shandong;
#出现报错，因此再翻转一次后:
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile All.GDPH.ESRD.CKD.DM.Raw.ASA --flip GDPH.Shandong-merge.missnp --make-bed --allow-no-sex --out All.GDPH.ESRD.CKD.DM.Raw.ASA.exclude.flip

###翻转后这里老是报错出现“Error: 11 variants with 3+ alleles present.”
#因此我直接去掉了   --exclude GDPH.Shandong-merge.missnp ;


/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile All.GDPH.ESRD.CKD.DM.Raw.AS.flip --exclude GDPH.Shandong-merge.missnp --make-bed --allow-no-sex --out All.GDPH.ESRD.CKD.DM.Raw.ASA.exclude

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile All.GDPH.ESRD.CKD.DM.Raw.ASA.exclude --bmerge Shandong.Control --make-bed  --allow-no-sex --out GDPH.Shandong;

##4) 和中山医肿瘤500个对照合并
#出现报错“Error: 290004 variants with 3+ alleles present.”
#####解决办法： --flip ${wd2}All.Merge.ASA-merge.missnp；担依然报错 Error: 273 variants with 3+ alleles present； 解决办法： -- exclude ${wd2}All.Merge.ASA-merge.missnp; 最后再merge

wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Merge.Raw.data/"
cp ${wd2}All.Merge.ASA-merge.missnp  ${wd1}All.Merge.ASA-merge.missnp
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}GDPH.Shandong --flip ${wd2}All.Merge.ASA-merge.missnp -make-bed --allow-no-sex --out  ${wd1}GDPH.Shandong.flip
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}GDPH.Shandong.flip --exclude ${wd2}All.Merge.ASA-merge.missnp -make-bed --allow-no-sex --out  ${wd1}GDPH.Shandong.flip.exclude

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}GDPH.Shandong.flip.exclude --bmerge ${wd1}Zhongzhong.Control --make-bed  --allow-no-sex --out ${wd2}All.Merge.ASA

