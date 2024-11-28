##############################以下脚本是根据risd号注释到cytoband请参考网站https://mp.weixin.qq.com/s/-jrwxJfazYoTXTRNip3bGQ。


####下载cytoBand 数据库
wd="/hwfssz5/ST_HEALTH/P20Z10200N0041/huyuxuan1/Software/annovar/";
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/humandb/"
perl ${wd}annotate_variation.pl --downdb --buildver hg38 cytoBand ${wd2}

######下载hg38的SNP数据:参考网站：https://annovar.openbioinformatics.org/en/latest/user-guide/download/#additional-databases
perl ${wd}annotate_variation.pl -buildver hg38 -downdb -webfrom annovar snp150 ${wd2}

#####将rsid 转换为可以annovar软件的输入格式
perl ${wd}convert2annovar.pl -format rsid /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/Meta.SNPlist.txt -dbsnpfile ${wd2}hg38_avsnp150.txt > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/XWAS.Meta.snplist.avinput


############进行cytoband注释
wd="/hwfssz5/ST_HEALTH/P20Z10200N0041/huyuxuan1/Software/annovar/";
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/humandb/"

perl ${wd}annotate_variation.pl -regionanno -dbtype cytoBand -out Meta.SNP.cytoBand -buildver hg38 /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/DBSNP/XWAS.Meta.snplist.avinput  ${wd2}

