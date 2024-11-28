#########该脚本是批量生成共定位分析R脚本的,分为step1 和step2两个步骤

##################Step1 从所有染色体中提取chrX的数据
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/GTEX.eQTL.all.assoc/"

ls ${wd}*allpairs.txt > ${wd}eQTL.file.name.txt; #####输出该路径下所有组织的eQTL数据文件名

sed -i "s|${wd}||g" ${wd}eQTL.file.name.txt; #### 去掉路径前缀

while IFS= read -r line; 
do filename="$line"
   new_filename="${filename%%.*}"
   grep ".*chrX_.*_b38.*" ${wd}${filename} > ${wd}${new_filename}.chrX.eQTL.hg38.txt
     done < ${wd}eQTL.file.name.txt;



##############Step2 以下是将所有组织的名称存入数组arr中，然后后续批量替换生成不同组织的共定位分析R脚本


awk -F "." '{print $1}' ${wd}eQTL.file.name.txt > ${wd}Tissue.txt; ######将所有感>兴趣的组织名称输入到文件Tissue.txt中。


########以下是将所有组织的名称存入数组arr中，方便后续批量替换生成不同>组织的共定位分析R脚本
declare -a arr=()
mapfile -t arr < <(cat ${wd}Tissue.txt)

########################################接下来是批量生成不同组织共定位分析的R脚本
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/01.Batch.Rscript/"
cd ${wd2};
rm *sh.e* *sh.o*;

for i in "${arr[@]}";
do sed  "s/Tissue/$i/g"  ${wd1}step16.Nearest.Gene.Coloc.analysis1.R > ${wd2}$i.Coloc.R;
   echo "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/anaconda3/bin/Rscript ${wd2}$i.Coloc.R" > ${wd2}$i.Coloc.sh;
    cd ${wd2};
    qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd2}$i.Coloc.sh;
done;

wd4="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/";
cd /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/12.Nearest.Genes.For.Coloc/Coloc.Results/;

rm ${wd4}tmp.txt;
cat ./*Coloc.Result.txt >> ${wd4}tmp.txt;
sed -i "s/\"//g" ${wd4}tmp.txt;
awk 'NF > 0' ${wd4}tmp.txt > ${wd4}ALL.coloc.final.Results.txt;
sed -i "/nsnps/d" ${wd4}ALL.coloc.final.Results.txt;
sed -i "1i nsnps PP.H0.abf PP.H1.abf PP.H2.abf PP.H3.abf PP.H4.abf gene tissue class SNP" ${wd4}ALL.coloc.final.Results.txt;
