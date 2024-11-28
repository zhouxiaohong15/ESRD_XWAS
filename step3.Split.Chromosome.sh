wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg37.data/";
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg37.data/";



####Step1.质控后的plink 转为vcf
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}All.QCed.ASA --recode vcf --out ${wd1}All.QCed.ASA.hg37;

#上述代码发出警告：Warning: At least one VCF allele code violates the official specification;
#other tools may not accept the file.  (Valid codes must either start with a
#'<', only contain characters in {A,C,G,T,N,a,c,g,t,n}, be an isolated '*', or
#represent a breakend.) 因此，我们必须将第四列对应的非上述字符的行删除：

awk -F'\t' '$4 != "D" && $4 != "I"' ${wd1}All.QCed.ASA.hg37.vcf > ${wd1}outputfile.txt;
mv ${wd1}outputfile.txt ${wd1}All.QCed.ASA.hg37.vcf;


################Step2. 将/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg37.data/All.QCed.ASA.hg37.vcf转换为PLINK形式,该vcf文件是经过了去除“I" "D"碱基的

#########All.QCed.hg37.recode对应的plink文件和All.QCed.hg37对应的plink文件区别是：去除了“I" "D">碱基的行，共400行

/hwfssz4/BC_PUB/Software/03.Soft_ALL/vcftools-0.1.16/bin/vcftools  --vcf ${wd1}All.QCed.ASA.hg37.vcf --plink  --out ${wd1}All.QCed.hg37.recode;

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --file ${wd1}All.QCed.hg37.recode  --make-bed --out ${wd1}All.QCed.hg37.recode;


########################Step3. Flip翻转：将BOT + ,TOP - 的位点翻转过来(参考文件：/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/ASA.Annotation.Files.txt；重点对应的是IlmnStrand和RefStrand两列。）
wd3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/"
/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}All.QCed.hg37.recode --flip ${wd3}ASA.BOT-Pos.TOP-Neg.txt  --make-bed --out ${wd1}All.QCed.hg37.recode.Flip;


##########################################Step4.转为vcf: 一定要用plink2转换为vcf
/hwfssz4/BC_PUB/Software/03.Soft_ALL/plink-2.0/plink2 --bfile ${wd1}All.QCed.hg37.recode.Flip  --recode vcf-4.2 --out ${wd1}All.QCed.hg37.recode.Flip --ref-from-fa --fa /ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference/human_g1k_v37/human_g1k_v37.fasta


###############Step5. 将${wd1}All.QCed.hg37.recode对应的PLink文件分割为22个染色体，注意：为防止链flip失败，因此用${wd1}All.QCed.hg37.recode而不是${wd1}All.QCed.hg37.recode.flip，来重新进行flip。

for chr in {1..23} 
 do /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile ${wd1}All.QCed.hg37.recode --chr $chr --flip /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/ASA.BOT-Pos.TOP-Neg.txt  --make-bed --out ${wd}All.QCed.hg37.Chr$chr
done;

##################Step5. 将上述的hg37坐标的22个染色体plink文件转换为vcf，一定要用PLINK2 --ref-from-fa, 不然会导致最终hg37转换到hg38坐标的时候，转换率过低<10%
for chr in $(seq 1 23)
  do /hwfssz4/BC_PUB/Software/03.Soft_ALL/plink-2.0/plink2 --bfile ${wd}All.QCed.hg37.Chr$chr --recode vcf-4.2 --out ${wd}All.QCed.hg37.Chr$chr --ref-from-fa --fa /ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference/human_g1k_v37/human_g1k_v37.fasta
done;

#####表头加上这一行 不然会报错
sed -i  "4i ##INFO=<ID=PR,Number=0,Type=Flag,Description=\"Provisional reference allele, may not be based on real reference genome\"> " /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg37.data/All.QCed.hg37.Chr*.vcf 


########有些文件有上述重复行，会报错，删除其中一行
for chr in  2 3 4 5 6 7 8 9 10 11 13
do sed -i '4d' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg37.data/All.QCed.hg37.Chr${chr}.vcf
done;
