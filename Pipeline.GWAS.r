

wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/ChinaMap.Imputation.Latest/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Important.file/"
wd3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/step1.Batch.GWAS/"
wd4="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Qsub.tasks/"
wd5="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Raw.Data/1.VCF.QCed.R2greater0.6/"
wd6="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Raw.Data/2.PLINK.UpdatedID.after.VCFQCed/"
wd7="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Raw.Data/3.PLINK.QCed.after.UpdatedID/"
wd8="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/GWAS.OF.ESRD/Final.Results/"


wdhg38="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg38.data/"
wd9="/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GWAS/Impute.Score/"
wd10="/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GWAS/Batch.Shell.Code/"
wd11="/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GWAS/Important.Files/"
wd12="/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GWAS/Raw.Data/QCed.PLINKdata.hg38.filtered.low.ER2/"
wd13="/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GWAS/Raw.Data/QCed.vcf.UpdatedID.hg38/"



path="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data/"
path1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center1.ESRD.CKD/"
path2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center2.GDPH.DM/"
path3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center3.GDPH.Health/"
path4="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center4.Zhongzhong.Control/"
path5="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center5.Shandong.Control/"



dir1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center1.ESRD.CKD/"
dir2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center2.GDPH.DM/"
dir3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center3.GDPH.Health/"
dir4="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center4.Zhongzhong.Control/"
dir5="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center5.Shandong.Control/"


#####mkdir 上述目录
for i in {1..5}; do
  eval dir=\$dir$i
  mkdir -p "$dir"
done


####################Alignment Strand
#####center1.ESRD.CKD
for chr in {1..22}; do
    plink --bfile ${path}ESRD_forward --chr $chr --make-bed --out ${path1}chr$chr
done;

######center2.GDPH.DM      
for chr in {1..22}; do
    plink --bfile ${path}GDPH.DM.Control --chr $chr --make-bed --out ${path2}chr$chr
done;

######center3.GDPH.Health
for chr in {1..22}; do
    plink --bfile ${path}GDPH.Control --chr $chr --make-bed --out ${path3}chr$chr
done;

######center4.Zhongzhong.Control
for chr in {1..22}; do
    plink --bfile ${path}Zhongzhong.Control --chr $chr --make-bed --out ${path4}chr$chr
done;

######center5.Shandong.Control
for chr in {1..22}; do
    plink --bfile ${path}Shandong.Control --chr $chr --make-bed --out ${path5}chr$chr
done;

########################################Alignment Strand;
#####https://databeauty.com/blog/tutorial/2017/02/20/GWAS-prephasing-and-imputation.html
###https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer

jar="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/GenotypeHarmonizer-1.4.25-SNAPSHOT/"
Ref="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/1000G_301CHN/"


###########Conver Reference panel to VCF files
for i in {1..22};
do
/share/app/bcftools-1.4/bin/bcftools convert --haplegendsample2vcf ${Ref}CHN.chr$i.hap,${Ref}CHN.chr$i.legend,${Ref}CHN.chr$i.samples  -o ${Ref}CHN.chr$i.vcf
bgzip ${Ref}CHN.chr$i.vcf
tabix -p vcf ${Ref}CHN.chr$i.vcf.gz
done;





################### flip  j reprensts five center, i represents chromosomes
# 定义路径数组
paths=(
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center1.ESRD.CKD/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center2.GDPH.DM/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center3.GDPH.Health/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center4.Zhongzhong.Control/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr/center5.Shandong.Control/"
)

dirs=(
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center1.ESRD.CKD/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center2.GDPH.DM/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center3.GDPH.Health/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center4.Zhongzhong.Control/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/00.Raw.data.split.chr.Aligment.Strand/center5.Shandong.Control/"
)

jar="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/software/GenotypeHarmonizer-1.4.25-SNAPSHOT/"
Ref="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/1000G_301CHN/"


#####################################Flip the data ***** Important action 遍历数组元素
for j in {0..4}; do
    for i in {1..22}; do
        input_file="${paths[j]}chr$i"
        output_file="${dirs[j]}chr$i.Aligned"
        java -Xmx60g -jar ${jar}GenotypeHarmonizer.jar --inputType PLINK_BED  --input "$input_file" --outputType PLINK_BED --output "$output_file"  --refType VCF --ref ${Ref}CHN.chr${i}
    done
done





##############################Checking unrunning file and re-runing
for j in {0..4}; do
      echo  ${dir[j]}chr*.Aligned* | grep -oE '[0-9]+' | sort -n > ${dir[j]}existing_numbers.txt
  done



for j in {0..4}; do
  fail=`uniq -c ${dir[j]}existing_numbers.txt | awk '$1 < 5 {print $2}'| tr '\n' ' '`
  echo ${fail}
  for i in ${fail}; do
  java -Xmx40g -jar ${jar}GenotypeHarmonizer.jar --inputType PLINK_BED  --input ${paths[j]}chr${i}  --outputType PLINK_BED --output ${dir[j]}chr${i}.Aligned  --refType VCF --ref ${Ref}CHN.chr${i}
   done
done 




#####Merge the plink data of chr1-chr22 of 5 centers
for j in {0..4}; do
    for i in {2..22}; do 
        plink --bfile ${dirs[j]}chr$i.Aligned --recode vcf --out ${dirs[j]}chr$i.Aligned
    done
done



for j in {0..4}; do
    echo "" > ${dirs[j]}tmp.list;
    for i in {2..22}; do 
        echo "${dirs[j]}chr$i.Aligned.bed ${dirs[j]}chr$i.Aligned.bim ${dirs[j]}chr$i.Aligned.fam" >> ${dirs[j]}tmp.list;
    done
done

for j in {0..4}; do
        chr1_file="${dirs[j]}chr1.Aligned"
        list="${dirs[j]}tmp.list"
        h=$((j + 1))
        plink --bfile  "$chr1_file" --allow-no-sex --noweb --merge-list "$list" --make-bed --out  ${dirs[j]}All.Merge.Center$h
done


for j in {0..4}; do
        h=$((j + 1))
        plink --bfile  ${dirs[j]}All.Merge.Center$h --recode vcf --out ${dirs[j]}All.Merge.Center$h
    done;


###############################################Flip the data
works=(
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/center1.ESRD.CKD/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/center2.GDPH.DM/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/center3.GDPH.Health/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/center4.Zhongzhong.Control/"
    "/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/center5.Shandong.Control/"
)

for j in {0..4}; do
 mkdir ${works[j]}
done



#############Step1. hg37 转为hg38, 从https://zhuanlan.zhihu.com/p/383252096上找到下载hg38.fa.gz的方式，并通过GATK建立hg38.dict
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/"
#建立字典
#java -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/picard-2.18.16/picard.jar CreateSequenceDictionary  R=${wd}Important.files/hg38.fa  O=${wd}Important.files/hg38.dict

#转换坐标，参考https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_tools/liftovervcf.html；

hg38="/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference_bigdatacompute/hg38_noalt_withrandom/"
  chain="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/Important.files/"
  Code="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Raw.data.hg38.Strand/Batch.Code/"
  picard="/share/app/picard/2.23.8/"

for j in {0..4}; do
   h=$((j + 1))
   for i in {1..22}; do
  echo "java -jar -Xmx200g ${picard}picard.jar LiftoverVcf I=${dirs[j]}chr$i.Aligned.vcf    O=${works[j]}chr$i.Aligned.hg38.vcf       CHAIN=${chain}b37ToHg38.over.chain  REJECT=${works[j]}rejected_variants.${i}.vcf R=${hg38}hg38.fa RECOVER_SWAPPED_REF_ALT=true" > ${Code}Center$h.Hg37toHg38.Chr${i}.sh
  done
done

for j in {0..4}; do
   h=$((j + 1))
   for i in {1..22};do
    nohup sh ${Code}Center$h.Hg37toHg38.Chr${i}.sh> ${Code}Center$h_output_${i}.log 2>&1 &
   done
done



###############################IMputation (ChinaMAP Imputation Server is a genotype imputation server utilizing the ChinaMAP reference panel constructed from the China Metabolic Analytics Project)
################(https://www.nature.com/articles/s41422-021-00564-z)



###########QC: INFO/DR2>=0.6 && INFO/MAF>=0.01
for i in {1..22};
do echo "/share/app/bcftools-1.4/bin/bcftools filter -i 'INFO/R2>=0.6 && INFO/MAF>=0.01' ${wd}chr${i}/chr${i}.dose.vcf.gz  -Ov -o ${wd5}/chr${i}.dose.R2above0.6.vcf" > ${wd3}Filter.chr${i}.sh; 
echo " qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd3}Filter.chr${i}.sh" > ${wd4}Filter.chr${i}.qsub.sh;
 cd ${wd4};
 sh ${wd4}Filter.chr${i}.qsub.sh; 
 done;

##########Updated IID
for i in {1..22};
do echo "plink --vcf ${wd5}/chr${i}.dose.R2above0.6.vcf --make-bed --double-id --update-ids ${wd2}IID.Update.txt --out ${wd6}/chr${i}.dose" > ${wd3}UpdatedID.chr${i}.sh;
echo " qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd3}UpdatedID.chr${i}.sh" > ${wd4}UpdatedID.chr${i}.qsub.sh; 
cd  ${wd4}; 
sh ${wd4}UpdatedID.chr${i}.qsub.sh;
done;

##########Caculate Impute score: ER2 for TYPED (filter out genotyped SNPs with low quality (empirical ER2<0.7)14 96 and re-impute;(https://www.biorxiv.org/content/10.1101/2021.10.19.464854v1.full.pdf))
echo "" > ${wd9}Impute.Score.txt



for i in {1..22};do
echo " grep -v '^#' ${wd5}/chr${i}.dose.R2above0.6.vcf | awk -F \"\t\" '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8}' > ${wd9}/Impute.Score.chr${i}.txt " >  ${wd10}Impute.Score.chr${i}.sh
echo " qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd10}Impute.Score.chr${i}.sh" >  ${wd10}Impute.Score.chr${i}.qsub.sh
cd  ${wd10}
sh ${wd10}Impute.Score.chr${i}.qsub.sh
done;

echo "" > ${wd11}Impute.Score.txt
cat ${wd9}/Impute.Score.chr*txt >> ${wd11}Impute.Score.txt

sed  '1iCHROM POS ID REF ALT QUAL FILTER INFO' ${wd11}Impute.Score.txt > ${wd11}Impute.Score.vcf


sed -i "s/\;/ /g" ${wd11}Impute.Score.vcf

grep "ER2=*" ${wd11}Impute.Score.vcf > ${wd11}TYPED.SNPs.vcf

sed -i "s/ER2=/ER2= /g" ${wd11}TYPED.SNPs.vcf

awk '{print NF; exit}' ${wd11}TYPED.SNPs.vcf

awk -F" " '$12  < 0.7' ${wd11}TYPED.SNPs.vcf | awk -F " " '{print $1":"$2}' > ${wd11}Low.Quality.Typed.SNPs.txt

sed  "s/ \+/;/g" ${wd11}Impute.Score.vcf > ${wd11}tmp.Score.vcf;
sed -i "s/=/;/g" ${wd11}tmp.Score.vcf;
awk -F ";"  '{print $3,$13}' ${wd11}tmp.Score.vcf > ${wd2}Impute.Score.txt; ####Extract column SNP and R2
rm ${wd11}tmp.Score.vcf;
sed -i "1iID R2" ${wd2}Impute.Score.txt;
less -S ${wd2}Impute.Score.txt;


#########Further QC: hwe 10e-4; Call rate; missing value
for i in {1..22};
do echo "plink --bfile ${wd6}/chr${i}.dose --hwe 0.0001 --mind 0.01 --maf 0.01 --geno 0.05 --make-bed --out ${wd7}chr${i}.dose.QCed" > ${wd3}QC.chr${i}.sh;
echo " qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd3}QC.chr${i}.sh" > ${wd4}QC.chr${i}.qsub.sh;
cd  ${wd4};
sh ${wd4}QC.chr${i}.qsub.sh;
done;


###########################Step2. Differnetial analysis to identify SNPs associated with phenotype

######Calculate covariants: 10 PCs
sed "" > ${wd2}Chr.list;
for i in {1..22};do echo "${wd7}chr${i}.dose.QCed" >> ${wd2}Chr.list; done;
 plink --merge-list ${wd2}Chr.list --make-bed --out ${wd7}Merged.ALL --allow-no-sex;
 plink --bfile ${wd7}Merged.ALL  --freq --out ${wd2}snp_freq;
 plink --bfile ${wd7}Merged.ALL --indep-pairwise 1500 150 0.2  --maf 0.01 --make-bed  --out ${wd7}ALL.Pruned 
 plink --bfile ${wd7}Merged.ALL --extract ${wd7}ALL.Pruned.prune.in  --make-bed --out ${wd7}Merged.ALL.Pruned
plink --bfile ${wd7}Merged.ALL.Pruned --pca 10 --out ${wd2}Merged.ALL.PCA.txt
##cp  ${wd2}GWAS.pheno.all.3794.txt  ${wd2}GWAS.pheno.all.3794.Backup.txt ###${wd2}GWAS.pheno.all.3794.Backup.txts是最原始的表型文件，不包含genomics PC信息
sh  ${wd2}Format.Pheno.file.sh
less -S ${wd2}GWAS.pheno.all.3794.txt

##################Defferential analysis
for i in {1..22};
do echo "plink --bfile ${wd7}chr${i}.dose.QCed --pheno ${wd2}GWAS.pheno.all.3794.txt  --pheno-name Group --covar ${wd2}GWAS.pheno.all.3794.txt --covar-name Age,Gender,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10  --allow-no-sex --assoc  --logistic --out ${wd8}chr${i}.Result" > ${wd3}GWAS.chr${i}.sh; 
echo " qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd3}GWAS.chr${i}.sh" > ${wd4}GWAS.chr${i}.qsub.sh; cd  ${wd4}; 
sh ${wd4}GWAS.chr${i}.qsub.sh; 
done;


#################Formatting the GWAS result
 echo "" > ${wd8}ESRD.GWAS.Results.txt;
  for i in {1..22};
  do grep "ADD" ${wd8}chr${i}.Result.assoc.logistic >> ${wd8}ESRD.GWAS.Results.txt; 
  done;

sed -i "s/\"//g" ${wd8}ESRD.GWAS.Results.txt
sed -i "s/ \+/ /g" ${wd8}ESRD.GWAS.Results.txt
sed -i '/^ CHR/d' ${wd8}ESRD.GWAS.Results.txt
sed -i '1iCHR SNP BP A1 TEST NMISS OR STAT P'  ${wd8}ESRD.GWAS.Results.txt
sed -i '2d' ${wd8}ESRD.GWAS.Results.txt
sed -i 's/^ \+//1' ${wd8}ESRD.GWAS.Results.txt
mv ${wd8}ESRD.GWAS.Results.txt ${wd8}ESRD.GWAS.org.Results.txt
cp ${wd8}ESRD.GWAS.org.Results.txt ${wd2}ESRD.GWAS.Results.txt
#####If the R2 of results of all SNPs are greater than 0.6 instead of 0.7, then excute the following code:
sh ${wd2}Filter.R2.0.7.SNPs.sh ##Save the result of SNPs (R2 > 0.7)


#################### Step4. Conditional & joint association analysis (COJO): to identify independent associations
sh  ${wd2}Format.COJO.input.file.sh

qsub -cwd -l vf=15G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd2}COJO.sh


#####Calculation LD between lead snp and other snps

##plink --bfile chr2.dose.QCed --r2 --ld-snp 2:24470601:G:C --ld-window-kb 5000 --ld-window 100000 --ld-window-r2 0.6 --out tmp.0.6



