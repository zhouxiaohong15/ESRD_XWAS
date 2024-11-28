cd /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/

###################Reference Panel: 是最新的1000GP hg38 : 参考网页：https://www.internationalgenome.org/data-portal/data-collection 和 https://pubmed.ncbi.nlm.nih.gov/36055201/。

##########Female
java -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/beagle-5.0/beagle.03Jul18.40b.jar gt=Female.chrX.hg38.NPR.Phasing.vcf.gz  ref=1000GP.ALL.Female.ChrX.NPR.final.vcf.gz  out=GDPH.Female.chrX.Imputed nthreads=16  ne=14269  impute=true

/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf GDPH.Female.chrX.Imputed.vcf.gz

#####Male

java -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/beagle-5.0/beagle.03Jul18.40b.jar gt=Male.chrX.hg38.NPR.Phasing.vcf.gz  ref=1000GP.ALL.Male.ChrX.NPR.final.vcf.gz  out=GDPH.Male.chrX.Imputed nthreads=16  ne=14269  impute=true

/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf GDPH.Male.chrX.Imputed.vcf.gz

#############################step2. QC: ER<0.3 以及MAF < 0.01 ( the cutoff of 0.3 is commonly used for post-imputation SNP filtering).
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools filter -i 'INFO/DR2>=0.3 && INFO/MAF>=0.01' GDPH.Male.chrX.Imputed.vcf.gz  -Oz -o GDPH.Male.chrX.Imputed.QCed.vcf.gz

/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools filter -i 'INFO/DR2>=0.3 && INFO/MAF>=0.01' GDPH.Female.chrX.Imputed.vcf.gz  -Oz -o GDPH.Female.chrX.Imputed.QCed.vcf.gz

################################去除中山医肿瘤和山东的对照样本（该step忽略）
bcftools view  -S ^/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/female.txt  GDPH.Female.chrX.Imputed.QCed.vcf.gz -Oz -o GDPH.Female.chrX.Imputed.QCed01.vcf.gz;
mv GDPH.Female.chrX.Imputed.QCed01.vcf.gz GDPH.Female.chrX.Imputed.QCed.vcf.gz;

bcftools view  -S ^/hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/male.txt  GDPH.Male.chrX.Imputed.QCed.vcf.gz -Oz -o GDPH.Male.chrX.Imputed.QCed01.vcf.gz;
mv GDPH.Male.chrX.Imputed.QCed01.vcf.gz GDPH.Male.chrX.Imputed.QCed.vcf.gz;


#########################转换为plink格式
plink --vcf GDPH.Male.chrX.Imputed.QCed.vcf.gz --remove /hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/male.txt --make-bed --double-id --out GDPH.Male.chrX.Imputed;

plink --vcf GDPH.Female.chrX.Imputed.QCed.vcf.gz  --make-bed --double-id --out GDPH.Female.chrX.Imputed;
plink --bfile GDPH.Female.chrX.Imputed --update-ids update.female.ids.txt --remove /hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/female.txt --make-bed --double-id --out GDPH.Female.chrX.Imputed;

##########################QC (男性不需要做hwe质控）
plink --bfile GDPH.Female.chrX.Imputed --hwe 0.000001 --mind 0.05 --MAF 0.05 --geno 0.05 --make-bed --out Female.Imputed

plink --bfile GDPH.Male.chrX.Imputed --make-bed --out Male.Imputed


#########################差异分析：XWAS

cp  /hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GDPH.Male.Pheno.new.txt  /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/


####PCA分析
plink --bfile Male.Imputed   --chr 23 --indep-pairwise 1500 150 0.2  --maf 0.01 --make-bed  --out Male.Imputed.Pruned

plink --bfile Male.Imputed --extract Male.Imputed.Pruned.prune.in  --make-bed --out Male.Imputed.Pruned


plink --bfile Male.Imputed.Pruned  -chr-set 23 --pca 10 --out Male.PCA ###额外的参数-chr-set 23

plink --bfile Male.Imputed --maf 0.05 --make-bed --out  Male.Imputed


####整理表型文件
awk -F ' ' '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' Male.PCA.eigenvec > PC.txt
sed -i '1i PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10' PC.txt
awk '{print $1,$2,$3,$4}' GDPH.Male.Pheno.new.txt > male.4col.txt
paste male.4col.txt PC.txt > all.Male.pheno.txt
sed -i "s/ /\t/g" all.Male.pheno.txt

plink --bfile Male.Imputed --pheno all.Male.pheno.txt  --pheno-name Group --covar all.Male.pheno.txt --covar-name Age  --allow-no-sex --assoc  --logistic --out ./Male.Result;

awk '$5=="ADD" {print}' Male.Result.assoc.logistic > Male.final.results.txt; sed -i '/NA/d' Male.final.results.txt;
sort -k9,9 -g Male.final.results.txt > Male.output.txt;
sed -i '1i CHR SNP BP A1 TEST NMISS OR STAT P' Male.output.txt; sed -i "s/ \+/ /g" Male.output.txt;


######################女性
cp  /hwfssz5/ST_HEALTH/P20Z10200N0041/zhouxiaohong/GDPH.Female.Pheno.new.txt  /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/


####PCA
plink --bfile Female.Imputed   --chr 23 --indep-pairwise 1500 150 0.2  --maf 0.01 --make-bed  --out Female.Imputed.Pruned

plink --bfile Female.Imputed --extract Female.Imputed.Pruned.prune.in  --make-bed --out Female.Imputed.Pruned


plink --bfile Female.Imputed.Pruned  -chr-set 23 --pca 10 --out Female.PCA ###额外的参数-chr-set 23

plink --bfile Female.Imputed --maf 0.05 --make-bed --out  Female.Imputed

####整理表型文件

awk -F ' ' '{print $3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' Female.PCA.eigenvec > PC.txt
sed -i '1i PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10' PC.txt
awk '{print $1,$2,$3,$4}' GDPH.Female.Pheno.new.txt > Female.4col.txt
paste Female.4col.txt PC.txt > all.Female.pheno.txt
sed -i "s/ /\t/g" all.Female.pheno.txt


plink --bfile Female.Imputed --pheno all.Female.pheno.txt  --pheno-name Group --covar all.Female.pheno.txt --covar-name Age  --allow-no-sex  --logistic --out ./Female.Result

awk '$5=="ADD" {print}' Female.Result.assoc.logistic > Female.final.results.txt; sed -i '/NA/d' Female.final.results.txt;
sort -k9,9 -g Female.final.results.txt > Female.output.txt;
sed -i '1i CHR SNP BP A1 TEST NMISS OR STAT P' Female.output.txt; sed -i "s/ \+/ /g" Female.output.txt;

