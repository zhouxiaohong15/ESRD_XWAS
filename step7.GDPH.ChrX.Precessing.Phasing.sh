
#########################################以下脚本是对于未质控的原始的X染色体:1)转换到hg38坐标 2)haploid转换为diploid 3)提取男性和女性的NPR区域（ASA芯片针对X染色体上的位点全部在NPR区域） 4)Phasing
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Merge.Raw.data/"

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink  --bfile ${wd}All.Merge.ASA --chr 23 --missing --impute-sex  --make-bed --out /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/ChromosomeX


wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg37.data/";
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg37.data/";

wd3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/"

/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink --bfile /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/ChromosomeX  --flip ${wd3}ASA.BOT-Pos.TOP-Neg.txt  --make-bed --out /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/chrX.Flip;

###############转为hg37坐标的vcf(很关键：这样才能使得后续有更高的转换成功率）
/hwfssz4/BC_PUB/Software/03.Soft_ALL/plink-2.0/plink2 --bfile /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/chrX.Flip --recode vcf-4.2 --out /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/chrX.hg37.Flip  --ref-from-fa --fa /ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference/human_g1k_v37/human_g1k_v37.fasta

awk -F'\t' '$4 != "D" && $4 != "I"' /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/chrX.hg37.Flip.vcf > /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/chrX.hg37.Flip.new.vcf

mv /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/chrX.hg37.Flip.new.vcf /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/chrX.hg37.Flip.vcf

############hg37转换为hg38
java -jar -Xmx100g /hwfssz4/BC_PUB/Software/03.Soft_ALL/picard-2.18.16/picard.jar LiftoverVcf I=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/chrX.hg37.Flip.vcf O=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/chrX.hg38.Flip.vcf REJECT=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/rejected_variants.chrX.vcf CHAIN=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/Important.files/b37ToHg38.over.chain R=/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference_bigdatacompute/hg38_noalt_withrandom/hg38.fa

#####转换为vcf.gz的格式
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -Oz -o chrX.hg38.Flip.vcf.gz chrX.hg38.Flip.vcf
###建立索引
/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf chrX.hg38.Flip.vcf.gz

#########################将文件统一变为diploid的形式
###### Convert haploid genotypes to homozygous diploids:
###Often chrX is represented as haploid genotypes for males, however, Beagle can only handle diploid genotypes.
### Fix the chromosome X ploidy to phased diploid
### Requires a ploidy.txt file containing 
### space-separated CHROM,FROM,TO,SEX,PLOIDY 
echo "chrX 1 156040895 M 2" > ploidy.txt
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools +fixploidy chrX.hg38.Flip.vcf.gz  -Ov -- -p ploidy.txt | sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | sed 's#1/0#1\|0#g;s#0/1#0\|1#g' | bcftools view -Oz -o  chrX.hg38.Flip.diploid.vcf.gz



####################提取 male-NPAR ,female-NPAR区域（不提取PAR区域，因为ASA芯片测的全是NPR区域）


###########female的NPR区域
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -S GDPH.Female.list chrX.hg38.Flip.diploid.vcf.gz  -Oz -o Female.chrX.hg38.NPR.vcf.gz;


############male的NPR区域
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -S GDPH.Male.list chrX.hg38.Flip.diploid.vcf.gz  -Oz -o Male.chrX.hg38.NPR.vcf.gz;


#################Phasing using Eagle
#######女性的数据进行phasing
/hwfssz4/BC_PUB/Software/03.Soft_ALL/Eagle-2.4.1/eagle          --vcf /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Female.chrX.hg38.NPR.vcf.gz --chrom chrX      --geneticMapFile /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/Genetic.Map.Eagle/genetic_map_hg38_withX.txt            --numThreads=8            --Kpbwt=20000            --outPrefix /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Female.chrX.hg38.NPR.Phasing.vcf.gz

########男性的数据进行phasing

/hwfssz4/BC_PUB/Software/03.Soft_ALL/Eagle-2.4.1/eagle          --vcf /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Male.chrX.hg38.NPR.vcf.gz --chrom chrX      --geneticMapFile /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/Genetic.Map.Eagle/genetic_map_hg38_withX.txt            --numThreads=8            --Kpbwt=20000            --outPrefix /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/Male.chrX.hg38.NPR.Phasing.vcf.gz


