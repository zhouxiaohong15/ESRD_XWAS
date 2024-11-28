
##########################以下脚本是对reference panel进行预处理，从而更快更好的进行填补分析，包括以下几个方面：1)检查是否有multiallelic sites，检查Sample ID 是否正确；2)将ploidy转换为diploidy 3)提取东亚人的数据 4)Multiple processing commands piped together (见下面的参数) 5)去除重复的SNP ID 6)提取1000GP的NPR区域，并分为男性和女性子集
####可参考网页https://www.protocols.io/run/genotype-imputation-workflow-v3-0-xbgfijw?step=11

#### Check for multiallelic sites

/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bin/bcftools norm -m -any /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz -Oz -o  /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/panel_phased_split_multiallelic_chrX.vcf


###Generate a list of the reference panel sample IDs
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/"

/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bin/bcftools  query -l ${wd}CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz > ${wd}panel_sample_IDs.txt


mv CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.v2.vcf.gz  1000GP.ChrX.ALL.hg38.panel.vcf.gz
/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf ${wd}1000GP.ChrX.ALL.hg38.panel.vcf.gz

# Multiple processing commands piped together
#The reference panel files should contain: 
#1.Remove the rare variants, here singletons and doubletons by setting AC threshold with 'bcftools view'.
#3.Keep only SNPs with 'bcftools view'. Here, the 1000GP data included a tag VT in the INFO field and data contain also structural variants which should be excluded.
#4.Align the variants to reference genome with 'bcftools norm' in order to have the REF and ALT alleles in the shortest possible representation and to confirm that the REF allele matches the reference genome, additionally remove duplicate variants (-d none).
#5.After alignment, remove multiallelic records with 'bcftools view', since these are formed during the alignment if the REF does not match with the reference genome.
#6.Finally, remove sites containing missing data with 'bcftools view'.
    /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -e 'INFO/AC<3 | INFO/AN-INFO/AC<3' 1000GP.ChrX.ALL.hg38.panel.vcf.gz -Ov |     /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools norm -f /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Impotant.files/hg38.fa -d none -Ov |       /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -m 2 -M 2 -v snps  -Ov | /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -g ^miss -Oz -o 1000GP.ALL.ChrX.QCed.vcf.gz

wd="jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/"
/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf ${wd}1000GP.ALL.ChrX.QCed.vcf.gz



################### Duplicate ID removal,这个步骤省略先 (本文件从官网下载下来就没有rsID号，因此考虑自己生成ID号：chr:BP:ALT:Ref)
#####:Remove duplicate IDs. If you wish to preserve all multiallelic sites, replace the ID column with a unique ID e.g. CHR_POS_REF_ALT (参考上述生成的1000GP_SNPID_chrX.vcf.gz)
#######Here, 1000GP did not contain multiallelic sites after AC filtering, and thus, RSIDs were preserved in the ID column. And since RSIDs are not always unique, duplicates should be removed.(参考上述生成的1000GP_chrX.ploidy.vcf.gz)
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX"

/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools query -f '%ID\n'1000GP.ALL.ChrX.QCed.vcf.gz |    sort | uniq -d > 1000GP_chrX.ploidy.dup_id;

    if [[ -s 1000GP_chrX.ploidy.dup_id ]]; then
        /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -e ID=@1000GP_chrX.ploidy.dup_id 1000GP.ALL.chrX.ploidy.vcf.gz   -Oz -o 1000GP.ALL.chrX.ploidy.Unique.vcf.gz
    else  mv 1000GP.ALL.chrX.ploidy.vcf.gz 1000GP.ALL.chrX.ploidy.Unique.vcf.gz
    fi

mv 1000GP.ALL.chrX.ploidy.Unique.vcf.gz 1000GP.ALL.Final.vcf.gz
/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf ${wd}1000GP.ALL.Final.vcf.gz

####If multiallelic sites are present in your data, in order to preserve them throughout the protocol, set ID field with unique IDs e.g. in format CHR_POS_REF_ALT. (RSIDs might contain duplicates, when the multiallelic sites are decomposed.) (我们的文件从下载开始就没有rsID,所以先自己转换为ID的格式e.g. in format CHR_POS_REF_ALT.）
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX"
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools annotate  --set-id '%CHROM\_%POS\_%REF\_%ALT' 1000GP.ALL.ChrX.QCed.vcf.gz     -Oz -o 1000GP.ALL.ChrX.final.vcf.gz

/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf 1000GP.ALL.ChrX.final.vcf.gz


#####Reference panel allele frequencies(这一步可以省略）
####Generate a tab-delimited file of the reference panel allele frequencies, one variant per line, with columns CHR, SNP (in format CHR_POS_REF_ALT), REF, ALT, AF (including the header line).
# Check if the VCF does NOT contain AF in the INFO field,
# and calculate it with bcftools +fill-tags plugin
#在多个来源或版本的 VCF 文件中，对于 AF 值的表示方式可能不一致。通过重新计算 AF 值，您可以确保在所有文件中 AF 值采用相同的格式，从而方便后续的数据整合和分析。(该步骤可以省略)
    /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools +fill-tags 1000GP.ALL.ChrX.final.vcf.gz -Oz -o 1000GP.ALL.ChrX.final.AF.vcf.gz -- -t AF 

# Generate a tab-delimited header(该步骤可以省略)
echo -e 'CHR\tSNP\tREF\tALT\tAF' > 1000GP_imputation_all.frq

# Query the required fields from the VCF file (该步骤可以省略)
# and append to the allele frequency file 
    /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools query  -f '%CHROM\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%INFO/AF\n'  1000GP.ALL.ChrX.final.AF.vcf.gz  >> 1000GP_imputation_all.frq

####### Convert each file to bref format (该步骤可以省略）
java -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/beagle-5.0/beagle.03Jul18.40b.jar gt=1000GP.ALL.ChrX.final.vcf.gz out=1000GP.ALL.ChrX.final.bref


###########提取1000GP的NPR区域，并分为男性和女性子集
cd /jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/


###提取NPR区域
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -i 'POS >= 2781479 && POS <= 155701382'  1000GP.ALL.ChrX.final.vcf.gz  -Oz  -o 1000GP.ALL.ChrX.NPR.final.vcf.gz

/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf 1000GP.ALL.ChrX.NPR.final.vcf.gz

###提取男性 (注意：需要转换为diploid,不然会报错)
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -S 1000GP.Male.list 1000GP.ALL.ChrX.NPR.final.vcf.gz -Oz -o 1000GP.ALL.Male.ChrX.NPR.final.vcf.gz

/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf  1000GP.ALL.Male.ChrX.NPR.final.vcf.gz

# Fix the chromosome X ploidy to phased diplo
# Requires a ploidy.txt file containing
# space-separated CHROM,FROM,TO,SEX,PLOIDY
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/06.Phasing.ChromosomeX/"
echo "chrX 1 156040895 M 2" > ${wd}ploidy.txt
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools +fixploidy ${wd}1000GP.ALL.Male.ChrX.NPR.final.vcf.gz -Ov -- -p ${wd}ploidy.txt |  sed 's#0/0#0\|0#g;s#1/1#1\|1#g' | /hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -Oz -o ${wd}1000GP.ALL.Male.ChrX.NPR.final01.vcf.gz

mv ${wd}1000GP.ALL.Male.ChrX.NPR.final01.vcf.gz ${wd}1000GP.ALL.Male.ChrX.NPR.final.vcf.gz

##提取女性
/hwfssz4/BC_PUB/Software/03.Soft_ALL/bcftools-1.9/bcftools view -S 1000GP.Female.list 1000GP.ALL.ChrX.NPR.final.vcf.gz --force-samples -Oz -o 1000GP.ALL.Female.ChrX.NPR.final.vcf.gz
/hwfssz4/BC_PUB/Software/03.Soft_ALL/htslib-2.1/tabix -f -p vcf 1000GP.ALL.Female.ChrX.NPR.final.vcf.gz



