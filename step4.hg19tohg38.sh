########坐标转换 , 注意：经常出现转换率不到10%，原因是flip出了问题，再就是我转换vcf的时候也出了问题 ，不是用的plink2 --ref-from-fa。
wd1="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg38.data/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/"

wd3="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/04.Split.Chrom.QCed.hg37.data/"
wd4="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg37.data/"



#############Step2. hg37 转为hg38, 从https://zhuanlan.zhihu.com/p/383252096上找到下载hg38.fa.gz的方式，并通过GATK建立hg38.dict
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/"
#建立字典
#java -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/picard-2.18.16/picard.jar CreateSequenceDictionary  R=${wd}Important.files/hg38.fa  O=${wd}Important.files/hg38.dict

#转换坐标，参考https://hemtools.readthedocs.io/en/latest/content/Bioinformatics_tools/liftovervcf.html；
for chr in {1..23}
do echo " java -jar -Xmx100g /hwfssz4/BC_PUB/Software/03.Soft_ALL/picard-2.18.16/picard.jar LiftoverVcf I=${wd3}All.QCed.hg37.Chr${chr}.vcf    O=${wd1}All.QCed.hg38.Chr${chr}.vcf     CHAIN=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/Important.files/b37ToHg38.over.chain  REJECT=rejected_variants.${chr}.vcf R=/ldfssz1/ST_BIGDATA/USER/st_bigdata/Sentieon/reference_bigdatacompute/hg38_noalt_withrandom/hg38.fa" > ${wd}05.Batch.Hg19toHg38/Hg37toHg38.Chr${chr}.sh
done;

########注意：！！！！ 在跑2  7 8 号染色体时 总是报错：ERROR	2023-08-22 11:00:54	LiftoverVcf	Encountered a contig, chr8_KI270821v1_alt that is not part of the target reference；；解决办法，参考https://gatk.broadinstitute.org/hc/en-us/community/posts/9389503444379-Issues-with-LiftoverVcf的Mohammad Deen Hayatu的方法，就是原因就是chain和hg38不配对，应该要使用gatk指定的hg38.fasta文件，而不是志那个hg38路径。官方指定的hg38.fa位于路劲/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Important.files里面

for chr in 2 7 8
do echo " java -Xmx100g -jar /hwfssz4/BC_PUB/Software/03.Soft_ALL/picard-2.18.16/picard.jar LiftoverVcf I=${wd3}All.QCed.hg37.Chr${chr}.vcf    O=${wd1}All.QCed.hg38.Chr${chr}.vcf     CHAIN=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg38.data/Important.files/b37ToHg38.over.chain  REJECT=rejected_variants.${chr}.vcf R=/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/Important.files/hg38.fa"> ${wd}05.Batch.Hg19toHg38/Hg37toHg38.Chr${chr}.sh
done;



########Step3.批量提交作业
wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/"
rm ${wd}05.Batch.Hg19toHg38/Batch.Hg37toHg38.qsub.sh;

for chr in {1..23};
do echo  "qsub -cwd -l vf=5G,p=1 -q st.q -P P18Z10200N0124  -binding linear:8 ${wd}05.Batch.Hg19toHg38/Hg37toHg38.Chr${chr}.sh " >> ${wd}05.Batch.Hg19toHg38/Batch.Hg37toHg38.qsub.sh;
done;


sh ${wd}05.Batch.Hg19toHg38/Batch.Hg37toHg38.qsub.sh
