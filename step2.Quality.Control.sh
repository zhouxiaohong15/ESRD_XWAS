wd="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/01.Merge.Raw.data/"
wd2="/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/tf11_data/ESRD_ASA/PhD.paper/02.ALL.QCed.hg.37.data/"



/jdfssz1/ST_HEALTH/P20Z10200N0041/Zhouxiaohong/conda_env/R4.0/lib/genetics.binaRies/bin/plink  --bfile ${wd}All.Merge.ASA --ci 0.95 --geno 0.05 --hwe 0.000001 --maf 0.01 --make-bed --mind 0.05 --not-chr Y,26,0  --missing --allow-no-sex  --out ${wd2}All.QCed.ASA
