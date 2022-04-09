# PCAUFEOPSD
PCA and TD based unsupervised FE with optimised SD

This is a sample souce code to perform analysis in the paper 
"Principal component analysis- and tensor decomposition- based unsupervised feature extraction to select more reasonable differentially methylated cytosines: Optimization of standard deviation versus state-of-art methods"
https://doi.org/10.1101/2022.04.02.486807

The following files should be downloaded and placed in the suitable directories

#./

GSE77965_series_matrix.txt.gz

https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77965/matrix/GSE77965_series_matrix.txt.gz

GPL13534-11288.txt.gz

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534

wgEncodeDukeDnaseFibrobl.fdr01peaks.hg19.bb

http://mirror.vcu.edu/vcu/encode/DNase-seq_Peaks/wgEncodeDukeDnaseFibrobl.fdr01peaks.hg19.bb

#./GSE77965

sample_table.csv 

GSM2062665_7796806124_R01C01_Grn.idat.gz
GSM2062665_7796806124_R01C01_Red.idat.gz
GSM2062666_7796806124_R02C01_Grn.idat.gz
GSM2062666_7796806124_R02C01_Red.idat.gz
GSM2062667_7796806124_R03C01_Grn.idat.gz
GSM2062667_7796806124_R03C01_Red.idat.gz
GSM2062668_7796806124_R04C01_Grn.idat.gz
GSM2062668_7796806124_R04C01_Red.idat.gz
GSM2062669_7796806124_R05C01_Grn.idat.gz
GSM2062669_7796806124_R05C01_Red.idat.gz
GSM2062670_7796806124_R06C01_Grn.idat.gz
GSM2062670_7796806124_R06C01_Red.idat.gz
GSM2062671_7796806124_R01C02_Grn.idat.gz
GSM2062671_7796806124_R01C02_Red.idat.gz
GSM2062672_7796806124_R02C02_Grn.idat.gz
GSM2062672_7796806124_R02C02_Red.idat.gz
GSM2062673_7796806124_R03C02_Grn.idat.gz
GSM2062673_7796806124_R03C02_Red.idat.gz
GSM2062674_7796806124_R04C02_Grn.idat.gz
GSM2062674_7796806124_R04C02_Red.idat.gz
GSM2062675_7796806124_R05C02_Grn.idat.gz
GSM2062675_7796806124_R05C02_Red.idat.gz
GSM2062676_7796806124_R06C02_Grn.idat.gz
GSM2062676_7796806124_R06C02_Red.idat.gz
GSM2062677_7786923078_R01C01_Grn.idat.gz
GSM2062677_7786923078_R01C01_Red.idat.gz
GSM2062678_7786923078_R01C02_Grn.idat.gz
GSM2062678_7786923078_R01C02_Red.idat.gz
GSM2062679_7786923078_R02C02_Grn.idat.gz
GSM2062679_7786923078_R02C02_Red.idat.gz
GSM2062680_7786923078_R02C01_Grn.idat.gz
GSM2062680_7786923078_R02C01_Red.idat.gz
GSM2062681_7878191101_R05C01_Grn.idat.gz
GSM2062681_7878191101_R05C01_Red.idat.gz
GSM2062682_7973201041_R04C02_Grn.idat.gz
GSM2062682_7973201041_R04C02_Red.idat.gz
GSM2062683_7786923078_R04C01_Grn.idat.gz
GSM2062683_7786923078_R04C01_Red.idat.gz
GSM2062684_7786923078_R03C01_Grn.idat.gz
GSM2062684_7786923078_R03C01_Red.idat.gz
GSM2062685_7786923078_R03C02_Grn.idat.gz
GSM2062685_7786923078_R03C02_Red.idat.gz
GSM2062686_7878191157_R03C02_Grn.idat.gz
GSM2062686_7878191157_R03C02_Red.idat.gz

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77965


#./GSE42308

sample_table.csv 

GSM1038308_5958154021_R01C01_Grn.idat.gz
GSM1038308_5958154021_R01C01_Red.idat.gz
GSM1038309_5958154021_R02C01_Grn.idat.gz
GSM1038309_5958154021_R02C01_Red.idat.gz
GSM1038310_5958154021_R03C01_Grn.idat.gz
GSM1038310_5958154021_R03C01_Red.idat.gz
GSM1038311_5958154021_R04C01_Grn.idat.gz
GSM1038311_5958154021_R04C01_Red.idat.gz
GSM1038312_5958154021_R05C01_Grn.idat.gz
GSM1038312_5958154021_R05C01_Red.idat.gz
GSM1038313_5958154021_R06C01_Grn.idat.gz
GSM1038313_5958154021_R06C01_Red.idat.gz

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42308

#./GSE34864_RAW

GSM856495_oocytes-BDF1-1.cpgs.txt.gz
GSM856496_oocytes-BDF1-2.cpgs.txt.gz
GSM856501_zygote-BDF1-BDF1-1.cpgs.txt.gz
GSM856502_zygote-BDF1-BDF1-2.cpgs.txt.gz

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE34864


