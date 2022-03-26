require(rtracklayer)
bb <- import("wgEncodeDukeDnaseFibrobl.fdr01peaks.hg19.bb")
y <- read.csv("GPL13534-11288.txt.gz",sep="\t",comment.char="#")
#GSE77965
x1 <- read.delim("GSE77965_series_matrix.txt.gz", comment.char="!")
x1p <- scale(x1[,2:13])
x1p[is.na(x1p)]<-0
PCA <- prcomp(x1p)
th <- function(sd){
    P2 <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=10000,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.1,th)$par
P <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
#FALSE   TRUE 
#397912  87665 
DMC <- x1[p.adjust(P,"BH")<0.01,1]

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[,1],y[,1])]!="",p.adjust(P,"BH")<0.01))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),p.adjust(P,"BH")<0.01))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),p.adjust(P,"BH")<0.01))

#COHCAP
require(COHCAP)
rownames(x1p) <- x1[,1]
write.table(file="GSE77965.txt",x1p,quote=F,sep="\t")
#modify first row as follows with adding "SiteID"
SiteID	GSM2062665	GSM2062666	GSM2062667	GSM2062668	GSM2062669	GSM2062670	GSM2062671	GSM2062672	GSM2062673	GSM2062674	GSM2062675	GSM2062676

beta.file = file.path("./","GSE77965.txt")
sample.file = file.path("./","GSE77965_sample.txt")
beta.table = COHCAP.annotate(beta.file, "GSE77965", "./",
                             platform="450k-UCSC")
filtered.sites = COHCAP.site(sample.file,beta.table,"GSE77965",
                             "./", ref="N")
genes <- read.csv("CpG_Site/GSE77965_CpG_site_filter.txt",sep="\t")
genes <- unique(unlist(strsplit(genes[,4],";")))

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[,1],y[,1])]!="",x1[,1] %in% genes[,1]))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),x1[, 1] %in% genes[, 1]))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),x1[, 1] %in% genes[, 1]))

#ChAMP
#excute in ./GSE77965
require(ChAMP)
myLoad <- champ.load("./",arraytype="450K")
champ.QC() 
myNorm <- champ.norm()
myDMP <- champ.DMP()
table(myDMP[[1]]$adj.P.Val<0.01)
#FALSE  TRUE 
#19486  8155 

DMC <- unique(myDMP[[1]]$gene[myDMP[[1]]$adj.P.Val<0.01])

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[, 1], y[, 1])] != "",x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))



#GSE42308
x1 <- read.csv("GSE42308_series_matrix.txt.gz", comment.char="!",sep="",skip=60)
x1[,-1] <- scale(x1[,-1])
x1[,-1][is.na(x1[,-1])]<-0
PCA <- prcomp(x1[,-1])
th <- function(sd){
    P2 <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=10000,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.5,th,method="Brent",lower=0,upper=1)$par
P <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)
#FALSE   TRUE 
#388486  97091 
DMC <- x1[p.adjust(P,"BH")<0.01,1]

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[,1],y[,1])]!="",p.adjust(P,"BH")<0.01))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),p.adjust(P,"BH")<0.01))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),p.adjust(P,"BH")<0.01))

#COHCAP
require(COHCAP)
x1p <- scale(x1[,-1])
rownames(x1p) <- x1[,1]
write.table(file="GSE42308.txt",x1p,quote=F,sep="\t")
#modify first row as follows with adding "SiteID"
SiteID	GSM1038308	GSM1038309	GSM1038310	GSM1038311	GSM1038312	GSM1038313

beta.file = file.path("./","GSE42308.txt")
sample.file = file.path("./","GSE42308_sample.txt")
beta.table = COHCAP.annotate(beta.file, "GSE42308", "./",
                             platform="450k-UCSC")
filtered.sites = COHCAP.site(sample.file,beta.table,"GSE42308",
                             "./", ref="N")
genes <- read.csv("CpG_Site/GSE42308_CpG_site_filter.txt",sep="\t")

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[,1],y[,1])]!="",x1[,1] %in% genes[,1]))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),x1[, 1] %in% genes[, 1]))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),x1[, 1] %in% genes[, 1]))

#ChAMP
#excute in ./GSE42308
require(ChAMP)

myLoad <- champ.load("./",arraytype="450K")
champ.QC() 
myNorm <- champ.norm()
myDMP <- champ.DMP()
table(myDMP[[1]]$adj.P.Val<0.01)
#FALSE  TRUE 
#18268 82915 

#coincidence with functional sites
fisher.test(table(y$DMR[match(x1[, 1], y[, 1])] != "",x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))
fisher.test(table(!is.na(y$Enhancer[match(x1[,1],y[,1])]),x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))
fisher.test(table(!is.na(y$DHS[match(x1[,1],y[,1])]),x1[,1] %in% rownames(myDMP[[1]])[myDMP[[1]]$adj.P.Val<0.01]))


#GSE34864
zygote_2 <- read.delim("GSE34864_RAW/GSM856502_zygote-BDF1-BDF1-2.cpgs.txt.gz", header=FALSE)
zygote_1 <- read.delim("GSE34864_RAW/GSM856501_zygote-BDF1-BDF1-1.cpgs.txt.gz", header=FALSE)
oocytes_2 <- read.delim("GSE34864_RAW/GSM856496_oocytes-BDF1-2.cpgs.txt.gz", header=FALSE)
oocytes_1 <- read.delim("GSE34864_RAW/GSM856495_oocytes-BDF1-1.cpgs.txt.gz", header=FALSE)

zygote_1 <- data.frame(paste(zygote_1[,1],zygote_1[,2],sep="_"),zygote_1[,4]/zygote_1[,3])
zygote_2 <- data.frame(paste(zygote_2[,1],zygote_2[,2],sep="_"),zygote_2[,4]/zygote_2[,3])
oocytes_1 <- data.frame(paste(oocytes_1[,1],oocytes_1[,2],sep="_"),oocytes_1[,4]/oocytes_1[,3])
oocytes_2 <- data.frame(paste(oocytes_2[,1],oocytes_2[,2],sep="_"),oocytes_2[,4]/oocytes_2[,3])
colnames(zygote_1)<-c("ID","value")
colnames(zygote_2)<-c("ID","value")
colnames(oocytes_1)<-c("ID","value")
colnames(oocytes_2)<-c("ID","value")
ID <- sort(unique(c(zygote_1[,1],zygote_2[,1],oocytes_1[,1],oocytes_2[,1])))
x <- data.frame(ID=ID,o1=oocytes_1[match(ID,oocytes_1[,1]),2],o2=oocytes_2[match(ID,oocytes_2[,1]),2],z1=zygote_1[match(ID,zygote_1[,1]),2],z2=zygote_2[match(ID,zygote_2[,1]),2])
x[,-1][is.na(x[,-1])]<-0
x<- x[rowSums(x[,-1]==0)<4,] #this was not done for #TTEST_old
PCA <- prcomp(scale(x[,-1]))

th <- function(sd){
    P2 <- pchisq((PCA$x[,3]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=1000,plot=F)
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.4,th)$par 
P2 <- pchisq((PCA$x[,3]/sd)^2,1,lower.tail=F)

table(p.adjust(P2,"BH")<0.01)
#FALSE    TRUE 
#677704 197208 
DMC <- x[p.adjust(P2,"BH")<0.01,1]

#coincidence with DHS
TTEST <- rep(list(NA),23)
ad <- c(0,cumsum(bb@seqnames@lengths))
a0 <- x[p.adjust(P2,"BH")<0.01,1]
CHRS <- unlist(bb@seqnames@values)
names(TTEST) <-CHRS
for (i in c(1:length(CHRS)))
{
    cat(i," ")
    CHR <- CHRS[i]
    CHR <- paste(CHR,"_",sep="")
    a <- a0[grep(CHR,a0)]
    if (length(a)>0)
    {
        a <- gsub(CHR,"",a)
        b <- subsetByOverlaps(IRanges(start=as.numeric(a),width=1),bb@ranges[(ad[i]+1):ad[i+1]])
        index <- x[grep(CHR,x[,1]),1] %in% paste(CHR,b@start,sep="")
        PP <- P2[grep(CHR,x[,1])]
        TTEST[[i]] <- t.test(PP[index],PP[!index])
    }
}

#DSS
zygote_2 <- read.delim("GSE34864_RAW/GSM856502_zygote-BDF1-BDF1-2.cpgs.txt.gz", header=FALSE)
zygote_1 <- read.delim("GSE34864_RAW/GSM856501_zygote-BDF1-BDF1-1.cpgs.txt.gz", header=FALSE)
oocytes_2 <- read.delim("GSE34864_RAW/GSM856496_oocytes-BDF1-2.cpgs.txt.gz", header=FALSE)
oocytes_1 <- read.delim("GSE34864_RAW/GSM856495_oocytes-BDF1-1.cpgs.txt.gz", header=FALSE)

require(DSS)
colnames(zygote_1) <- c("chr","pos","N","X")
colnames(zygote_2) <- c("chr","pos","N","X")
colnames(oocytes_1) <- c("chr","pos","N","X")
colnames(oocytes_2) <- c("chr","pos","N","X")
BSobj <- makeBSseqData( list(zygote_1, zygote_2, oocytes_1, oocytes_2), c("C1","C2", "N1", "N2") )

dmlTest = DMLtest(BSobj, group1=c("C1", "C2"), group2=c("N1", "N2"))
table(dmlTest$fdr<0.01)
#FALSE    TRUE 
#1062323   57560 

#coincidence with DHS
a0 <- paste(dmlTest$chr,dmlTest$pos,sep="_")[dmlTest$fdr<0.01]
P2 <- dmlTest$pval
x <-  matrix(paste(dmlTest$chr,dmlTest$pos,sep="_"),ncol=1)
CHRS <- unlist(bb@seqnames@values)
names(TTEST) <-CHRS
for (i in c(1:length(CHRS)))
{
    cat(i," ")
    CHR <- CHRS[i]
    CHR <- paste(CHR,"_",sep="")
    a <- a0[grep(CHR,a0)]
    if (length(a)>0)
    {
        a <- gsub(CHR,"",a)
        b <- subsetByOverlaps(IRanges(start=as.numeric(a),width=1),bb@ranges[(ad[i]+1):ad[i+1]])
        index <- x[grep(CHR,x[,1]),1] %in% paste(CHR,b@start,sep="")
        PP <- P2[grep(CHR,x[,1])]
        TTEST[[i]] <- t.test(PP[index],PP[!index])
    }
}

#DMRcate
BSobj <- renameSeqlevels(BSobj, mapSeqlevels(seqlevels(BSobj), "UCSC"))

tissue <- factor(c("C","C","N","N"))

design <- model.matrix(~tissue)

colnames(design) <- gsub("tissue", "", colnames(design)) 
colnames(design)[1] <- "Intercept"
rownames(design) <- colnames(BSobj)

methdesign <- edgeR::modelMatrixMeth(design)
methdesign

cont.mat <- limma::makeContrasts(C_vs_N=Intercept-N,
                                 levels=methdesign)
cont.mat
seq_annot <- sequencing.annotate(BSobj, methdesign, all.cov = TRUE,
                                 contrasts = TRUE, cont.matrix = cont.mat,
                                 coef = "C_vs_N", fdr=0.05)

table(seq_annot@ranges@elementMetadata@listData$is.sig)
#FALSE 
#492683 

#metilene
chrs <- unique(unlist(lapply(strsplit(x[,1],"_"),"[",1)))
chrs <- chrs[order(as.numeric(gsub("chr","",chrs)))]
x_all <-NULL
for (i in c(1:length(chrs))){
    cat(i," ")
    index <- grep(paste(chrs[i],"_",sep=""),x[,1])
    pos <- as.numeric(unlist(lapply(strsplit(x[index,1],"_"),"[",2)))
    x_all <- rbind(x_all,x[index,][order(pos),])
}
xx <- strsplit(x_all[,1],"_")
xx <- data.frame(unlist(lapply(xx,"[",1)),unlist(lapply(xx,"[",2)))
x_all <- data.frame(xx,x_all[,-1])
colnames(x_all)<- c("chr","pos","g1_A","g1_B","g2_A","g2_B")
write.table(file="GSE34864_metilene.csv",x_all,row.names=F,sep="\t",quote=F)

#----  excute in linux ---
#./metilene_v0.2-8/metilene -a g1 -b g2 -f 3 GSE34864_metilene.csv > output.dat #DMC
#./metilene_v0.2-8/metilene -a g1 -b g2  GSE34864_metilene.csv > output.dat #DMR
#---- excute in linux -----
output <- read.delim("output.dat", header=FALSE)
table(output[,4]<0.05)

#DMC
#FALSE 
#528636 

#DMR
#FALSE  TRUE 
#5179    27 

#edmr
pos <- data.frame(granges(bis_1072))[,1:2]
require(emdr)
x_Diff <- data.frame(chr=pos[,1],start=as.numeric(pos[,2]),end=as.numeric(pos[,2]),
                     strand="+",pvalue=P2,qvalue=p.adjust(P2,"BH"),meth.diff=PCA$x[,2]/sd)
mydmr =edmr(x_Diff, ACF=F, DMR.methdiff=1,DMC.methdiff=1)
mysigdmr=filter.dmr(mydmr,DMR.qvalue=0.01,mean.meth.diff=1)
#GRanges object with 12662 ranges and 5 metadata columns:
    
    
#EH1072
require(bsseq)
library(ExperimentHub)
eh <- ExperimentHub()
bis_1072 <- eh[["EH1072"]]
meth <- getMeth(bis_1072)
PCA <- prcomp(scale(meth))
th <- function(sd){
    P2 <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
    hc<- hist(1-P2,breaks=1000,plot=F) 
    return(sd(hc$count[1:sum(hc$breaks<1-min(P2[p.adjust(P2,"BH")>0.01]))]))
}
sd <- optim(0.25,th)$par
P2 <- pchisq((PCA$x[,2]/sd)^2,1,lower.tail=F)
table(p.adjust(P2,"BH")<0.01)
#FALSE     TRUE 
#19680799  2186751 
aa <- paste("chr",rep(bis_1072@rowRanges@seqnames@values,bis_1072@rowRanges@seqnames@lengths,sep=""),"_",bis_1072@rowRanges@ranges@start,sep="")
DMC <- aa[p.adjust(P2,"BH")<0.01]

#coincidence with DHS
aa <- paste("chr",rep(bis_1072@rowRanges@seqnames@values,bis_1072@rowRanges@seqnames@lengths,sep=""),"_",bis_1072@rowRanges@ranges@start,sep="")
a0 <- aa[p.adjust(P2,"BH")<0.01]
TTEST <- rep(list(NA),23)
ad <- c(0,cumsum(bb@seqnames@lengths))
CHRS <- unlist(bb@seqnames@values)
names(TTEST) <-CHRS
for (i in c(1:length(CHRS)))
{
    cat(i," ")
    CHR <- CHRS[i]
    CHR <- paste(CHR,"_",sep="")
    a <- a0[grep(CHR,a0)]
    if (length(a)>0)
    {
        a <- gsub(CHR,"",a)
        b <- subsetByOverlaps(IRanges(start=as.numeric(a),width=1),bb@ranges[(ad[i]+1):ad[i+1]])
        index <- aa[grep(CHR,aa)] %in% paste(CHR,b@start,sep="")
        if (sum(index)>3)
        {
            PP <- P2[grep(CHR,aa)]
            TTEST[[i]] <- t.test(PP[index],PP[!index])
        }
    }
}

#DSS
dmlTest = DMLtest(bis_1072, group1=c("Lymph-N-Tcon-R1","Lymph-N-Tcon-R2","Lymph-N-Tcon-R3"), group2=c("Lymph-N-Treg-R1","Lymph-N-Treg-R2","Lymph-N-Treg-R3"))

#DMRcate
require(SummarizedExperiment)
bsseq::pData(bis_1072) <- data.frame(replicate=gsub(".*-", "", colnames(bis_1072)),
                                     tissue=substr(colnames(bis_1072), 1,
                                                   nchar(colnames(bis_1072))-3),
                                     row.names=colnames(bis_1072))
colData(bis_1072)$tissue <- gsub("-", "_", colData(bis_1072)$tissue)
as.data.frame(colData(bis_1072))
bis_1072 <- renameSeqlevels(bis_1072, mapSeqlevels(seqlevels(bis_1072), "UCSC"))

tissue <- factor(pData(bis_1072)$tissue)
tissue <- relevel(tissue, "Liver_Treg")

design <- model.matrix(~tissue)

colnames(design) <- gsub("tissue", "", colnames(design))
colnames(design)[1] <- "Intercept"
rownames(design) <- colnames(bis_1072)

methdesign <- edgeR::modelMatrixMeth(design)
methdesign

cont.mat <- limma::makeContrasts(treg_vs_tcon=Lymph_N_Treg-Lymph_N_Tcon,
                                 fat_vs_ln=Fat_Treg-Lymph_N_Treg,
                                 skin_vs_ln=Skin_Treg-Lymph_N_Treg,
                                 fat_vs_skin=Fat_Treg-Skin_Treg,
                                 levels=methdesign)
cont.mat
seq_annot <- sequencing.annotate(bis_1072, methdesign, all.cov = TRUE,
                                 contrasts = TRUE, cont.matrix = cont.mat,
                                 coef = "treg_vs_tcon", fdr=0.05)

table(seq_annot@ranges@elementMetadata@listData$is.sig)
#FALSE     TRUE 
#18246866     5812 

#coincidence with DHS
aa <- paste(rep(seq_annot@ranges@seqnames@values,seq_annot@ranges@seqnames@lengths,sep=""),"_",seq_annot@ranges@ranges@start,sep="")
a0 <- aa[seq_annot@ranges@elementMetadata@listData$is.sig]
P2<- seq_annot@ranges@elementMetadata@listData$ind.fdr
TTEST <- rep(list(NA),23)
ad <- c(0,cumsum(bb@seqnames@lengths))
CHRS <- unlist(bb@seqnames@values)
names(TTEST) <-CHRS
for (i in c(1:length(CHRS)))
{
    cat(i," ")
    CHR <- CHRS[i]
    CHR <- paste(CHR,"_",sep="")
    a <- a0[grep(CHR,a0)]
    if (length(a)>0)
    {
        a <- gsub(CHR,"",a)
        b <- subsetByOverlaps(IRanges(start=as.numeric(a),width=1),bb@ranges[(ad[i]+1):ad[i+1]])
        index <- aa[grep(CHR,aa)] %in% paste(CHR,b@start,sep="")
        if (sum(index)>3)
        {
            PP <- P2[grep(CHR,aa)]
            TTEST[[i]] <- t.test(PP[index],PP[!index])
        }
    }
}

#metilene

pos <- data.frame(granges(bis_1072))[,1:2]
pos <- data.frame(paste("chr",pos[,1],sep=""),pos[,2],meth[,10:15])
colnames(pos)<- c("chr","pos","g1_A","g1_B","g1_C","g2_A","g2_B","g2_C")
write.table(file="bis_1072_metilene.csv",pos,row.names=F,sep="\t",quote=F)

#----- execute in linux ---
./metilene_v0.2-8/metilene -a g1 -b g2 -f 3 -t 10 bis_1072_metilene.csv > output.dat #DMC
./metilene_v0.2-8/metilene -a g1 -b g2 -t 10  bis_1072_metilene.csv > output.dat #DMR
#----- execute in linux ---
output <- read.delim("output.dat", header=FALSE)
table(output[,4]<0.05)
#DMC
#FALSE 
#514490 

#DMR
#FALSE  TRUE 
#8678 12478 


#edmr
pos <- data.frame(granges(bis_1072))[,1:2]
require(emdr)
x_Diff <- data.frame(chr=pos[,1],start=as.numeric(pos[,2]),end=as.numeric(pos[,2]),
                     strand="+",pvalue=P2,qvalue=p.adjust(P2,"BH"),meth.diff=PCA$x[,2]/sd)
mydmr =edmr(x_Diff, ACF=F, DMR.methdiff=1,DMC.methdiff=1)
mysigdmr=filter.dmr(mydmr,DMR.qvalue=0.01,mean.meth.diff=1)
#GRanges object with 12662 ranges and 5 metadata columns:
    
    
#TD for EH1072
Z <- array(NA,c(dim(meth)[1],5,3))
Z[,1,] <- data.matrix(meth[,1:3])
Z[,2,] <- data.matrix(meth[,4:6])
Z[,3,] <- data.matrix(meth[,7:9])
Z[,4,] <- data.matrix(meth[,10:12])
Z[,5,] <- data.matrix(meth[,13:15])
Z <- apply(Z,2:3,scale)
HOSVD <- hosvd(as.tensor(Z),c(10,5,3))
HOSVD$Z@data[,2,1]

#[1]   -10.798544 -2949.303662     5.825865     5.244385   -16.020706
#[6]    44.489423    -7.521644    31.130885   -22.490978    11.688110
