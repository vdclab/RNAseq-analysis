library(NOISeq)
library(edgeR)
library(statmod)
library(org.Hs.eg.db)
library(VennDiagram)
library('GO.db')
library(xlsx)

#10 oct 2016 and 15 nov 2021
#author: Cécile Pereira and Pedro Salguero
#Analysis Leticia data

#sourcing data
mycounts=read.table("counts/table_reads_counts.csv",sep="\t",header=TRUE) 
mycounts2=read.table("counts/table_reads_counts_dropbox.csv",sep="\t",header=TRUE)

normalize_dropbox=read.table("counts/table_after_normalization_dropbox.csv",sep=" ")

#Compare Dropbox and my own read-count-file
mycounts_df <- mycounts[order(row.names(mycounts)),]
mycounts2_df <- mycounts2[order(row.names(mycounts2)), colnames(mycounts_df)]

write.xlsx(mycounts_df,'counts/table_before_normalization.xlsx')
write.xlsx(mycounts2_df,'counts/table_before_normalization_dropbox.xlsx')

mycounts <- mycounts_df
mycounts2 <- mycounts2_df

comparation <- rowSums(mycounts_df != mycounts2_df)
comparation <- comparation[comparation!=0]

setwd("./results/")
mycounts['tgt',] 

#control: number of reads for the gene tgt
#             vc_sv_m1lb_rz_ec_S20_L002_ vc_sv_t1ni_rz_ec_S18_L002_ vc_sv_m2lb_rz_ec_S21_L002_ vc_sv_m3lb_rz_ec_S12_L002_ vc_sv_t3lb_rz_ec_S15_L002_
# tgt                        973                        261                        776                       2966                        182
# vc_sv_t2lb_rz_ec_S14_L002_ vc_sv_t1lb_rz_ec_S13_L002_ vc_sv_m3ni_rz_ec_S17_L002_ vc_sv_t3ni_rz_ec_S19_L002_ vc_sv_m2ni_rz_ec_S22_L002_
# tgt                        160                        261                        781                        381                        801
# vc_sv_t2ni_rz_ec_S23_L002_ vc_sv_m1ni_rz_ec_S16_L002_
# tgt                        259                       1029

#4 conditions: 
#t-lb
#t-ni
#m-lb
#m-ni

lst_to_delete <- c("__ambiguous", "__no_feature", "__alignment_not_unique","__not_aligned", "__too_low_aQual")
otherFeatures <- mycounts[lst_to_delete,]
mycounts <- mycounts[!rownames(mycounts) %in% lst_to_delete,]
#filter for low counts (see page 73)
A=rowSums(mycounts)
dataFilt = mycounts[A>10,]

#groups
tgt=grep('_t[1-3][a-z]*_',colnames(mycounts))
m=grep('_m[1-3][a-z]*_',colnames(mycounts))
lb=grep('_[a-z][1-3]lb_',colnames(mycounts))
ni=grep('_[a-z][1-3]ni_',colnames(mycounts))

#groups
# tgt=grep('_t[1-3][a-z]*_',colnames(mycounts2))
# m=grep('_m[1-3][a-z]*_',colnames(mycounts2))
# lb=grep('_[a-z][1-3]lb_',colnames(mycounts2))
# ni=grep('_[a-z][1-3]ni_',colnames(mycounts2))

groups=c()
del=c()
mil=c()
for(i in 1:12){
  g=''
  if(i %in% tgt){
    g='tgt_'
    del=c(del,'tgt')
  }
  else{
    g='m_'
    del=c(del,'m')
  }
  if(i %in% lb){
    mil=c(mil,'lb')
    g=paste(g,'lb',sep='')
  }
  else{
    mil=c(mil,'ni')
    g=paste(g,'ni',sep='')
  }
  groups=c(groups,g)
}
df=data.frame(mil,del,groups)

gfffile=read.table('../genome/GCF_000005845.2_ASM584v2_genomic.gff',sep='\t',quote="",stringsAsFactors = FALSE)
gffgenes=gfffile[gfffile[,3]=='gene',]
length=gffgenes[,5]-gffgenes[,4]#length of the gene (end-beginning position)
tmpname=gffgenes[,9]
#extraction of the name of the gene
namelength=c()
for(t in tmpname){
  tmp=strsplit(as.vector(t),';')
  tmp=tmp[[1]][3]
  tmp=gsub("Name=",'',tmp)
  #print(tmp)
  namelength=c(namelength,tmp)
}
names(length)=namelength

mydata=readData(dataFilt, factors = df, length=length)
mydata

# ExpressionSet (storageMode: lockedEnvironment)
# assayData: 4386 features, 12 samples 
# element names: exprs 
# protocolData: none
# phenoData
# sampleNames: vc_sv_t1ni_rz_ec_S18_L002_ vc_sv_m2ni_rz_ec_S22_L002_ ... vc_sv_t1lb_rz_ec_S13_L002_ (12 total)
# varLabels: mil del groups
# varMetadata: labelDescription
# featureData
# featureNames: aaaD aaaE ... zwf (4386 total)
# fvarLabels: Length
# fvarMetadata: labelDescription
# experimentData: use 'experimentData(object)'
# Annotation: 

myPCA=dat(mydata,type='PCA')
pdf("PCA_beforeARSyN.pdf", width = 5, height = 5)
#par(mfrow=c(1,2))
explo.plot(myPCA,factor='groups')
#problem with m+lb?
dev.off()

mydata2corr2 = NOISeq::ARSyNseq(mydata, factor = "groups", batch=FALSE, norm='rpkm', logtransf = FALSE) #removing the noize from the data

#head(assayData(mydata2corr2)$exprs)
normalize_result <- assayData(mydata2corr2)$exprs
write.table(normalize_result,'table_after_normalization.csv')
write.xlsx(normalize_result,'table_after_normalization.xlsx')

normalize_dropbox <- normalize_dropbox[order(row.names(normalize_dropbox)),colnames(normalize_result)]
write.xlsx(normalize_dropbox,'table_after_normalization_dropbox.xlsx')

sub_normalize_result <- normalize_result[rownames(normalize_dropbox),]
write.table(sub_normalize_result,'subtable_after_normalization.csv')
write.xlsx(sub_normalize_result,'subtable_after_normalization.xlsx')

comparation_after <- rowSums(normalize_result != normalize_dropbox)
comparation_after <- comparation_after[comparation_after==0]

#PCA
myPCA=dat(mydata2corr2,type='PCA')
pdf("PCA_afterARSyN.pdf", width = 5, height = 5)
#par(mfrow=c(1,2))
explo.plot(myPCA,factor='groups')
#explo.plot(myPCA,factor='TimeTreat', samples = c(1,3))
dev.off()

#a) voom + limma
gr=factor(df$groups)
design=model.matrix(~0+gr)
colnames(design)=levels(gr)
vo=voom(mydata2corr2,design,plot=TRUE)
# This  converts  the  counts  to  log-counts  per  million  with  associated precision weights
#The limma-voom method assumes that rows with zero or very low counts have been removed.

fit=lmFit(vo,design)
mycontrasts=makeContrasts(
  mediaWT=m_ni-m_lb,
  mediaTGT=tgt_ni-tgt_lb,
  WT_TGT_LB=tgt_lb-m_lb,
  WT_TGT_NI=tgt_ni-m_ni,
  levels=design
)
fit3=contrasts.fit(fit,mycontrasts)
fit3=eBayes(fit3)
result=decideTests(fit3,method="global")
summary(result)

### NEW
#           mediaWT mediaTGT WT_TGT_LB WT_TGT_NI
# Down      1202     1025       434       644
# NotSig    2279     2571      3407      2938
# Up         905      790       545       804

write.table(result,'results_DE.csv')
write.xlsx(result,'results_DE.xlsx')

for(d in colnames(mycontrasts)){
  filename=paste("toptreat_",d,".csv",sep='')
  write.table(topTreat(fit3,coef=d,number=Inf),filename)
  filename=paste("toptreat_",d,".xlsx",sep='')
  write.xlsx(topTreat(fit3,coef=d,number=Inf),filename)
  pdf(paste("volcanoplot_",d,".pdf",sep=''))
  volcanoplot(fit3,coef=d,names=fit3$genes$NAME,main=d)
  dev.off()
}
#counts (before normalization), table after the normalization, log2fold change
complete <- as.data.frame(fit3)
rownames(complete) <- rownames(fit3)

write.table(complete,'complete_fit_table.csv')
write.xlsx(complete,'complete_fit_table.xlsx')

plotMD(fit3)
abline(0,0,col="blue")

qqt(fit3$t,df=fit3$df.prior+fit3$df.residual,pch=16,cex=0.2)
abline(1,1)

#features deleted
write.table(as.data.frame(otherFeatures), file = 'htseq-count_features.csv')
write.xlsx(as.data.frame(otherFeatures), file = 'htseq-count_features.xlsx')
