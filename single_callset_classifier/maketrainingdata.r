library("plyr")

args<-commandArgs(T)
bedfile<-args[[1]]
allelebalancefile<-args[[2]]
coveragefile<-args[[3]]
indelproxfile<-args[[4]]
sitefeaturefile<-args[[5]]
illuminafreqfile<-args[[6]]
rtgfreqfile<-args[[7]]
gatkfreqfile<-args[[8]]
outfilewithNA<-args[[9]]
outfilenoNA<-args[[10]]


overlaps=read.table(bedfile)
colnames(overlaps)=c("CHR","START","END","NCalledIn","CalledIn","illumina","rtg","gatk")

ab=read.table(allelebalancefile,header=F)
colnames(ab)=c("CHR","START","END","AB")

coverage=read.table(coveragefile,header=F)
colnames(coverage)=c("CHR","START","END","TotalCoverage")

indelprox=read.table(indelproxfile,header=F)
colnames(indelprox)=c("CHR","START","END","INDELPROX")

sitefeatures=read.table(sitefeaturefile,header=T)

illuminaf=read.table(illuminafreqfile,header=F,skip=1)
colnames(illuminaf)=c("ib37CHR","END","illuminaNALLELES","illuminaNCHR","illuminaREF","illuminaF")
rtgf=read.table(rtgfreqfile,header=F,skip=1)
colnames(rtgf)=c("rb37CHR","END","rtgNALLELES","rtgNCHR","rtgREF","rtgF")
gatkf=read.table(gatkfreqfile,header=F,skip=1)
colnames(gatkf)=c("gb37CHR","END","gatkNALLELES","gatkNCHR","gatkREF","gatkF")

alldata=join_all(list(overlaps,sitefeatures,ab,coverage,indelprox,illuminaf,rtgf,gatkf))

selectcolumns=c("CHR","START","END","NCalledIn","CalledIn","illumina","rtg","gatk","QUAL","AC","Dels","FS","HaplotypeScore","MQ","QD","BaseQRankSum","MQRankSum","ReadPosRankSum","GC","AB","TotalCoverage","INDELPROX","illuminaF","rtgF","gatkF","Coverage")

savedata=alldata[,selectcolumns]
write.table(savedata,file=outfilewithNA,quote=F,row.names=F)

savedata[is.na(savedata)]=0
write.table(savedata,file=outfilenoNA,quote=F,row.names=F)
