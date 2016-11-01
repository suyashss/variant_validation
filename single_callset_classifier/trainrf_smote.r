library("randomForest",quietly=T)
library("ROCR",quietly=T)
library("DMwR")
library("plyr",quietly=T)
args<-commandArgs(T)
calltype<-args[[1]]
featurefile<-args[[2]]
omniposfile<-args[[3]]
omninegfile<-args[[4]]
features=read.table(featurefile,header=T)
omnipos=read.table(omniposfile,header=F,stringsAsFactors=F)
omnineg=read.table(omninegfile,header=F,stringsAsFactors=F)
colnames(omnipos)<-c("CHR","START","END")
colnames(omnineg)<-c("CHR","START","END")

cat("Calltype is",calltype,"\n")
seqcalls=features[features[,c(calltype)]==1,]

pos=join(seqcalls,omnipos,type="inner")
neg=join(seqcalls,omnineg,type="inner")

trnpos=pos
trnneg=neg

cat("Number of positive examples=",dim(trnpos)[1],"\n")
cat("Number of negative examples=",dim(trnneg)[1],"\n")

trainpos=cbind(trnpos,Res=rep("2",dim(trnpos)[1]))
trainneg=cbind(trnneg,Res=rep("1",dim(trnneg)[1]))

train=rbind(trainpos,trainneg)

nmin=min(dim(trainneg)[1],dim(trainpos)[1])
nmax=max(dim(trainneg)[1],dim(trainpos)[1])

seqcallF=paste(calltype,"F",sep="")
train1=train[,c("Res","AB","Coverage","TotalCoverage","INDELPROX",seqcallF,"QUAL","Dels","FS","MQ","QD","HaplotypeScore","AC","GC","BaseQRankSum","MQRankSum","ReadPosRankSum")]
#train1=train[,c("Res","AB","Coverage","TotalCoverage","INDELPROX","IlluminaAF")]
rfformula=as.formula(paste("Res~AB+Coverage+TotalCoverage+INDELPROX+",seqcallF,"+QUAL+Dels+FS+MQ+QD+HaplotypeScore+AC+GC"))
print(rfformula)
train.imputed=SMOTE(rfformula,data=train1,perc.over=2000,perc.under=300)
cat("Computed new data frame\n")
sites=randomForest(rfformula,data=train.imputed,ntrees=1000,importance=T,na.action=na.omit)
print(importance(sites))
write.table(importance(sites),file=paste(calltype,".importance",sep=""),quote=F)
print(sites)

sites.pr=predict(sites,type="prob")[,1]

sites.pr_obj=prediction(sites.pr,train.imputed$Res[!apply(train.imputed,1,function(x){any(is.na(x))})])

sites.perf.auc=performance(sites.pr_obj,"auc")
#print(sites.perf.auc)
sites.perf=performance(sites.pr_obj,"tpr","fpr")
cutoffs<-data.frame(cut=sites.perf@alpha.values[[1]],fpr=sites.perf@x.values[[1]],tpr=sites.perf@y.values[[1]])
write.table(x=cutoffs,file=paste(calltype,"_rocinfo.txt",sep=""),quote=F,row.names=F,col.names=T)
chosen=cutoffs[cutoffs$fpr<0.05,]
cutoffval=chosen[dim(chosen)[1],1]
tprval=chosen[dim(chosen)[1],3]
print(tail(chosen))
cat("Finished\n")

library("sqldf")
notinomnipos<-sqldf('select * from seqcalls except select * from pos')
notinomniposneg<-sqldf('select * from notinomnipos except select * from neg')

notomnipred<-predict(sites,notinomniposneg,"prob")
cat("Number of sites",dim(notomnipred)[1],"\n")
cat("Positive sites",sum(notomnipred[,1]>cutoffval,na.rm=T),"\n")

seqcallspred<-predict(sites,seqcalls,"prob")
cat("Number of sites",dim(seqcalls)[1],"\n")
npospred=sum(seqcallspred[,1]>cutoffval,na.rm=T);cat("Positive sites",npospred,"\n")
cat("NA  sites",sum(is.na(seqcallspred[,1])),"\n")

resultvec=c(tprval,npospred,dim(seqcalls)[1])
write(resultvec,file=paste(calltype,"_classifer.results",sep=""))

cat(levels(seqcalls$CalledIn),"\n")
scoretable=cbind(levels(seqcalls$CalledIn)[seqcalls$CalledIn],seqcallspred[,1],seqcalls$END)
head(scoretable)
write.table(scoretable,file=paste(calltype,"_classifier.sitescores",sep=""),row.names=F,col.names=F,quote=F)
