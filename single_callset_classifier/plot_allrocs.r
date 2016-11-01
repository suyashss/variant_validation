args<-commandArgs(T)
outfile<-args[[1]]
calltypevec=c("Illumina","RTG","GATK")

df=c()
for(calltype in calltypevec){
	f=read.table(paste(tolower(calltype),"_rocinfo.txt",sep=""),header=T)
	df=rbind(df,cbind(method=rep(calltype,dim(f)[1]),f[,2:3] ) )
}
colnames(df)<-c("Method","fpr","tpr")

library("ggplot2")

p<-ggplot(df,aes(x=fpr,y=tpr,colour=Method))+geom_line(size=1.5)
p<-p+theme(axis.title.x = element_text(face="bold", size=20),
           axis.text.x  = element_text(size=16,angle=0,vjust=0.5,colour='black'))
p<-p+theme(axis.title.y = element_text(face="bold", size=20),
           axis.text.y  = element_text(size=16,colour='black'))
p<-p+ theme(legend.title = element_text(size=16),
        legend.text = element_text(size = 16))
p<-p+xlab("False Positive Rate")+ylab("True Positive Rate")
png(outfile)
print(p)
dev.off()
