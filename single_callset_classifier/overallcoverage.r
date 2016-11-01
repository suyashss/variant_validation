args<-commandArgs(T)
outfile<-args[[1]]

flist=list.files("./coverage_stats",full.names=T)
cat("Read files",flist[[1]],"\n")
overallcoverage=read.table(flist[[1]])
for(i in 2:length(flist)){
	temp=read.table(flist[[i]])
	overallcoverage[,4]=overallcoverage[,4]+temp[,4]
	if(i%%70==0){cat(i,"\n")}
}
write.table(overallcoverage,file=outfile,quote=F,row.names=F,col.names=F)
