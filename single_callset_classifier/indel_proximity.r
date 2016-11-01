args<-commandArgs(T)
bedfile=args[[1]]
indelfile=args[[2]]
outfile=args[[3]]

bed=read.table(bedfile,header=F)
indels=scan(indelfile)
ob=apply(bed,1,function(x){0+(min(abs(as.numeric(x[3])-indels))<10)})
res=cbind(bed,ob)
write.table(x=res,file=outfile,quote=F,sep="\t",row.names=F,col.names=F)
