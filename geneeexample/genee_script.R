args = commandArgs(trailingOnly=TRUE)
chr = args[2]
trait = args[3]
library("genee")

load(paste("chr",chr,".RData",sep = ""))
load(paste(trait,"_",chr,"_betas.RData",sep = ""))
mydata[,3]=mydata[,2]

myresult=genee(mydata,ld,alpha=0.5,upper=50000,lower=50000)
write.table(myresult,paste(trait,"_",chr,"_genee_result.txt",sep = ""),sep = '\t')