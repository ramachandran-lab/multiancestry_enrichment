args = commandArgs(trailingOnly=TRUE)
trait = args[1]
library("UpSetR")
mydata=read.csv(paste("../gene_manhattan_plot/",trait,".upset.input.csv",sep=''))
sets = dim(mydata)[2]-1
pdf(paste(trait,'.upset.pdf',sep=''),onefile = F)
upset(mydata,nsets=sets, order.by = "freq", point.size = 3.5,mainbar.y.label = "Gene Intersections", sets.x.label = "Significant Genes",text.scale=c(1))
dev.off()