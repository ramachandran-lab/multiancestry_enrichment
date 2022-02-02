rm(list = ls())
#set directory
setwd("~")

#if to plot unstructured results, set property == 'oo', otherwise 'wo'
property = 'oo'

#load in data
load(paste("gwas_power_", property,"_nobuffer.RData", sep=""))
load(paste("gwas_fdr_", property,"_nobuffer.RData", sep=""))
load(paste("gwas_power_se_", property,"_nobuffer.RData", sep=""))
load(paste("gwas_fdr_se_", property,"_nobuffer.RData", sep=""))

#making figure
pdf(paste("~/Desktop/manuscripts/round1/chr1/gwas_power_plot/gwas_power_", property,"_nobuffer.pdf", sep=""), width = 10, height = 4)
#barplot
barCenters = barplot(mydata_power,  beside = TRUE, xlab = 'Scenarios', ylab = 'Power', ylim = c(0,1), pch = 16, bty = 'n', col = c('orange', 'royalblue', 'darkgreen', 'tomato'))
#add axis
axis(1, at = c(1:8)*5-1.5,cex.axis = 0.8,labels = c(expression(paste('h'^2, '=0.2/n=5000', sep='')), expression(paste('h'^2, '=0.2/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.6/n=5000', sep='')), expression(paste('h'^2, '=0.6/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.2/n=5000', sep='')), expression(paste('h'^2, '=0.2/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.6/n=5000', sep='')), expression(paste('h'^2, '=0.6/n=10000', sep=''))))

#add error bar, divide by 10 is to derive standard errors as myerrors store the standard deviations, and the number of replicates is 10
arrows(barCenters, mydata_power - myerror_power/10, barCenters, mydata_power + myerror_power/10,length=0.05, angle=90, code=3, lwd = 1.5)
abline(v = 20.5, lwd = 1.5, lty = 2, col = "grey")
#add legend
legend('topright', c(expression(paste('gene-', epsilon)), 'RSS', 'SKAT', 'GWAS'), fill = c('orange', 'royalblue', 'darkgreen', 'tomato'), cex =0.8, bty = 'n')
legend(5, 1.0, c('Sparsity = 1 percent'), bty = 'n')
legend(25, 1.0, c('Sparsity = 10 percent'), bty = 'n')
dev.off()


#making figure

pdf(paste("~/Desktop/manuscripts/round1/chr1/gwas_power_plot/gwas_FDR_", property,"_nobuffer.pdf", sep=""), width = 10, height = 4)
#barplot
barCenters = barplot(mydata_fdr, beside = TRUE, xlab = 'Scenarios', ylab = 'FDR', ylim = c(0,1), pch = 16, bty = 'n', col = c('orange', 'royalblue', 'darkgreen', 'tomato'))
#add axis
axis(1, at = c(1:8)*5-1.5,cex.axis = 0.8,labels = c(expression(paste('h'^2, '=0.2/n=5000', sep='')), expression(paste('h'^2, '=0.2/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.6/n=5000', sep='')), expression(paste('h'^2, '=0.6/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.2/n=5000', sep='')), expression(paste('h'^2, '=0.2/n=10000', sep='')),
                                                    expression(paste('h'^2, '=0.6/n=5000', sep='')), expression(paste('h'^2, '=0.6/n=10000', sep=''))))
#add error bar, divide by 10 is to derive standard errors as myerrors store the standard deviations, and the number of replicates is 10
arrows(barCenters, mydata_fdr - myerror_fdr/10, barCenters, mydata_fdr + myerror_fdr/10,length=0.05, angle=90, code=3, lwd = 1.5)
abline(v = 20.5, lwd = 1.5, lty = 2, col = "grey")
#add legend
legend('topright', c(expression(paste('gene-', epsilon)), 'RSS', 'SKAT', 'GWAS'), fill = c('orange', 'royalblue', 'darkgreen', 'tomato'), cex =0.8, bty = 'n')
legend(5, 1.0, c('Sparsity = 1 percent'), bty = 'n')
legend(25, 1.0, c('Sparsity = 10 percent'), bty = 'n')
dev.off()
