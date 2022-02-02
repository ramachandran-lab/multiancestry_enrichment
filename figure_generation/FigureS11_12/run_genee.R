library(genee)
args <- commandArgs(TRUE)
chr=1
num_sam=2000
p_genes=as.numeric(args[1])
if(p_genes == 0.01){
  p_snps = 0.00125
}else{
  p_snps = 0.03
}
ncausal_intergenic=5
pve=as.numeric(args[2])
repl=as.numeric(args[3])
load(paste("/gpfs/data/sramacha/ukbiobank_jun17/multi_ancestry/2000/raw_",chr,"_",num_sam,".RData", sep=""))
load(paste("/gpfs/data/sramacha/ukbiobank_jun17/multi_ancestry/2000/PCs_",chr,"_",num_sam,".RData", sep=""))
X=raw[[1]]
print(dim(X))
bim=raw[[2]]
map=raw[[3]]
#get the information about snps in each gene
glist.hg19 <- read.delim("/gpfs/data/sramacha/WTCCC/glist-hg19", header=FALSE)
all_genes<-as.character(glist.hg19[,4])##all possible genes
genes_num_snps<-rep(0,length = length(all_genes))#number of snps in each possible gene, could be 0
genes_snps<-list()#only contain genes which have snps and store the snps in those genes
j=0
my_genes<-c()#genes have snps
for (i in 1:length(all_genes)) {
  temp_gene<-all_genes[i]
  if(as.character(glist.hg19[i,1])!="X" && as.character(glist.hg19[i,1])!="Y")
  {
    temp_chr=as.numeric(as.character(glist.hg19[i,1]))
    start=glist.hg19[i,2]-50000
    end=glist.hg19[i,3]+50000
    genes_num_snps[i]=length(which(map[,1]==temp_chr & start<=map[,4] & end>=map[,4]))
    if(genes_num_snps[i]!=0 && genes_num_snps[i]!=1){
      j=j+1
      genes_snps[[j]]=which(map[,1]==temp_chr & start<=map[,4] & end>=map[,4])
      my_genes<-c(my_genes,temp_gene)
    }
  }
}

unique_snps<-function(gene, my_genes, genes_snps){
  ngenes=length(genes_snps)
  gene_index=match(gene, my_genes)
  gene.id=c(1:ngenes)
  gene.id_rest=gene.id[-gene_index]
  dif_snps<-list()##different snps just in this gene not in other genes
  for(i in gene.id){
    dif_snps[[i]]=setdiff(genes_snps[[gene_index]],genes_snps[[i]])
  }
  dif_snps=dif_snps[-gene_index]
  result=Reduce(intersect, dif_snps)
  return(result)
}

genes_snps_unique<-lapply(my_genes, my_genes=my_genes, genes_snps=genes_snps, unique_snps)
ngenes_snps_unique<-sapply(genes_snps_unique, length)


####random assign causal genes and causal snps only in genes with unique snps
assign_random_unique<-function(genes_snps, genes_snps_unique, ngenes_snps_unique, my_genes, p_genes, p_snps, map, ncausal_intergenic){
  nsnp=dim(map)[1]
  ngenes=length(my_genes)
  ncausal_genes=round(ngenes*p_genes+0.5)#number of causal genes
  qualified_genes.id=which(ngenes_snps_unique!=0)
  stopifnot(length(qualified_genes.id)>ncausal_genes)###assert that qualified genes are more than the number of causal genes
  causal_genes=sample(qualified_genes.id, ncausal_genes)#id of causal genes
  genes.id=c(1:ngenes)
  noncausal_genes=genes.id[-causal_genes]
  ncausal_snps_in_genes=round(nsnp*p_snps+0.5)#number of causal snps in genes
  snps_in_genes<-c()#all snps in genes
  for(i in 1:length(my_genes)){
    snps_in_genes<-c(snps_in_genes, genes_snps[[i]])
  }
  unique_snps_in_causal_genes<-c()#all unique snps in causal genes
  for(i in 1:length(causal_genes)){
    unique_snps_in_causal_genes<-c(unique_snps_in_causal_genes, genes_snps_unique[[causal_genes[i]]])
  }
  snps.id=c(1:nsnp)#id of all snps
  snps_in_genes=sort(unique(snps_in_genes))
  unique_snps_in_causal_genes=sort(unique(unique_snps_in_causal_genes))
  stopifnot(length(unique_snps_in_causal_genes)>ncausal_snps_in_genes)###assert that qualified genes are more than causal genes
  snps_intergenic=snps.id[-snps_in_genes]#id of intergenic snps
  causal_snps_in_genes<-sample(unique_snps_in_causal_genes, ncausal_snps_in_genes)#id of causal snps in genes
  stopifnot(length(snps_intergenic)>ncausal_intergenic)###assert that the number intergenic snps are more than causal snps intergenic region
  causal_snps_intergenic<-sample(snps_intergenic, ncausal_intergenic)#id of causal snps in intergenic
  causal_snps=sort(unique(c(causal_snps_in_genes, causal_snps_intergenic)))
  noncausal_snps=snps.id[-causal_snps]
  ncausal_gene=rep(0,length(my_genes))
  for (i in 1:length(my_genes)) {
    ncausal_gene[i]=length(intersect(genes_snps[[i]], causal_snps))
  }
  causal_genes=which(ncausal_gene!=0)
  noncausal_genes=which(ncausal_gene==0)
  all_genes=list()
  all_genes$causal_genes=causal_genes
  all_genes$noncausal_genes=noncausal_genes
  all_snps=list()
  all_snps$causal_snps=causal_snps
  all_snps$causal_snps_in_genes=causal_snps_in_genes
  all_snps$causal_snps_intergenic=causal_snps_intergenic
  all_snps$noncausal_snps=noncausal_snps
  assign_result=list()
  assign_result$all_genes=all_genes
  assign_result$all_snps=all_snps
  return(assign_result)
}


########################################################
########################################################
###############phenotype simulation#####################

pheno_simu<-function(pve, X, assign_result, PCs){
  num_snps=dim(X)[2]
  num_inds=dim(X)[1]
  Xmarginal=X[,assign_result$all_snps$causal_snps]
  Xnoise=X[,assign_result$all_snps$noncausal_snps]
  #marginal effects
  beta_marginal=rnorm(dim(Xmarginal)[2], 0, sd=1)
  y_marginal=Xmarginal%*%beta_marginal
  beta_marginal=beta_marginal*as.numeric(sqrt(pve/var(y_marginal)))
  y_marginal=Xmarginal%*%beta_marginal
  #error effects
  y_err=rnorm(num_inds)
  y_err=y_err*as.numeric(sqrt((1-pve)/var(y_err)))
  #
  y=y_marginal+y_err
  true_betas=vector(length=num_snps)
  true_betas[assign_result$all_snps$causal_snps]=beta_marginal
  pheno_result=list()
  pheno_result$beta_marginal=beta_marginal
  pheno_result$true_betas=true_betas
  pheno_result$y=y
  pheno_result$y_err=y_err
  pheno_result$y_marginal=y_marginal
  ###pcs
  ### Define the effects of the PCs ###
  beta_pc=rnorm(dim(PCs)[2])
  y_pcs=PCs%*%beta_pc
  beta_pc=beta_pc*as.numeric(sqrt(0.1/var(y_pcs)))
  y_pcs=PCs%*%beta_pc
  #error effects
  y_err_wPCs=rnorm(num_inds)
  y_err_wPCs=y_err_wPCs*as.numeric(sqrt((1-pve-0.1)/var(y_err_wPCs)))
  y_wPCs=y_marginal+y_err_wPCs+y_pcs
  pheno_result$wPCs$beta_marginal_wPCs=beta_marginal
  pheno_result$wPCs$true_betas_wPCs=true_betas
  pheno_result$wPCs$y_wPCs=y_wPCs
  pheno_result$wPCs$y_pcs=y_pcs
  pheno_result$wPCs$y_err_wPCs=y_err
  pheno_result$wPCs$y_marginal_wPCs=y_marginal
  return(pheno_result)
}



########################################################
########################################################
###############Linear Regression #######################
summary_stats<-function(X, y, PCs){
  num_snps=dim(X)[2]
  pval=vector(length = num_snps)
  beta_est=vector(length = num_snps)
  beta_se=vector(length = num_snps)
  for (i in 1:num_snps) {
    relation=lm(y~X[,i])
    beta_est[i]=summary(relation)$coefficient[2,1]
    beta_se[i]=summary(relation)$coefficient[2,2]
    pval[i]=summary(relation)$coefficient[2,4]
  }
  pval_wPCs=vector(length = num_snps)
  beta_est_wPCs=vector(length = num_snps)
  beta_se_wPCs=vector(length = num_snps)
  for (i in 1:num_snps) {
    relation_wPCs=lm(y~X[,i]+PCs)
    beta_est_wPCs[i]=summary(relation_wPCs)$coefficient[2,1]
    beta_se_wPCs[i]=summary(relation_wPCs)$coefficient[2,2]
    pval_wPCs[i]=summary(relation_wPCs)$coefficient[2,4]
  }
  gwas_result=list()
  gwas_result$beta_est=beta_est
  gwas_result$beta_se=beta_se
  gwas_result$pval=pval
  gwas_result$beta_est_wPCs=beta_est_wPCs
  gwas_result$beta_se_wPCs=beta_se_wPCs
  gwas_result$pval_wPCs=pval_wPCs
  return(gwas_result)
}

########################################################
###############My simulation#######################
###############My simulation#######################
my_simulator<-function(X, p_snps, p_genes, ncausal_intergenic, glist.hg19,map, bim, my_genes, genes_snps, genes_snps_unique, ngenes_snps_unique, pve, PCs){
  num_inds=dim(X)[1]
  num_snps=dim(X)[2]
  assign_result=assign_random_unique(genes_snps, genes_snps_unique, ngenes_snps_unique, my_genes, p_genes, p_snps, map, ncausal_intergenic)
  pheno_result=pheno_simu(pve, X, assign_result, PCs)
  gwas_result=summary_stats(X, pheno_result$y, PCs)
  gwas_result_wPCs_pheno=summary_stats(X, pheno_result$wPCs$y_wPCs, PCs)
  ########################################################
  ###############Summary Statistics#######################
  
  ##Shrinkage
  #ash_result<-ash(gwas_result$beta_est, gwas_result$beta_se)$result
  #shrink_result=list()
  #shrink_result$ash_betas=ash_result[,9]
  #shrink_result$ash_se=ash_result[,10]
  simu_result=list()
  simu_result$assign_result=assign_result
  simu_result$pheno_result=pheno_result
  simu_result$gwas_result=gwas_result
  simu_result$gwas_result_wPCs_pheno=gwas_result_wPCs_pheno
  #simu_result$shrink_result=shrink_result
  ######
  mydata = bim[,c(1,2,4)]
  mydata = cbind(mydata, gwas_result$beta_est)
  mydata_wPCs = bim[,c(1,2,4)]
  mydata_wPCs = cbind(mydata_wPCs, gwas_result_wPCs_pheno$beta_est_wPCs)
  simu_result$mydata = mydata
  simu_result$mydata_wPCs = mydata_wPCs
  return(simu_result)
}

simu_result=my_simulator(X, p_snps, p_genes, ncausal_intergenic, glist.hg19,map, bim, my_genes, genes_snps, genes_snps_unique, ngenes_snps_unique, pve, PCs)




setwd("/gpfs/data/sramacha/ukbiobank_jun17/multi_ancestry/genee/")

load('/gpfs/data/sramacha/ukbiobank_jun17/multi_ancestry/2000/cor_1_2000.RData')
ld = as.matrix(ld)
fit = genee(simu_result$mydata, ld, alpha = 0.5)
fit_wPCs = genee(simu_result$mydata_wPCs, ld, alpha = 0.5)
all_data = list()
all_data$simu_result = simu_result
all_data$fit = fit
all_data$fit_wPCs = fit_wPCs
#
save(all_data, file = paste("/gpfs/data/sramacha/ukbiobank_jun17/multi_ancestry/2000/res/res_",pve, "_", p_genes, "_",repl, ".RData" ))




