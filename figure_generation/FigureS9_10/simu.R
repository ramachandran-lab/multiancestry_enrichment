
args <- commandArgs(TRUE)
p_snps=as.numeric(args[1])
p_genes=as.numeric(args[2])
pve=as.numeric(args[3])
p_share=as.numeric(args[4])
repl = as.integer(args[5])
buffer = 50000

X_afr = as.matrix(read.table("african_sub.txt"))
PC_afr = as.matrix(read.table("PCs_african_sub.txt"))
X_bri = as.matrix(read.table("british_sub.txt"))
PC_bri = as.matrix(read.table("PCs_british_sub.txt"))
bim = read.table("african_sub.bim")

#get the information about snps in each gene
glist.hg19 <- read.delim("/users/wcheng8/data/WTCCC/glist-hg19", header=FALSE)
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
    start=glist.hg19[i,2]-buffer
    end=glist.hg19[i,3]+buffer
    genes_num_snps[i]=length(which(bim[,1]==temp_chr & start<=bim[,4] & end>=bim[,4]))
    if(genes_num_snps[i]!=0 && genes_num_snps[i]!=1){
      j=j+1
      genes_snps[[j]]=which(bim[,1]==temp_chr & start<=bim[,4] & end>=bim[,4])
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
assign_random_unique<-function(genes_snps, genes_snps_unique, ngenes_snps_unique, my_genes, p_genes, p_snps, bim, p_share){
  nsnp=dim(bim)[1]
  ngenes=length(my_genes)
  ncausal_genes=round(ngenes*p_genes+0.5)#number of causal genes
  nsnps_per_gene = unlist(lapply(genes_snps, length))
  qualified_genes.id = which(nsnps_per_gene>30)
  # qualified_genes.id=which(ngenes_snps_unique!=0)
  # stopifnot(length(qualified_genes.id)>ncausal_genes)###assert that qualified genes are more than the number of causal genes
  causal_genes=sample(qualified_genes.id, ncausal_genes)#id of causal genes
  genes.id=c(1:ngenes)
  noncausal_genes=genes.id[-causal_genes]
  #
  causal_snps_A = c()
  causal_snps_B = c()
  for (g in causal_genes) {
    #ncausal_in_gene = round(p_snps*length(genes_snps[[g]])+0.5)
    ncausal_in_gene = round(p_snps*length(genes_snps[[g]])-0.5)
    n_total = round((2.0-p_share)*ncausal_in_gene)
    #n_total = round((2.0-p_share)*ncausal_in_gene + 0.5)
    temp_causal_all = sample(genes_snps[[g]], n_total)
    L = length(temp_causal_all)
    to_be_share = round(L*p_share)
    if(p_share == 0){
      causal_snps_A = c(causal_snps_A, temp_causal_all[1:to_be_share])
      causal_snps_B = c(causal_snps_B, temp_causal_all[(to_be_share+1):L]) 
    }else if(p_share == 0.25 || p_share == 0.5){
      L2 = round((L - to_be_share)/2)
      causal_snps_A = c(causal_snps_A, temp_causal_all[1:to_be_share], temp_causal_all[(to_be_share+1):(L2 + to_be_share)])
      causal_snps_B = c(causal_snps_B, temp_causal_all[1:to_be_share], temp_causal_all[(L2 + to_be_share + 1): L])
    }else if(p_share == 1.0){
      causal_snps_A = c(causal_snps_A, temp_causal_all)
      causal_snps_B = c(causal_snps_B, temp_causal_all)
    }
  }
  causal_snps_A = unique(sort(causal_snps_A))
  causal_snps_B = unique(sort(causal_snps_B))
  noncausal_snps_A = c(1:nsnp)[-causal_snps_A]
  noncausal_snps_B = c(1:nsnp)[-causal_snps_B]
  all_genes=list()
  all_genes$causal_genes=causal_genes
  all_genes$noncausal_genes=noncausal_genes
  all_snps=list()
  res_A = list()
  res_B = list()
  res_A$causal_snps=causal_snps_A
  res_A$noncausal_snps=noncausal_snps_A
  res_B$causal_snps=causal_snps_B
  res_B$noncausal_snps=noncausal_snps_B
  all_snps[[1]] = res_A
  all_snps[[2]] = res_B
  assign_result=list()
  assign_result$all_genes=all_genes
  assign_result$all_snps=all_snps
  return(assign_result)
}


########################################################
########################################################
###############phenotype simulation#####################

pheno_simu<-function(pve, X, assign_result, PCs, idx){
  num_snps=dim(X)[2]
  num_inds=dim(X)[1]
  Xmarginal=X[,assign_result$all_snps[[idx]]$causal_snps]
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
my_simulator<-function(X_afr, X_bri, p_snps, p_genes, glist.hg19,bim, my_genes, genes_snps, genes_snps_unique, ngenes_snps_unique, pve, PC_afr, PC_bri, p_share){
  assign_result=assign_random_unique(genes_snps, genes_snps_unique, ngenes_snps_unique, my_genes, p_genes, p_snps, bim, p_share)
  pheno_result_afr=pheno_simu(pve, X_afr, assign_result, PC_afr, 1)
  pheno_result_bri=pheno_simu(pve, X_bri, assign_result, PC_bri, 2)
  # gwas_result_afr=summary_stats(X_afr, pheno_result_afr$y, PC_afr)
  gwas_result_wPCs_pheno_afr=summary_stats(X_afr, pheno_result_afr$wPCs$y_wPCs, PC_afr)
  # gwas_result_bri=summary_stats(X_bri, pheno_result_bri$y, PC_bri)
  gwas_result_wPCs_pheno_bri=summary_stats(X_bri, pheno_result_bri$wPCs$y_wPCs, PC_bri)
  
  ########################################################
  ###############genee#######################
  mydata_afr = bim[,c(1,2,4)]
  mydata_afr = cbind(mydata_afr, gwas_result_wPCs_pheno_afr$beta_est_wPCs)
  mydata_bri = bim[,c(1,2,4)]
  mydata_bri = cbind(mydata_bri, gwas_result_wPCs_pheno_bri$beta_est_wPCs)
  simu_result=list()
  simu_result$assign_result=assign_result
  simu_result$pheno_result_afr=pheno_result_afr
  simu_result$pheno_result_bri=pheno_result_bri
  simu_result$gwas_result_wPCs_pheno_afr=gwas_result_wPCs_pheno_afr
  simu_result$gwas_result_wPCs_pheno_bri=gwas_result_wPCs_pheno_bri
  simu_result$mydata_afr = mydata_afr
  simu_result$mydata_bri = mydata_bri
  return(simu_result)
}

simu_result=my_simulator(X_afr, X_bri, p_snps, p_genes, glist.hg19,bim, my_genes, genes_snps, genes_snps_unique, ngenes_snps_unique, pve, PC_afr, PC_bri, p_share) 
save(simu_result, file = paste("/users/wcheng8/data/wcheng8/simulation/multi_ancestry/simu_data/simu_", pve, "_", p_snps, "_", p_genes, "_", p_share,"_repl",repl,".RData", sep=""))

#
# for(i in 1:50){
#   set.seed(round(runif(1, min=1, max=10^8)))
#   simu_result=my_simulator(X_afr, X_bri, p_snps, p_genes, glist.hg19,bim, my_genes, genes_snps, genes_snps_unique, ngenes_snps_unique, pve, PC_afr, PC_bri, p_share) 
#   save(simu_result, file = paste("/users/wcheng8/data/wcheng8/simulation/multi_ancestry/simu_data/simu_", pve, "_", p_snps, "_", p_genes, "_", p_share,"_repl",repl,".RData", sep=""))
# }



