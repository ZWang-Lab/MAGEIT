
####################################################################################
#########################    Load R libraries and files   ##########################
####################################################################################

library(Rcpp)
library(RcppArmadillo)
library(CompQuadForm)
library(foreach)
library(iterators)
library(doParallel)
library(gtx)
library(MASS)
library(corpcor)

sourceCpp("MAGEIT.cpp")
source("MAGEIT.R")

####################################################################################
###########################         Data generation          #######################
####################################################################################

registerDoParallel(20)  ## Set the cores for parallel computing
set.seed(666)  ## Set the seed

ns=10000  ## Set the subject number
haplotype = read.table("haplotype.txt",header=F)  ## Read the 0-1 haplotype file
n.1 = dim(haplotype)[1]
haplotype = as.matrix(haplotype)
sample1=sample(1:n.1,ns,replace=TRUE)
sample2=sample(1:n.1,ns,replace=TRUE)
genotype=haplotype[sample1,]+haplotype[sample2,]  ## Generate the SNP data
nonraredup.common=genotype[,which(colSums(genotype)>2*ns*0.05)]  ## Choose common SNPs with rare allel frequency >0.05
nonraredup.common=as.matrix(nonraredup.common)
n.2 = dim(nonraredup.common)[2]
nonraredup.rare=genotype[,which(colSums(genotype)>2*ns*0.005 & colSums(genotype)<2*ns*0.05)]  ## Choose rare SNPs with rare allel frequency >0.005 & <0.05
nonraredup.rare=as.matrix(nonraredup.rare)
n.3 = dim(nonraredup.rare)[2]  

####################################################################################
###########################        Initialization        ###########################
####################################################################################

n = 1000  ## Iteration numbers
ss = 5000  ## Sample size
sn_common = 10  ## Common SNP numbers
sn_rare = 40  ## Rare SNP numbers
cn_common = 2  ## Causal common SNP numbers
cn_rare = 8  ## Causal rare SNP numbers
n_common.main = 2  ## Number of common SNPs which have genetic main effects
n_rare.main = 8  ## Number of rare SNPs which have genetic main effects
alpha_common_main = 0.09  ## Mean of genetic main effects for common SNPs
alpha_rare_main = 0.17  ## Mean of genetic main effects for rare SNPs
alpha_common = 0.19  ## Mean of interaction effects for common SNPs
alpha_rare = 0.59  ## Mean of interaction effects for rare SNPs
 
####################################################################################
###########################          Simulation          ###########################
####################################################################################

method.name = c("p_MAGEIT_RAN", "p_MAGEIT_FIX")
cat(c(method.name,"\n"), file="power_result.txt",append=TRUE )

oper = foreach(i = 1:n, .combine = rbind) %dopar%
  {

    x1 = rnorm(ss, mean=62.4, sd=11.5)  ## Mimicking age							
    x2 = rbinom(ss, prob=0.52, size=1)  ## Mimicking sex							
    e1 = rbinom(ss, prob=0.5, size=1)  ## Environmental factor
    indicator = sample(1:10000, ss, replace=F)

    indi_common = sample(1:n.2,sn_common,replace = F)
    none_indi_common = sample(1:sn_common,cn_common)
    G_common = rep(10,sn_common*ss)
    G_common = matrix(G_common,ncol=sn_common)
    SNP=NULL
    for (j in 1:sn_common)
    {
      SNP = nonraredup.common[,indi_common[j]]										
      G_common[,j] = SNP[indicator];
    }

    indi_rare = sample(1:n.3,sn_rare,replace = F)
    none_indi_rare = sample(1:sn_rare,cn_rare)
    G_rare = rep(10,sn_rare*ss)
    G_rare = matrix(G_rare,ncol=sn_rare)
    SNP=NULL
    for (j in 1:sn_rare)
    {
      SNP = nonraredup.rare[,indi_rare[j]]										
      G_rare[,j] = SNP[indicator];
    }

    indi.common.main = sample(1:sn_common, n_common.main, replace=F)
    G.common.main = G_common[, indi.common.main]

    indi.rare.main = sample(1:sn_rare, n_rare.main, replace=F)
    G.rare.main = G_rare[, indi.rare.main]

	G = cbind(G_common, G_rare)														
    epsilon = rnorm(ss, mean=0, sd=1.5)
    SNP_common = G_common[,none_indi_common]
    SNP_rare = G_rare[,none_indi_rare]
    SNP = cbind(SNP_common, SNP_rare)

	alpha1 = runif(n_common.main, min=alpha_common_main-0.02, max=alpha_common_main+0.02)
	alpha2 = runif(n_rare.main, min=alpha_rare_main-0.02, max=alpha_rare_main+0.02)
	beta1 = runif(cn_common, min=alpha_common-0.02, max=alpha_common+0.02)
	beta2 = runif(cn_rare, min=alpha_rare-0.02, max=alpha_rare+0.02)
	SNP.sign = rep(c(1,1),(cn_common+cn_rare))

	y_common = y_rare = y_common.main = y_rare.main = rep(0,ss)
	for (i in 1:cn_common){
		y_common = y_common + beta1[i]*SNP[,i]*SNP.sign[i]*e1
	}
	for (j in (cn_common+1):(cn_common+cn_rare)){
		y_rare = y_rare + beta2[j-cn_common]*SNP[,j]*SNP.sign[j]*e1
	}

	if (n_common.main == 1){
		y_common.main = y_common.main + alpha1*G.common.main*SNP.sign[1]
	} else {
		for (k in 1:n_common.main){
			y_common.main = y_common.main + alpha1[k]*G.common.main[,k]*SNP.sign[k]
		}
	}

	if(n_rare.main == 1){
		y_rare.main = y_rare.main + alpha2*G.rare.main*SNP.sign[10]
	} else{
		for (t in 1:n_rare.main){
			y_rare.main = y_rare.main + alpha2[t]*G.rare.main[,t]*SNP.sign[t+n_common.main]
		}
	}

	y = 0.05*x1+0.057*x2+0.64*e1+ epsilon + y_common + y_rare + y_common.main + y_rare.main

	p_MAGEIT_RAN = MAGEIT_RAN.C(y, x1, x2, e1, G, 2)
	p_MAGEIT_FIX = MAGEIT_FIX.C(y, x1, x2, e1, G)

	tt = cbind(p_MAGEIT_RAN, p_MAGEIT_FIX)
	cat(c(tt,"\n"), file="power_result.txt",append=TRUE )
}



