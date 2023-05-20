library(Rcpp,lib.loc="/home/slcxding/R_libs")
library(RcppArmadillo,lib.loc="/home/slcxding/R_libs")
library(CompQuadForm,lib.loc="/home/slcxding/R_libs")
library(rareGE,lib.loc="/home/slcxding/R_libs")
library(iSKAT)
library(foreach,lib.loc="/home/slcxding/R_libs")
library(iterators,lib.loc="/home/slcxding/R_libs")
library(doParallel,lib.loc="/home/slcxding/R_libs")
library(gtx,lib.loc="/home/slcxding/R_libs")
library(MASS,lib.loc="/home/slcxding/R_libs")
library(corpcor,lib.loc="/home/slcxding/R_libs")
library(MiSTi,lib.loc="/home/slcxding/R_libs")
library(aGE,lib.loc="/home/slcxding/R_libs")
library(truncnorm,lib.loc="/home/slcxding/R_libs")

sourceCpp("/home/slcxding/File/MAGEA.cpp")
source("/home/slcxding/File/ADABFGE.R")


 
pval_MAGE = function(y, X, e1, G, pre){

sub.n = dim(G)[1]  
maf = colSums(G)/(sub.n*2)
weight.num = dbeta(maf,1,25)
weight = diag(weight.num)

L = matrix(rep(NA,ss*n_mcmc),ncol=n_mcmc)
for(t in 1:n_mcmc)
{
  ### Save the labels ###
  cc = y
  
  ### Find the Number of Cases and Controls ###
  n.cases = sum(cc==1)
  n.controls = sum(cc==0)
  
  ### Set the Threshold ###
  thresh=qnorm(1-pre,mean=0,sd=1)
  
  ### Bernoulli Distributed Case and Control Data ###
  Pheno=rep(NA,length(cc));
  Pheno[cc==0] = rtruncnorm(n.controls,b=thresh)
  Pheno[cc==1] = rtruncnorm(n.cases,a=thresh)
  L[,t] = Pheno
}

y = rowMeans(L)

    y = y - mean(y)
    S = diag(e1)%*%G%*%weight
    G = Standard(G)
    S = Standard(S)
    M = ComputeProj_Gfixx(X, e1, G)
    y.star = M%*%y
    S.star = M%*%S  
    S.star = GetLinearKernel(S.star)
    
    S1 = GetS_est_Gfix(S.star,M)
    q = GetQ_est_Gfix(S.star,M,y.star)
    delta = ComputeDelta(S1,q)
    lambda = comeigen_fixed_Gfix(S1, q, S.star, M) 
    S_null = S1[2,2]
    q_null = q[2]    
    delta_null = 1/S_null*q_null
    V = delta_null*M
    V_diag = diag(V)
    eigen_value = lambda*V_diag
    davies_val = davies(delta[1],eigen_value,lim=500000, acc=1e-8)
    p_mage = davies_val$Qq

    return(p_mage)
}

pval_GESAT = function(y, X, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = as.matrix(X)
    p_gesat = GESAT(G,y,E,X, out_type="D")$pvalue
  
    return(p_gesat)
}


pval_iSKAT = function(y, X, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = as.matrix(X)
    p_iSKAT = iSKAT(G,y,E,X, out_type="D",MAF_cutoff=1)$pvalue

    return(p_iSKAT)
}

pval_rare = function(y, X, e1, G){

    X_rare = cbind(e1,X)
    X_rare = as.matrix(X_rare)     
    rare_result = rareGE(y,G,X_rare,family = "binomial")
    pval_intfix = rare_result$pINT_FIX
    pval_intran = rare_result$pINT_RAN
    p = cbind(pval_intfix, pval_intran)

    return(p)

}

pval_MiSTi = function(y, X, e1, G){

n = ncol(G)
k = ncol(X)
data = cbind(y,e1,X,G)
data = as.data.frame(data)
result = MiSTi(data,p = n,m=1,d=k,outcome_type = "Binary")$pvalue

return(result)
}


pval_aGE = function(y, X, e1, G){

cov = cbind(e1,X)
cov = as.matrix(cov)
result = aGE(y,G,cov,model="gaussian")
p_aGE = result[7:8]

return(p_aGE)
}


pval_ADABF = function(y, X, e1, G){
  
  cov = X
  p_adabf = ADABFGE(y, G, e1, Y.Type="D", E.Type="D", Cov=cov, Sig=0.00001, Precision.P=1)
  return(p_adabf)
  
}


