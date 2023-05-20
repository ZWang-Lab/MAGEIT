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


sourceCpp("/home/slcxding/File/MAGEA.cpp")
source("/home/slcxding/File/ADABFGE.R")

pval_GESAT = function(y, X, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = as.matrix(X)
    p_gesat = GESAT(G,y,E,X)$pvalue
  
    return(p_gesat)
}


pval_iSKAT = function(y, X, e1, G){

    G = as.matrix(G)
    y = as.matrix(y)
    E = as.matrix(e1)
    X = as.matrix(X)
    p_iSKAT = iSKAT(G,y,E,X)$pvalue

    return(p_iSKAT)
}

pval_rare = function(y, X, e1, G){

    X_rare = cbind(e1,X)
    X_rare = as.matrix(X_rare)     
    rare_result = rareGE(y,G,X_rare)
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
result = MiSTi(data,p = n,m=1,d=k)$pvalue

return(result)
}


pval_aGE = function(y, x1, x2, e1, G){

cov = cbind(e1,x1,x2)
cov = as.matrix(cov)
result = aGE(y,G,cov,model="gaussian")
p_aGE = result[7:8]

return(p_aGE)
}


pval_ADABF = function(y, X, e1, G){
  
  cov = X
  p_adabf = ADABFGE(y, G, e1, Y.Type="C", E.Type="D", Cov=cov, Sig=0.00001, Precision.P=1)
  return(p_adabf)
  
}
















