
MAGE_RAN.body = function(y.star, G.star, S.star, M, n.fix)
{
  
	S.kernel = GetLinearKernel(S.star)
	G.kernel = GetLinearKernel(G.star)
	S1 = GetS_est(G.kernel,S.kernel,M, n.fix)
	S1.inv = mat_inv_3(S1)
	q = GetQ_est(G.kernel,S.kernel,M,y.star)
	delta = ComputeDelta(S1,q)

	sigma.hat = delta[2]

	S1_null_num = c(S1[1,1],S1[1,3],S1[3,1],S1[3,3])
	S1_null = matrix(S1_null_num, byrow=T, ncol=2)
	q_null = q[c(1,3)]

	delta_null = ComputeDelta(S1_null,q_null)

	omega.hat = delta_null[1]
	tau.hat = delta_null[2]

	H = S1.inv[2,1]*G.kernel + S1.inv[2,2]*S.kernel + S1.inv[2,3]*diag(dim(M)[1])
	GM = omega.hat*G.kernel + tau.hat*M
	HM = MatMult(GM,H)

	lambda = comeigen(HM)

	davies_val = davies(sigma.hat,lambda,lim=500000, acc=1e-8)
	p_mage = davies_val$Qq
  
  return(p_mage) 
}


MAGE_FIX.body = function(y.star, S.star, M)
{
  
	S.kernel = GetLinearKernel(S.star)
	S1 = GetS_est_Gfix(S.kernel, M)
	S1.inv = MatInv(S1)
	q = GetQ_est_Gfix(S.kernel,M,y.star)
	delta = ComputeDelta(S1,q)

	sigma.hat = delta[1]

	S1_null = S1[2,2]
	q_null = q[2]

	tau.hat = 1/S1_null*q_null

	H = S1.inv[1,1]*S.kernel + S1.inv[1,2]*diag(dim(M)[1])

	GM = tau.hat*M
	HM = MatMult(GM,H)

	lambda = comeigen(GM,H)

	davies_val = davies(sigma.hat,lambda,lim=500000, acc=1e-8)
	p_mage = davies_val$Qq
  
  return(p_mage) 
}



















trans.ytype = function(y, n_mcmc, k)
{

    ### Save the labels ###
    cc = y
    
    ### Find the Number of Cases and Controls ###
    n.cases = sum(cc==1)
    n.controls = sum(cc==0)
    
    ### Set the Threshold ###
    thresh=qnorm(1-k,mean=0,sd=1)
    
    ### Bernoulli Distributed Case and Control Data ###
    Pheno=rep(NA,length(cc));
    Pheno[cc==0] = etruncnorm(b=thresh)
    Pheno[cc==1] = etruncnorm(a=thresh)

  continuous.y = Pheno 
  return(continuous.y)
}

MAGE.continous.common = function(y, x1, x2, e1, G)
{
  y = y - mean(y)
  S = diag(e1)%*%G
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G)
  p.continuous.common = MAGE.body(M,y,S)
  return(p.continuous.common)
}

MAGE_RAN.C = function(y, x1, x2, e1, G, n.cov)
{
	n.fix = n.cov+2
	sub.n = dim(G)[1]  
	maf = colSums(G)/(sub.n*2)
	weight.num = dbeta(maf,1,25)
	weight = diag(weight.num)
	M = ComputeProj(x1, x2, e1)
	G.weight = MatMult(G,weight)
	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	G.star = MatMult(M,G.weight)
	S.star = MatMult(M,S)

  p.continuous.ran = MAGE_RAN.body(y.star, G.star, S.star, M, n.fix)
  return(p.continuous.ran)
}


MAGE_FIX.C = function(y, x1, x2, e1, G)
{
	sub.n = dim(G)[1]  
	maf = colSums(G)/(sub.n*2)
	weight.num = dbeta(maf,1,25)
	weight = diag(weight.num)
	M = ComputeProj_Gfix(x1, x2, e1, G)
	G.weight = MatMult(G,weight)
	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

  p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
  return(p.continuous.fix)
}


MAGE_BUR.C = function(y, x1, x2, e1, G)
{
	sub.n = dim(G)[1]  
	maf = colSums(G)/(sub.n*2)
	weight.num = dbeta(maf,1,25)
	weight = diag(weight.num)
	M = ComputeProj_Gfix(x1, x2, e1, G%*%weight.num)
	G.weight = MatMult(G,weight)
	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

  p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
  return(p.continuous.fix)
}









MAGE.binary.common = function(y, x1, x2, e1, G)
{
  y = trans.ytype(y, n_mcmc=10000, k=0.2)
  y = y - mean(y)
  S = diag(e1)%*%G
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.binary.common = MAGE.body(M,y,S)
  return(p.binary.common)
}

MAGE.binary.rare = function(y, x1, x2, e1, G, maf)
{
  weight.num = dbeta(maf,1,25)
  weight = diag(weight.num)
  y = trans.ytype(y, n_mcmc=10000, k=0.2)
  y = y - mean(y)
  S = diag(e1)%*%G%*%weight
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.binary.rare = MAGE.body(M,y,S)
  return(p.binary.rare)
}

MAGE.binary.weight = function(y, x1, x2, e1, G, maf, pre)
{
  weight.num = dbeta(maf,1,25)
  weight = diag(weight.num)
  y = trans.ytype(y, n_mcmc=10000, k=pre)
  y = y - mean(y)
  S = diag(e1)%*%G%*%weight
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.binary.weight = MAGE.body(M,y,S)
  return(p.binary.weight)
}

MAGE.continuous.weight = function(y, x1, x2, e1, G, maf)
{
  weight.num = dbeta(maf,1,25)
  weight = diag(weight.num)
  y = y - mean(y)
  S = diag(e1)%*%G%*%weight
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.continous.weight = MAGE.body(M,y,S)
  return(p.continous.weight)
}

MAGE.binary.original = function(y, x1, x2, e1, G, pre)
{
  y = trans.ytype(y, n_mcmc=10000, k=pre)
  y = y - mean(y)
  S = diag(e1)%*%G
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.binary.original = MAGE.body(M,y,S)
  return(p.binary.original)
}


MAGE.continuous.original = function(y, x1, x2, e1, G)
{
  y = y - mean(y)
  S = diag(e1)%*%G
  G = Standard(G)
  S = Standard(S)
  M = ComputeProj_Gfix(x1, x2, e1, G) 
  p.continous.original = MAGE.body(M,y,S)
  return(p.continous.original)
}

MAGE.continuous = function(y, x1, x2, e1, G)
{
  sub.n = dim(G)[1]  
  maf = colSums(G)/(sub.n*2)
  common.index = maf>=0.05
  rare.index = maf<0.05
  
  G.common = G[,common.index]
  G.rare = G[,rare.index]
  maf.rare = maf[rare.index]
  
  p.continuous.common = MAGE.continous.common(y, x1, x2, e1, G.common)
  p.continuous.rare = MAGE.continous.rare(y, x1, x2, e1, G.rare, maf.rare)
  ACAT.T = tan((0.5-p.continuous.common)*pi)+tan((0.5-p.continuous.rare)*pi)
  p.acat = 0.5 - atan(ACAT.T/2)/pi
  p.min = min(p.continuous.common, p.continuous.rare)
  p.weight = MAGE.continuous.weight(y, x1, x2, e1, G, maf)
  p.original = MAGE.continuous.original(y, x1, x2, e1, G)
  pvalue = c(p.acat, p.min, p.weight, p.original)
  names(pvalue) = c("p.acat", "p.min", "p.weight", "p.original")
  return(pvalue)
}


MAGE.binary = function(y, x1, x2, e1, G, pre)
{
  sub.n = dim(G)[1]  
  maf = colSums(G)/(sub.n*2)
  #common.index = maf>=0.05
  #rare.index = maf<0.05
  
  #G.common = G[,common.index]
  #G.rare = G[,rare.index]
  #maf.rare = maf[rare.index]
  
  #p.binary.common = MAGE.binary.common(y, x1, x2, e1, G.common, pre)
  #p.binary.rare = MAGE.binary.rare(y, x1, x2, e1, G.rare, maf.rare, pre)
  #ACAT.T = tan((0.5-p.binary.common)*pi)+tan((0.5-p.binary.rare)*pi)
  #p.acat = 0.5 - atan(ACAT.T/2)/pi
  #p.min = min(p.binary.common, p.binary.rare)
  p.weight = MAGE.binary.weight(y, x1, x2, e1, G, maf, pre)
  p.original = MAGE.binary.original(y, x1, x2, e1, G, pre)
  pvalue = c(p.weight, p.original)
  names(pvalue) = c("p.weight", "p.original")
  return(pvalue)
}

pval_MAGE = function(y, x1, x2, e1, G, y.type, pre)
{
  if(y.type == "C")
  {
    pvalue = MAGE.continuous(y, x1, x2, e1, G)
  }
  else if(y.type == "D")
  {
    pvalue = MAGE.binary(y, x1, x2, e1, G, pre)
  }
  else
  {
    print("have to specity the outcome type")
  }
  return(pvalue)
}













