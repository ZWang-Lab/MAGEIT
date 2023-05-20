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
	lambda = Re(comeigen(HM))

	davies_val = davies(sigma.hat,lambda,lim=500000, acc=1e-8)
	if (davies_val$ifault == 0){
		p_mage = davies_val$Qq
	} else{
		p_mage = liu(sigma.hat, lambda)$Qq
		print("mage_liu")
	} 
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

	lambda = Re(comeigen(HM))

	davies_val = davies(sigma.hat,lambda,lim=500000, acc=1e-8)
	p_mage = davies_val$Qq
	if (davies_val$ifault == 0){
		p_mage = davies_val$Qq
	} else{
		p_mage = liu(sigma.hat, lambda)$Qq
		print("mage_liu")
	}
	return(p_mage) 
}





trans.ytype = function(y, none_matrix)
{
	myprobit = glm(y~none_matrix, family=binomial(link='probit'))
	coeff = myprobit$coefficients

	mean_none = coeff[1]
	for (i in 2:length(coeff)){
		mean_none = mean_none + coeff[i]*none_matrix[,i-1]
	}
	Pheno=rep(0,length(y))
	
	Pheno[y==0] = etruncnorm(mean=mean_none, sd=1, b=0)[y==0]
	Pheno[y==1] = etruncnorm(mean=mean_none, sd=1, a=0)[y==1]  

	return(Pheno)
}



MAGE_RAN.C = function(y, x1, x2, e1, G, n.cov)
{
	n.fix = n.cov+2
	sub.n = dim(G)[1]
	snp.n = dim(G)[2] 
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj(x1, x2, e1)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj(x1, x2, e1)
		G.weight = MatMult(G,weight)
	}

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
	snp.n = dim(G)[2]  
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)
	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj_Gfix(x1, x2, e1, G)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj_Gfix(x1, x2, e1, G)
		G.weight = MatMult(G,weight)
	}

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

	p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
	return(p.continuous.fix)
}




MAGE_RAN.B = function(y, x1, x2, e1, G, n.cov)
{
 	n.fix = n.cov+2 
	sub.n = dim(G)[1]
	snp.n = dim(G)[2]  
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj(x1, x2, e1)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj(x1, x2, e1)
		G.weight = MatMult(G,weight)
	}

	none_matrix = cbind(x1,x2,e1)

	y = trans.ytype(y, none_matrix)

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	G.star = MatMult(M,G.weight)
	S.star = MatMult(M,S)

	p.continuous.ran = MAGE_RAN.body(y.star, G.star, S.star, M, n.fix)
	return(p.continuous.ran)
}



MAGE_FIX.B = function(y, x1, x2, e1, G)
{
	sub.n = dim(G)[1]
	snp.n = dim(G)[2]  
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj_Gfix(x1, x2, e1, G)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj_Gfix(x1, x2, e1, G)
		G.weight = MatMult(G,weight)
	}

	none_matrix = cbind(x1,x2,e1)

	y = trans.ytype(y, none_matrix)

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

	p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
	return(p.continuous.fix)
}



MAGE_RAN_REAL.C = function(y, X, e1, G, n.cov)
{
	n.fix = n.cov+2
	sub.n = dim(G)[1]
	snp.n = dim(G)[2] 
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProjxx(X, e1)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProjxx(X, e1)
		G.weight = MatMult(G,weight)
	}

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	G.star = MatMult(M,G.weight)
	S.star = MatMult(M,S)

	p.continuous.ran = MAGE_RAN.body(y.star, G.star, S.star, M, n.fix)
	return(p.continuous.ran)
}



MAGE_FIX_REAL.C = function(y, X, e1, G)
{
	sub.n = dim(G)[1]
	snp.n = dim(G)[2]  
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj_Gfixx(X, e1, G)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj_Gfixx(X, e1, G)
		G.weight = MatMult(G,weight)
	}
	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

	p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
	return(p.continuous.fix)
}



MAGE_FIX_REAL.B = function(y, X, e1, G)
{
	sub.n = dim(G)[1]
	snp.n = dim(G)[2]  
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProj_Gfixx(X, e1, G)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProj_Gfixx(X, e1, G)
		G.weight = MatMult(G,weight)
	}

	none_matrix = cbind(X,e1)

	y = trans.ytype(y, none_matrix)

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	S.star = MatMult(M,S)

	p.continuous.fix = MAGE_FIX.body(y.star, S.star, M)
	return(p.continuous.fix)
}




MAGE_RAN_REAL.B = function(y, X, e1, G, n.cov)
{
	n.fix = n.cov+2
	sub.n = dim(G)[1]
	snp.n = dim(G)[2] 
	maf = colSums(G)/(sub.n*2)

	weight.num = rep(NA, snp.n)
	index.common = maf>=0.05
	index.rare = maf<0.05
	scaler = dbeta(0.05,1,25)/dbeta(0.05,0.5,0.5)
	weight.num[index.common] = scaler*dbeta(maf[index.common],0.5,0.5)
	weight.num[index.rare] = dbeta(maf[index.rare],1,25)

	if (length(weight.num)==1){
		weight = weight.num
		M = ComputeProjxx(X, e1)
		G.weight = G*weight
	}else{
		weight = diag(weight.num)
		M = ComputeProjxx(X, e1)
		G.weight = MatMult(G,weight)
	}

	none_matrix = cbind(X,e1)

	y = trans.ytype(y, none_matrix)

	S = MatMult(diag(e1),G.weight)
	y.star = M%*%y
	G.star = MatMult(M,G.weight)
	S.star = MatMult(M,S)

	p.continuous.ran = MAGE_RAN.body(y.star, G.star, S.star, M, n.fix)
	return(p.continuous.ran)
}
























