

// load Rcpp
#include <RcppArmadillo.h>
#include <omp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]



// [[Rcpp::export]] 
arma::mat MatMult(arma::mat A, arma::mat B){
	mat C;

	C = A*B;

	return C;
}


// [[Rcpp::export]]
arma::mat Compute5PCs(arma::mat X,int top = 5){
    mat U;
    vec s;
    mat V;
    svd(U,s,V,X);
    
       int n = U.n_cols;
         int m = s.n_elem;
         vec z(n-m, fill::zeros);
         vec c;
         c = join_cols(s,z);

    mat PCs = U*diagmat(c);
    return PCs.cols(0,top-1);
}



// [[Rcpp::export]]
arma::mat ComputePCs(arma::mat X,int top = 10){
    mat U;
    vec s;
    mat V;
    svd(U,s,V,X);
    
       int n = U.n_cols;
         int m = s.n_elem;
         vec z(n-m, fill::zeros);
         vec c;
         c = join_cols(s,z);

    mat PCs = U*diagmat(c);
    return PCs.cols(0,top-1);
}


// [[Rcpp::export]]
arma:: vec comeigen_fixed_Gfix(arma:: mat& S,arma::vec& q,arma::mat& Sc,arma::mat& M)
{
  arma::mat Sinv = arma::zeros(2,2);
  Sinv = arma::inv(S);
  arma::mat Vsigma = Sinv(0,0)*Sc+Sinv(0,1)*M;
  arma::vec eigval1;
  arma::mat eigvec1;
  arma::eig_sym(eigval1,eigvec1,Vsigma);
  return eigval1;
}









// [[Rcpp::export]]
arma:: vec comeigen_fixed(arma:: mat& S,arma::vec& q,arma::mat& Gc,arma::mat& Sc,arma::mat& M)
{
  arma::mat Sinv = arma::zeros(3,3);
  Sinv = arma::inv(S);
  arma::mat Vsigma = Sinv(1,0)*Gc+Sinv(1,1)*Sc+Sinv(1,2)*M;
  arma::vec eigval1;
  arma::mat eigvec1;
  arma::eig_sym(eigval1,eigvec1,Vsigma);
  return eigval1;
}







// [[Rcpp::export]]
arma::mat GetLinearKernel(arma::mat G){ //The data matrix is n*p
  double p = G.n_cols;
  return G*G.t()/p; //Notice G is n*p, then G.t() is p*n; Thus, this return a n*n kernel.
  //In our first case, n =1000, p =2
}

// [[Rcpp::export]] 
arma:: mat CovComb(arma:: vec X1, vec X2){//Assume both n*n
  int n = X1.n_elem;
  mat cov = zeros(n,2);
  cov.col(0) = X1;
  cov.col(1) = X2;
  return cov;
  
}
// [[Rcpp::export]] 
arma:: mat mat_inv_3(arma::mat S){
  mat Sinv = zeros(3,3);
  double det = S(0,0)*(S(1,1)*S(2,2) - S(1,2)*S(2,1)) - S(0,1)*(S(1,0)*S(2,2)-S(1,2)*S(2,0))+S(0,2)*(S(1,0)*S(2,1)-S(1,1)*S(2,0));
  Sinv(0,0) = S(1,1)*S(2,2) - S(1,2)*S(2,1);
  Sinv(0,1) = -(S(1,0)*S(2,2) - S(1,2)*S(2,0));
  Sinv(0,2) = S(1,0)*S(2,1) - S(1,1)*S(2,0);
  Sinv(1,0) = -(S(0,1)*S(2,2) - S(0,2)*S(2,1));
  Sinv(1,1) = S(0,0)*S(2,2) - S(0,2)*S(2,0);
  Sinv(1,2) = -(S(0,0)*S(2,1) - S(2,0)*S(0,1));
  Sinv(2,0) = S(0,1)*S(1,2) - S(0,2)*S(1,1);
  Sinv(2,1) = -(S(0,0)*S(1,2) - S(0,2)*S(1,0));
  Sinv(2,2) = S(0,0)*S(1,1) - S(1,0)*S(0,1);
  Sinv = Sinv.t();
  Sinv = Sinv.t()/det;
  return Sinv;
  
}


// [[Rcpp::export]]
arma::mat GetS_est_Gfix(arma::mat Sc, mat M){
  mat S = zeros(2,2);//The data matrix is n*p
  S(0,0) = trace(Sc*Sc);
  S(0,1) = trace(Sc*M);
  S(1,0) = S(0,1);
  S(1,1) = trace(M*M);
  return S;
}


// [[Rcpp::export]]
arma::mat GetS_est(arma::mat Gc,mat Sc, mat M){
  mat S = zeros(3,3);//The data matrix is n*p
  S(0,0) = trace(Gc*Gc);
  S(0,1) = trace(Gc*Sc);
  S(0,2) = trace(Gc*M);
  S(1,0) = S(0,1);
  S(1,1) = trace(Sc*Sc);
  S(1,2) = trace(Sc*M);
  S(2,0) = S(0,2);
  S(2,1) = S(1,2);
  S(2,2) = trace(M*M);
  return S;
}


// [[Rcpp::export]]
arma::mat GetS_est4(arma::mat Gc,mat Sc, mat M, mat Uc){
  mat S = zeros(4,4);//The data matrix is n*p
  S(0,0) = trace(Gc*Gc);
  S(0,1) = trace(Gc*Sc);
  S(0,2) = trace(Gc*Uc);
  S(0,3) = trace(Gc*M);
  S(1,0) = S(0,1);
  S(1,1) = trace(Sc*Sc);
  S(1,2) = trace(Sc*Uc);
  S(1,3) = trace(Sc*M);
  S(2,0) = S(0,2);
  S(2,1) = S(1,2);
  S(2,2) = trace(Uc*Uc);
  S(2,3) = trace(Uc*M);
  S(3,0) = S(0,3);
  S(3,1) = S(1,3);
  S(3,2) = S(2,3);
  S(3,3) = trace(M*M); 
  return S;
}

// [[Rcpp::export]]
arma::mat GetQ_est_Gfix(arma::mat Sc, mat M,vec yc){
  vec q = zeros(2);//The data matrix is n*p
  q(0) = as_scalar(yc.t()*Sc*yc);//Represents gene main effect
  q(1) = as_scalar(yc.t()*M*yc);//Represents error term
  return q;
}


// [[Rcpp::export]]
arma::mat GetQ_est(arma::mat Gc,mat Sc, mat M,vec yc){
  vec q = zeros(3);//The data matrix is n*p
  q(0) = as_scalar(yc.t()*Gc*yc);//Represents ge interaction component
  q(1) = as_scalar(yc.t()*Sc*yc);//Represents gene main effect
  q(2) = as_scalar(yc.t()*M*yc);//Represents error term
  return q;
}


// [[Rcpp::export]]
arma::mat GetQ_est4(arma::mat Gc,mat Sc, mat Uc, mat M,vec yc){
  vec q = zeros(4);//The data matrix is n*p
  q(0) = as_scalar(yc.t()*Gc*yc);//Represents gene main effect
  q(1) = as_scalar(yc.t()*Sc*yc);//Represents interaction effect
  q(2) = as_scalar(yc.t()*Uc*yc);// Represents covariate effect
  q(3) = as_scalar(yc.t()*M*yc);//Represents error term
  return q;
}


// [[Rcpp::export]] 
arma::mat ComputeProj_Gfix(arma::vec X1, vec X2, vec E, mat Gc){
  int n = E.n_elem;
  int m = Gc.n_cols;
  mat b = zeros(n,m+4);
  b.col(0) = ones<vec>(n); 
b.col(1) = X1;
b.col(2) = X2;
b.col(3) = E;
	for(int i = 4; i<m+4; ++i){
	b.col(i) = Gc.col(i-4);
	}
  mat M = eye<mat>(n,n)-b*(pinv(b.t()*b))*b.t(); 
  return M;
}


// [[Rcpp::export]] 
arma::mat ComputeProj_Gfixx(arma::mat X, vec E, mat Gc){
  int n = E.n_elem;
  int m = Gc.n_cols;
  int k = X.n_cols;
  mat b = zeros(n,m+k+2);
  b.col(0) = ones<vec>(n); 

	for(int num = 1; num<k+1; ++num){
	b.col(num) = X.col(num-2);
	}

   b.col(k+1) = E;

	for(int i = k+2; i<m+k+2; ++i){
	b.col(i) = Gc.col(i-4);
	}
  mat M = eye<mat>(n,n)-b*(pinv(b.t()*b))*b.t(); 
  return M;
}









// [[Rcpp::export]] 
arma::mat ComputeProj(arma::vec X1, vec X2, vec E){
  int n = E.n_elem;
  mat b = zeros(n,4);
  b.col(0) = ones<vec>(n); b.col(1) = X1; b.col(2) = X2; b.col(3) = E;
  mat M = eye<mat>(n,n)-b*(inv(b.t()*b))*b.t(); 
  return M;
}

// [[Rcpp::export]] 
arma::mat ComputeProj_randomX(arma:: vec E){
  int n = E.n_elem;
  mat b = zeros(n,2);
  b.col(0) = ones<vec>(n); b.col(1) = E;
  mat M = eye<mat>(n,n)-b*(inv(b.t()*b))*b.t(); 
  return M;
}
// [[Rcpp::export]] 
arma::vec ComputeDelta(arma::mat S, vec q){
  mat Sinv = inv(S);
  vec delta = Sinv*q;
  return delta;
}
// [[Rcpp::export]] 

arma::mat Standard(arma::mat S){
  int p = S.n_cols;
  int i ;
  vec x; 
  for (i=0; i<p; i++){
    x = S.col(i);
    S.col(i) = (x - mean(x))/stddev(x);
  }
  return S;
}
// [[Rcpp::export]] 
arma::mat Standard_vec(arma::vec S){
  vec x = (S - mean(S))/stddev(S);
  return S;
}
// [[Rcpp::export]] 

arma::mat Center(arma::mat S){
  int p = S.n_cols;
  int i ;
  vec x; 
  for (i=0; i<p; i++){
    x = S.col(i);
    S.col(i) = (x - mean(x));
  }
  return S;
}


//Following the instruction by Zhou's paper, the SNPs should be standardized.
// [[Rcpp::export]]  
arma:: mat Sigma_mapit(arma::mat S,vec yc, vec delta, mat Gc, mat Sc, mat M){
  mat Vy = delta(0)*Gc+delta(1)*Sc+delta(2)*M;
  mat q = zeros(3,3); 
  q(0,0) = as_scalar(2*yc.t() * Gc* Vy* Gc*yc);
  q(0,1) = as_scalar(2*yc.t() * Gc* Vy* Sc*yc);
  q(0,2) = as_scalar(2*yc.t() * Gc* Vy* M*yc);
  q(1,0) =q(0,1);
  q(1,1) = as_scalar(2*yc.t() * Sc* Vy* Sc*yc);q(1,2)=as_scalar(2*yc.t() * Sc* Vy* M *yc);
  q(2,0) = q(0,2);
  q(2,1) = q(1,2);
  q(2,2) = as_scalar(2*yc.t() * M* Vy* M*yc);
  mat Sinv = zeros(3,3);
  Sinv = inv(S);
  mat sigma = Sinv*q*Sinv.t();
  return sigma;
}
// [[Rcpp::export]]
arma:: mat Sigma_mapit_rancov(arma::mat S,vec yc, vec delta, mat Gc, mat Sc, mat Uc, mat M){
  mat Vy = delta(0)*Gc+delta(1)*Sc+delta(2)*Uc+delta(3)*M;
  mat q = zeros(4,4); 
  q(0,0) = as_scalar(2*yc.t() * Gc* Vy* Gc*yc);
  q(0,1) = as_scalar(2*yc.t() * Gc* Vy* Sc*yc);
  q(0,2) = as_scalar(2*yc.t() * Gc* Vy* Uc*yc);
  q(0,3) = as_scalar(2*yc.t() * Gc* Vy* M*yc);
  q(1,0) =q(0,1);
  q(1,1) = as_scalar(2*yc.t() * Sc* Vy* Sc*yc);
  q(1,2)=as_scalar(2*yc.t() * Sc* Vy* Uc *yc);
  q(1,3)=as_scalar(2*yc.t() * Sc* Vy* M *yc);
  q(2,0) = q(0,2);
  q(2,1) = q(1,2);
  q(2,2) = as_scalar(2*yc.t() * Uc* Vy* Uc*yc);
  q(2,3) = as_scalar(2*yc.t() * Uc* Vy* M*yc);
  q(3,0) = q(0,3);
  q(3,1) = q(1,3);
  q(3,2) = q(2,3);
  q(3,3) = as_scalar(2*yc.t() * M* Vy* M*yc);
  mat Sinv = zeros(4,4);
  Sinv = inv(S);
  mat sigma = Sinv*q*Sinv.t();
  return sigma;
}

// [[Rcpp::export]]  
arma:: vec eig_null(arma:: mat S, vec q,mat Gc, mat Sc, mat M)
{
  vec q_null = zeros(2);
  q_null(0) = q(0);
  q_null(1) = q(2); 
  mat Sinv = inv(S);
  mat S_null = zeros(2,2);
  S_null(0,0) = S(0,0);
  S_null(0,1) = S(0,2);
  S_null(1,0) = S(2,0);
  S_null(1,1) = S(2,2);
  vec delta_null = inv(S_null)*q_null;// 2 by 1 vector
  mat Vy_null = delta_null(0)*Gc+delta_null(1)*M; //Submatrix can be used drectly. 
  vec eigval1;
  mat eigvec1; 
  eig_sym(eigval1,eigvec1,Vy_null);
  mat poseig = eigvec1.cols(find(eigval1>0));
  mat v_half = (poseig)*diagmat(sqrt(eigval1(find(eigval1>0))))*trans(poseig); //has to be positive definite
  mat Vsigma = Sinv(1,0)*Gc+Sinv(1,1)*Sc+Sinv(1,2)*M;
  vec eigval2;
  eig_sym(eigval2,v_half*Vsigma*v_half);
  return eigval2;
}


