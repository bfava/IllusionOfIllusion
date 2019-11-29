
#include <RcppArmadilloExtensions/sample.h>
#include <iostream>
#include <algorithm>
// [[Rcpp::depends(RcppArmadillo)]]



using namespace Rcpp;
//using namespace arma;

// [[Rcpp::export]]
NumericVector pqr(float sigma2, int k, NumericVector lq, NumericVector l1q, NumericVector lr2, NumericVector l1r2, NumericVector c1, NumericVector area, float cterm, NumericVector z, int a = 1, int b = 1, int A = 1, int B = 1) {
  
  int sumz = sum(z);
  
  NumericVector lres = (- (cterm * c1) / sigma2) +
                       ((3 * (sumz / 2) + a - 1) * lq) +
                       ((k - sumz + b - 1) * l1q) +
                       ((A - 1 - sumz / 2) * lr2) +
                       ((sumz / 2 + B - 1) * l1r2);
                        
  double max_lres = max(lres);
  NumericVector res = exp(lres - max_lres);
  res = res * area;
  
  return res;
  
}


// [[Rcpp::export]]
arma::vec pz_t(arma::vec z, arma::vec zp, float r2, float q, float lambda2_bar, arma::vec y, arma::mat X, arma::mat U, arma::vec phi, int Tn, int k){
  
  float lq = log(q);
  float l1q = log(1 - q);
  float l2 = log(2);
  float lg = log(r2) - (log(k) + lq + log(lambda2_bar) + log(1 - r2));
  float g2 = exp(lg);
  
  arma::vec ytil = y - U * phi;
  
  int tauz = sum(z);
  arma::mat Xtil = X.cols(find(z == 1));
  arma::vec ones(tauz);
  ones = ones.ones();
  ones = ones / g2;
  arma::mat Wtil = Xtil.t() * Xtil + diagmat(ones);
  arma::vec betatilhat = inv(Wtil) * Xtil.t() * ytil;
  
  int tauzp = sum(zp);
  arma::mat Xtilp = X.cols(find(zp == 1));
  arma::vec onesp(tauzp);
  onesp = onesp.ones();
  onesp = onesp / g2;
  arma::mat Wtilp = Xtilp.t() * Xtilp + diagmat(onesp);
  arma::vec betatilhatp = inv(Wtilp) * Xtilp.t() * ytil;
  
  //double ldet = real(log_det(Wtilp));
  //double ldet = 2 * sum(log(diagmat(chol(Wtilp))));
  
  arma::vec ans = exp(
    (tauzp * lq) +
      ((k - tauzp) * l1q) -
      ((tauzp / 2) * lg) -
      (real(log_det(Wtilp)) / 2) -
      ((Tn / 2) * (log(ytil.t() * ytil - betatilhatp.t() * Wtilp * betatilhatp) - l2)) -
      
      ((tauz * lq) +
      ((k - tauz) * l1q) -
      ((tauz / 2) * lg) -
      (real(log_det(Wtil)) / 2) -
      ((Tn / 2) * (log(ytil.t() * ytil - betatilhat.t() * Wtil * betatilhat) - l2)))
  );
  
  return ans;
  
}


// [[Rcpp::export]]
arma::vec pz_t_nphi(arma::vec z, arma::vec zp, float r2, float q, float lambda2_bar, arma::vec y, arma::mat X, int Tn, int k){
  
  float lq = log(q);
  float l1q = log(1 - q);
  float l2 = log(2);
  float lg = log(r2) - (log(k) + lq + log(lambda2_bar) + log(1 - r2));
  float g2 = exp(lg);
  
 //float g2 = (1 / (k * q * lambda2_bar)) * (r2 / (1 - r2));
  arma::vec ytil = y;
  
  int tauz = sum(z);
  arma::mat Xtil = X.cols(find(z == 1));
  arma::vec ones(tauz);
  ones = ones.ones();
  ones = ones / g2;
  arma::mat Wtil = Xtil.t() * Xtil + diagmat(ones);
  arma::vec betatilhat = inv(Wtil) * Xtil.t() * ytil;
  
  int tauzp = sum(zp);
  arma::mat Xtilp = X.cols(find(zp == 1));
  arma::vec onesp(tauzp);
  onesp = onesp.ones();
  onesp = onesp / g2;
  arma::mat Wtilp = Xtilp.t() * Xtilp + diagmat(onesp);
  arma::vec betatilhatp = inv(Wtilp) * Xtilp.t() * ytil;
  

  
  //double ldet = real(log_det(Wtilp));
  //double ldet = 2 * sum(log(diagmat(chol(Wtilp))));
  
  arma::vec ans = exp(
    (tauzp * lq) +
      ((k - tauzp) * l1q) -
      ((tauzp / 2) * lg) -
      (real(log_det(Wtilp)) / 2) -
      ((Tn / 2) * (log(ytil.t() * ytil - betatilhatp.t() * Wtilp * betatilhatp) - l2)) -
      
      ((tauz * lq) +
      ((k - tauz) * l1q) -
      ((tauz / 2) * lg) -
      (real(log_det(Wtil)) / 2) -
      ((Tn / 2) * (log(ytil.t() * ytil - betatilhat.t() * Wtil * betatilhat) - l2)))
  );
  
  return ans;
  
}


// [[Rcpp::export]]
arma::vec refz(arma::vec z, float r2, float q, arma::vec y, arma::mat X, arma::mat U, arma::vec phi, int Tn, int k){
  
  arma::vec zp(z.size());
  arma::vec probzp(1);
  IntegerVector samp(2);
  NumericVector probs(2);
  probs(0) = 1;
  IntegerVector chg(1);
  int n = z.size();
  IntegerVector posv(n);
  std::iota(posv.begin(), posv.end(), 0);
  std::random_shuffle(posv.begin(), posv.end());
  int pos;
  
  for (int j = 0; j < n; j++){
    
    pos = posv(j);
    if(!((sum(z) == 1) & (z(pos) == 1))){
      
      zp = z;
      zp(pos) = 1 - zp(pos);
      probzp = pz(z, zp, r2, q, y, X, U, phi, Tn, k);
      samp(0) = z(pos);
      samp(1) = zp(pos);
      probs(1) = probzp(0);
      
      chg = RcppArmadillo::sample(samp, 1, false, probs);
      z(pos) = chg(0);
    
    }
    
  }
  
  return z;

}
