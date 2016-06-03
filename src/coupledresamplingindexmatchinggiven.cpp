#include <RcppEigen.h>
#include "systematic.h"

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerVector indexmatching_given_cpp(int ndraws, NumericVector w1, NumericVector w2, NumericVector uniforms,
                                      IntegerVector ancestors_ref){
  int nparticles = w1.size();
  double unif_systematic1 = uniforms(ndraws);
  
  NumericVector nu(nparticles);
  NumericVector couplingprobabilities(nparticles);
  
  double alpha = 0;
  for (int i = 0; i < nparticles; i ++){
    nu(i) = std::min(w1(i), w2(i));
    alpha += nu(i);
  }
  for (int i = 0; i < nparticles; i ++){
    couplingprobabilities(i) = nu(ancestors_ref(i)) / w1(ancestors_ref(i));
  }
  // minorization proba
  NumericVector mu = nu / alpha;
  // residuals
  NumericVector R2 = (w2 - nu) / (1 - alpha);
  IntegerVector ancestors2(ndraws);
  //
  LogicalVector coupled(ndraws);
  int ncoupled = 0;
  for (int i = 0; i < ndraws; i ++){
    if (uniforms(i) < couplingprobabilities(i)){
      coupled(i) = true;
      ncoupled ++;
    } else {
      coupled(i) = false;
    }
  }
  IntegerVector a_R2;
  if (ncoupled < ndraws){
    a_R2 = systematic_resampling_n_(R2, ndraws - ncoupled, unif_systematic1);
  }
  int coupledcounter = 0;
  int uncoupledcounter = 0;
  for (int i = 0; i < ndraws; i ++){
    if (coupled(i)){
      ancestors2(i) = ancestors_ref(i);
      coupledcounter ++;
    } else {
      ancestors2(i) = a_R2(uncoupledcounter);
      uncoupledcounter ++;
    }
  }
  return ancestors2;
}

