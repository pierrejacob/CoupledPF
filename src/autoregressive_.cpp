#include <RcppEigen.h>
using namespace Rcpp;
using namespace std;
using namespace Eigen;


// [[Rcpp::export]]
NumericVector ar_generate_randomness_cpp(int nparticles, int datalength, int dimension){
  RNGScope scope;
  NumericVector normal_draws = rnorm((1 + datalength) * nparticles * dimension, 0, 1);
  return normal_draws;
}

// [[Rcpp::export]]
NumericVector ar_perturb_randomness_cpp(const NumericVector & randomness, double rho, int dimension){
  RNGScope scope;
  int l = randomness.size();
  NumericVector newrand(l);
  double v = sqrt(1.0 - rho*rho);
  NumericVector normal_draws = rnorm(l, 0, 1);
  for (int i = 0; i < l; i ++){
    newrand(i) =  rho * randomness(i) + v * normal_draws(i);
  }
  return newrand;
}

// [[Rcpp::export]]
NumericMatrix ar_rinit_rcpp(int nparticles, NumericVector theta, NumericVector rand, int dimension){
  NumericMatrix xparticles(dimension, nparticles);
  for (int iparticle = 0; iparticle < nparticles; iparticle ++){
    for (int istate = 0; istate < dimension; istate ++){
      xparticles(istate, iparticle) =  rand((istate * nparticles) + iparticle);
    }
  }
  return xparticles;
}

// [[Rcpp::export]]
NumericMatrix ar_rtransition_rcpp(const NumericMatrix & xparticles, NumericVector theta,
                                     int time, NumericVector rand, int dimension, NumericMatrix & A){
  int nparticles = xparticles.ncol();
  NumericMatrix new_x(dimension, nparticles);
  std::fill(new_x.begin(), new_x.end(), 0.);
  for (int iparticle = 0; iparticle < nparticles; iparticle ++){
    for (int istate = 0; istate < dimension; istate ++){
      for (int j = 0; j < dimension; j ++){
        new_x(istate, iparticle) +=  A(istate, j) * xparticles(j,iparticle);
      }
      new_x(istate, iparticle) += rand(dimension * nparticles * time + (iparticle * dimension) + istate);
    }
  }
  return new_x;
}
