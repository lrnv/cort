#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector dcbCopula(const NumericMatrix u,
                              const NumericMatrix x,
                              const NumericVector m){

  int n_obs = u.nrow();
  int dim = x.ncol();
  int n_pts = x.nrow();
  Rcpp::NumericMatrix floor_x(x);
  Rcpp::NumericMatrix floor_u(u);
  Rcpp::NumericVector rez(n_obs);

  // coerce grid :
  for(int i = 0; i < n_pts; i++){
    floor_x(i,_) = floor(x(i,_)*m);
  }

  for (int j = 0; j < n_obs; j++){
    // coerce points to inf_points :
    floor_u(j,_) = floor(u(j,_)*m);
    for(int i = 0; i < n_pts; i++){
      if(is_true(all(floor_u(j,_) == floor_x(i,_)))){
        rez(j) += 1.0/n_pts;
      }
    }
  }
  for(int d =0; d < dim; d++){
    rez = rez * m(d);
  }
  return(rez);
}

//
// x = as.matrix(copula@data)
//   seuil_inf = boxes_from_points(x,copula@m)
//   seuil_inf_u = boxes_from_points(u,copula@m)
//   return(sapply(1:nrow(u), function(i){mean(apply(t(seuil_inf) == seuil_inf[i,],2,prod))}))
//
