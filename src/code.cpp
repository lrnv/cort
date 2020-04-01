#include <Rcpp.h>
using namespace Rcpp;

// This a simple rewritting of my bootstrap function
double compute_stat(const NumericMatrix z,
                    const NumericMatrix min,
                    const NumericMatrix max,
                    const int d_rem) {

    int d = z.nrow();
    int n = z.ncol();
    int D = min.ncol(); //  == 2^d
    bool temp_l, temp_k;
    double check, statistic, lambda_l, lambda_k, f_l, f_k;

    statistic = 0.0;
    for (int n_leave = 0; n_leave < D; n_leave++){
      f_l = 0;
      f_k = 0;
      lambda_l = 1;
      lambda_k = 1;
      for (int dim = 0; dim < d; dim++){
        lambda_l *= max(dim,n_leave) - min(dim,n_leave);
      }
      lambda_k = lambda_l / (max(d_rem,n_leave) - min(d_rem,n_leave));
      for (int n_obs = 0; n_obs < n; n_obs++){
        temp_l = true;
        temp_k = true;
        for(int dim = 0; dim < d; dim++){
          check = (min(dim,n_leave) <= z(dim,n_obs)) && (z(dim,n_obs) < max(dim,n_leave));
          temp_l = temp_l && check;
          if(d_rem != dim){
            temp_k = temp_k && check;
          }
        }
        if(temp_l){
          f_l += 1.0/n;
        }
        if(temp_k){
          f_k += 1.0/n;
        }
      }
      statistic += (f_l*f_l)/lambda_l - 2 * (f_l*f_k)/lambda_k;
    }
    return(statistic);
  }


// Now that we have a way of computing the statistic, we need to bootstrap the thing.

// [[Rcpp::export]]
Rcpp::NumericMatrix cortMonteCarlo(const NumericMatrix z,
                                   const NumericMatrix min,
                                   const NumericMatrix max,
                                   const int N) {
  // For each d_rem in 1:d, compute the statistic:
  int d = z.nrow();
  int n = z.ncol();
  NumericMatrix z_boot(d,n);
  NumericMatrix result(N+1,d);
  result.fill(0.0);

  for (int d_rem = 0; d_rem < d; d_rem++){
    z_boot = Rcpp::clone(z);
    result(0,d_rem) += compute_stat(z,min,max,d_rem);

    for (int n_boot = 1; n_boot <= N; n_boot++){
      z_boot(d_rem,_) = runif(n);
      result(n_boot,d_rem) = compute_stat(z_boot,min,max,d_rem);
    }
  }
  return(result);
}


// [[Rcpp::export]]
double quadProd(const NumericMatrix a,
                const NumericMatrix b,
                const NumericVector kern,
                const NumericMatrix other_a,
                const NumericMatrix other_b,
                const NumericVector other_kern) {

  int d = a.ncol();
  int n = a.nrow();
  int other_n = other_a.nrow();
  double temp;
  double rez = 0.0;

  for (int i=0; i<n; i++){
    for (int j=0; j<other_n; j++){
      temp = kern(i)*other_kern(j);
      for (int dim=0; dim<d; dim++){
        temp *= std::max(std::min(b(i,dim),other_b(j,dim)) - std::max(a(i,dim),other_a(j,dim)),0.0);
      }
      rez += temp;
    }
  }
  return(rez);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix normMatrix(const List as,
                               const List bs,
                               const List kernels){

  int n = as.length();
  Rcpp::NumericMatrix result(n,n);
  for(int i=0; i< n; i++){
    for(int j=0; j< n; j++){
      if(i <= j){
        result(i,j) = quadProd(as[i], bs[i], kernels[i], as[j], bs[j], kernels[j]);
        result(j,i) = result(i,j);
      }

    }
  }
  return(result);
}


































