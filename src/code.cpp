#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericVector cortMonteCarlo(const NumericMatrix z,
                                   const NumericMatrix min,
                                   const NumericMatrix max,
                                   const int N) {
  int d = z.nrow();
  int n = z.ncol();
  int D = min.ncol(); //  == 2^d
  double f_l;
  NumericVector lambda_k(D);
  NumericVector lambda_l(D);
  NumericVector f_k(D);
  LogicalMatrix core_checks(D,n);
  NumericMatrix z_boot(d,n);
  NumericVector observed_stat(d);
  NumericMatrix bootstraped_stat(N,d);
  NumericVector p_values(d);

  observed_stat.fill(0.0);
  bootstraped_stat.fill(0.0);
  //p_values.fill(0.0);

  for (int d_rem = 0; d_rem < d; d_rem++){

    f_k.fill(0.0);
    lambda_k.fill(1.0);
    lambda_l.fill(1.0);
    core_checks.fill(true);

    // Compute the statistic :
    for (int n_leave = 0; n_leave < D; n_leave++){
      f_l = 0.0;
      for (int dim = 0; dim < d; dim++){
        lambda_l(n_leave) *= max(dim,n_leave) - min(dim,n_leave);
      }
      lambda_k(n_leave) = lambda_l(n_leave) / (max(d_rem,n_leave) - min(d_rem,n_leave));

      for (int n_obs = 0; n_obs < n; n_obs++){
        for(int dim = 0; dim < d; dim++){
          if(d_rem != dim){
            core_checks(n_leave,n_obs) = core_checks(n_leave,n_obs) && (min(dim,n_leave) <= z(dim,n_obs)) && (z(dim,n_obs) < max(dim,n_leave));
          }
        }
        if(core_checks(n_leave,n_obs)){
          f_k(n_leave) += 1.0/n;
        }
        if(core_checks(n_leave,n_obs) & (min(d_rem,n_leave) <= z(d_rem,n_obs)) && (z(d_rem,n_obs) < max(d_rem,n_leave))){
          f_l += 1.0/n;
        }
      }
      observed_stat(d_rem) += (f_l*f_l)/lambda_l(n_leave) - 2 * (f_l*f_k(n_leave))/lambda_k(n_leave);
    }

    // Now bootstrap the same thing :
    z_boot = Rcpp::clone(z);
    for (int n_boot = 0; n_boot < N; n_boot++){
      z_boot(d_rem,_) = runif(n); // The bootstrap is here.
      for (int n_leave = 0; n_leave < D; n_leave++){
        f_l = 0.0;
        for (int n_obs = 0; n_obs < n; n_obs++){
          if(core_checks(n_leave,n_obs) && (min(d_rem,n_leave) <= z_boot(d_rem,n_obs)) && (z_boot(d_rem,n_obs) < max(d_rem,n_leave))){
            f_l += 1.0/n;
          }
        }
        bootstraped_stat(n_boot,d_rem) += (f_l*f_l)/lambda_l(n_leave) - 2 * (f_l*f_k(n_leave))/lambda_k(n_leave);
      }
    }
  }
  for(int d_rem =0; d_rem < d; d_rem++){
    p_values(d_rem) = mean(observed_stat(d_rem) <= bootstraped_stat(_,d_rem));
  }
  return(p_values);
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

// [[Rcpp::export]]
double lossFunc(const NumericVector bp,
                const NumericMatrix bin_repr,
                const NumericMatrix z){

  int D = bin_repr.ncol();
  int d = bin_repr.nrow();
  int n = z.ncol();
  double n_pts, vol, loss=0;
  Rcpp::NumericVector min(d);
  Rcpp::NumericVector max(d);

  for (int n_box=0; n_box < D; n_box++){
    vol = 1.0;
    n_pts= 0.0;
    // Compute sides and volumes :
    for (int dim=0; dim < d; dim++){
      min(dim) = bp(dim)*bin_repr(dim,n_box);
      max(dim) = std::pow(bp(dim),1.0-bin_repr(dim,n_box));
      vol = vol*(max(dim) - min(dim));
    }
    if(vol > 0){
      // check if observations are inside the box :
      for(int n_obs=0; n_obs<n; n_obs++){
        if(is_true(all(min<=z(_,n_obs))) && is_true(all(z(_,n_obs)<max))){
          n_pts += 1.0;
        }
      }
      loss -=std::pow(n_pts/n,2)/vol;
    }
  }
  return(loss);
}






