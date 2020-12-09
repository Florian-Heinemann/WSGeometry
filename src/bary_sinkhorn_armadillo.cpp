#include <RcppArmadillo.h>
#include <ctime>
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
//[[Rcpp::plugins(openmp)]]
// [[Rcpp::depends("RcppArmadillo")]]


//[[Rcpp::export]]
Rcpp::List bary_sinkhorn_arma(arma::mat weights,arma::mat frechet,int maxIter, double lambda,arma::mat C,double thresh, int threads) {
  int N = weights.n_cols;
  int M = weights.n_rows;
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  
  arma::mat K= exp((-1*C)/lambda);
  arma::mat KT=trans(K);
  arma::mat b=arma::ones(M,N);
  arma::mat a=arma::ones(M,N);
  arma::mat p=arma::ones(M,N);
  arma::mat pvec=arma::ones(M,1);
  arma::mat pOld=arma::ones(M,1)*1000;
  int iter=0;
  double p_change=1000;
  while ((iter<maxIter) && (p_change>thresh)){
    //Rcout << iter+1<< ",";
    //a update loop
    a=weights/(K*b);
    a=KT*a;
    //p update loop
    p=a;
    pvec=pow(p.col(0),frechet(0));
    for (int i=1;i<N;i++){
      pvec=pvec%pow(p.col(i),frechet(i));
    }
    //update b loop
    #pragma omp parallel for 
    for (int i=0;i<N;i++){
      b.col(i)=pvec/a.col(i);
    }
    p_change=arma::accu(abs((pvec-pOld)));
    pOld=pvec;
    iter+=1;
  }
  //Rcout << "\n";
  return  Rcpp::List::create(Rcpp::Named("bary") =pvec,Rcpp::Named("iterations")=iter);
}

