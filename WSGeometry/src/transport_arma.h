#include <RcppArmadillo.h>
#ifdef _OPENMP
# include <omp.h>
#endif
#include <iostream>
#include <vector>
#include "network_simplex_simple.h"
#include <stdio.h>
using namespace Rcpp;
using namespace std;
using namespace lemon;


//[[Rcpp::export]]
arma::mat transport_network_dual_arma(arma::mat a, arma::mat b, arma::mat C){
  struct TsFlow {
    int from, to;
    double amount;
  };
  typedef FullBipartiteDigraph Digraph;
  //DIGRAPH_TYPEDEFS(FullBipartiteDigraph);
  
  
  int64_t n1 = a.n_rows;
  int64_t n2 = b.n_rows; 
  std::vector<double> weights1(n1), weights2(n2);
  
  Digraph di(n1, n2);
  NetworkSimplexSimple<Digraph, double, double, int64_t> net(di, true, n1 + n2, n1*n2);
  
  int64_t idarc = 0;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      double d =C(i,j);
      net.setCost(di.arcFromId(idarc), d);
      idarc++;
    }
  }
  
  for (int i = 0; i < n1; i++) {
    weights1[di.nodeFromId(i)] = a(i,0);
  }
  for (int i = 0; i < n2; i++) {
    weights2[di.nodeFromId(i)] = (-1)*b(i,0);
  }
  net.supplyMap(&weights1[0], n1, &weights2[0], n2);
  net.run();
  arma::mat Tpot=arma::zeros(n1,1);
  for (int64_t i = 0; i < (n1); i++) {
    Tpot(i,0)=(-1)*net.potential(i);
  }
  return(Tpot);
}

//[[Rcpp::export]]

arma::mat transport_network_primal_arma(arma::mat a, arma::mat b, arma::mat C){
  struct TsFlow {
    int from, to;
    double amount;
  };
  typedef FullBipartiteDigraph Digraph;
  //DIGRAPH_TYPEDEFS(FullBipartiteDigraph);
  
  
  int64_t n1 = a.n_rows;
  int64_t n2 = b.n_rows; 
  std::vector<double> weights1(n1), weights2(n2);
  
  Digraph di(n1, n2);
  NetworkSimplexSimple<Digraph, double, double, int64_t> net(di, true, n1 + n2, n1*n2);
  
  int64_t idarc = 0;
  for (int i = 0; i < n1; i++) {
    for (int j = 0; j < n2; j++) {
      double d =C(i,j);
      net.setCost(di.arcFromId(idarc), d);
      idarc++;
    }
  }
  
  for (int i = 0; i < n1; i++) {
    weights1[di.nodeFromId(i)] = a(i,0);
  }
  for (int i = 0; i < n2; i++) {
    weights2[di.nodeFromId(i)] = (-1)*b(i,0);
  }
  net.supplyMap(&weights1[0], n1, &weights2[0], n2);
  net.run();
  std::vector<TsFlow> flow;
  flow.reserve(n1 + n2 - 1);
  arma::mat Tplan=arma::zeros(n1,n2);
  int count=0;
  for (int64_t i = 0; i < n1; i++) {
    for (int64_t j = 0; j < n2; j++)
    {
      TsFlow f;
      f.amount = net.flow(di.arcFromId(i*n2 + j));
      Tplan(i,j)=f.amount;
      count+=1;
    }
  }
  return(Tplan);
}
