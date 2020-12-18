#include <RcppArmadillo.h>
#include <ctime>
#ifdef _OPENMP
# include <omp.h>
#endif
#include "transport_arma.h"
using namespace Rcpp;
using namespace std;
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(openmp)]]


//Helper tools
//[[Rcpp::export]]
double rand_cxx(){
  double u= unif_rand();
  return(u);
 // random_device rd; 
//  mt19937 mt_rand(rand());
//  return(mt_rand()/static_cast<double>(mt_rand.max()));
}
//[[Rcpp::export]]
int index_cxx(arma::mat x,double y){
  arma::mat sum_vec=arma::cumsum(x)*(1/arma::sum(x))-y;
  const arma::Mat<unsigned int> ind_vec=find(sum_vec >= 0,1,"first");
  //arma::colvec ind  =find(vec > 0,0,"first").;
  //arma::mat ind=find(vec>0,0,"first");
  return(ind_vec(0));
}
//[[Rcpp::export]]
double norm(arma::mat x, arma::mat y){
  arma::mat out=arma::zeros(1,1);
  out=arma::sum(abs((x-y)));
  return(out(0));
}
//[[Rcpp::export]]
arma::mat expo(arma::mat a){
  return(exp(a));
}


//[[Rcpp::export]]
arma::mat gen_cost(const arma::mat A,const arma::mat B){
  arma::mat x=arma::zeros(1,A.n_rows);
  arma::mat z=arma::zeros(1,B.n_rows);
  int nA=A.n_rows;
  int nB=B.n_rows;
  for (int i=0;i<nA;i++){
    x(i)=arma::accu(A.row(i)%A.row(i));
  }
  for (int i=0;i<nB;i++){
    z(i)=arma::accu(B.row(i)%B.row(i));
  }
  arma::mat Cmat=(arma::trans(z)*arma::ones(1,A.n_rows)+arma::ones(B.n_rows,1)*x-(2.0*B*arma::trans(A)));
  return Cmat;
}

//[[Rcpp::export]]
Rcpp::List wsbary_cxx_armaP(const Rcpp::List weightsR,arma::mat positions1,const Rcpp::List positionssetR,const arma::mat frechet_weights,const bool fixed_support,const int maxIter,const int weights_maxIter,const int pos_maxIter ,const double stepsize,const double thresh,const bool headstart, const int headstartlength, int threads) {
  double n=weightsR.size();
  int n2=weightsR.size();
  const double m=positions1.n_rows;
  arma::mat a=arma::ones(m,1)/m;
  arma::mat duals=arma::zeros(m,n2);
  arma::mat s=arma::zeros(n2,1);
  arma::mat ahat=arma::ones(m,1)/m;
  arma::mat atilde=arma::ones(m,1)/m;
  arma::mat alpha=arma::ones(m,1)/m;
  arma::mat aold=arma::ones(m,1)*100.0;
  arma::mat positionsOld=arma::ones(m,positions1.n_cols)*10000.0;
  arma::mat achange=arma::zeros(1,1);
  arma::mat pchange=arma::zeros(1,1);
  arma::mat z;
  arma::mat x;
  arma::mat onesA;
  arma::mat onesB;
  arma::mat Cmat;
  arma::mat positions2;
  arma::mat b;
  arma::mat frechet=arma::diagmat(frechet_weights);
  arma::mat aVeryOld=arma::ones(m,1)*100.0;
  arma::mat positionsVeryOld=arma::ones(m,positions1.n_cols)*10000.0;
  Rcpp::List CostMats;
  vector<arma::mat> TransportMats;
  Rcpp::List subWeights;
  arma::mat Tplan=arma::zeros(positions1.n_rows,positions1.n_cols);
  vector<arma::mat> Positions;
  vector<arma::mat> Weights;
  double beta=stepsize;
  double atsum;
  double full_change=10000000;
  int iter=0;
  int weights_iter=0;
  int pos_iter=0;
  int head_iter=0;
  int head_i;
  int L;
  int dataCount=0;
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  arma::mat positionsBary=positions1;
  for (int i=0;i<n2;i++){
    Positions.push_back((Rcpp::as<arma::mat>(positionssetR.at(i))));
    Weights.push_back((Rcpp::as<arma::mat>(weightsR.at(i))));
    TransportMats.push_back(arma::zeros(m,m));
  }
  if (fixed_support){
    positions2=Positions.at(0);
    z=arma::diagvec((positionsBary*arma::trans(positionsBary)));
    x=arma::diagvec(positions2*arma::trans(positions2));
    onesB=arma::ones(1,positions2.n_rows);
    onesA=arma::ones(positionsBary.n_rows,1);
    Cmat=(z*onesB+onesA*arma::trans(x)-(2.0*positionsBary*arma::trans(positions2)));
    
    //Weights Headstart
    if (weights_maxIter>0){
      if (headstart==TRUE){
        while (head_iter<(headstartlength*n2)){
          beta=(head_iter+2.0)/2.0;
          a=(1-(1/beta))*ahat+((1/beta)*atilde);
          head_i=index_cxx(frechet_weights,rand_cxx());
          positions2=Positions.at(head_i);
          b=Weights.at(head_i);
          duals.col(head_i)=transport_network_dual_arma(a,b,Cmat);
          s(head_i,0)=arma::sum(duals.col(head_i))/m;
          duals.col(head_i)=duals.col(head_i)-(arma::ones(m,1)*s(head_i,0));
          alpha=duals.col(head_i);
          atilde=atilde%exp(((-1)*(beta)*alpha));
          atsum=accu(atilde);
          atilde=atilde/atsum;
          ahat=(1-(1/beta))*ahat+(atilde*(1/beta));
          head_iter++;
        }
      }
      head_iter=0;
      //Weights Iterations
      
      achange(0)=10000;
      weights_iter=0;
      while( (weights_iter<weights_maxIter) && (achange(0)>thresh)){
        beta=(weights_iter+2.0)/2.0;
        a=(1-(1/beta))*ahat+((1/beta)*atilde);
#pragma omp parallel for
        for (int i=0; i<n2;i++){
          duals.col(i)=transport_network_dual_arma(a,Weights.at(i),Cmat);
          s(i,0)=arma::sum(duals.col(i))/m;
          duals.col(i)=duals.col(i)-(arma::ones(m,1)*s(i,0));
        }
        alpha=arma::sum(duals*frechet,1);
        atilde=atilde%exp(((-1)*beta*alpha));
        atsum=accu(atilde);
        atilde=atilde/atsum;
        ahat=(1-(1/beta))*ahat+(atilde*(1/beta));
        achange=arma::sum(abs((a-aold)));
        aold=a;
        weights_iter++;
      }
      
    }
    return(Rcpp::List::create(positionsBary,a,iter));
  }
  
  
  while ((full_change>thresh)&(iter<maxIter)){
    
    if (pos_maxIter>0){
      //Position Headstart
      beta=stepsize;
      if (headstart==TRUE){
        while (head_iter<(headstartlength*n2)){
          Tplan=arma::zeros(positionsBary.n_rows,positionsBary.n_cols);
          head_i=index_cxx(frechet_weights,rand_cxx());
          positions2=Positions.at(head_i);
          b=Weights.at(head_i);
          z=arma::diagvec((positionsBary*arma::trans(positionsBary)));
          x=arma::diagvec(positions2*arma::trans(positions2));
          onesB=arma::ones(1,positions2.n_rows);
          onesA=arma::ones(positionsBary.n_rows,1);
          Cmat=(z*onesB+onesA*arma::trans(x)-(2.0*positionsBary*arma::trans(positions2)));
          Tplan=((arma::diagmat(arma::ones(m,1)/a))*transport_network_primal_arma(a,b,Cmat)*positions2-positionsBary);
          positionsBary=positionsBary+(stepsize*Tplan);
          head_iter++;
        }
      }
      head_iter=0;
      //Position Update Iterations
      pchange(0)=10000;
      pos_iter=0;
      while ((pos_iter<pos_maxIter)&&(pchange(0)>thresh)){
        Tplan=arma::zeros(positionsBary.n_rows,positionsBary.n_cols);
        #pragma omp parallel for
        for (int i=0; i<n2;i++){
          arma::mat positions2=Positions.at(i);
          arma::mat z=arma::diagvec((positionsBary*arma::trans(positionsBary)));
          arma::mat x=arma::diagvec(positions2*arma::trans(positions2));
          arma::mat onesB=arma::ones(1,positions2.n_rows);
          arma::mat onesA=arma::ones(positionsBary.n_rows,1);
          arma::mat Cmat=(z*onesB+onesA*arma::trans(x)-(2.0*positionsBary*arma::trans(positions2)));
          TransportMats.at(i)=((arma::diagmat(arma::ones(m,1)/a))*transport_network_primal_arma(a,Weights.at(i),
                                              Cmat)*(Positions.at(i))-positionsBary);
        }
        for (int i=0;i<n2;i++){
          Tplan=Tplan+(frechet_weights(i,0)*TransportMats.at(i));
        }
        //Tplan=Tplan/n;
        positionsBary=positionsBary+(stepsize*Tplan);
        // Tplan=(arma::diagmat(arma::ones(m,1)/a))*(Tplan/n);
        // positionsBary=(1-beta)*positionsBary+beta*Tplan;
        pchange=arma::sum(abs((positionsBary-positionsOld)));
        positionsOld=positionsBary;
        pos_iter++;
      }
    }
    
    
    //Weights Headstart
    if (weights_maxIter>0){
      if (headstart==TRUE){
        while (head_iter<(headstartlength*n2)){
          beta=(head_iter+2)/2;
          a=(1-(1/beta))*ahat+((1/beta)*atilde);
          head_i=index_cxx(frechet_weights,rand_cxx());
          positions2=Positions.at(head_i);
          b=Weights.at(head_i);
          z=arma::diagvec((positionsBary*arma::trans(positionsBary)));
          x=arma::diagvec(positions2*arma::trans(positions2));
          onesB=arma::ones(1,positions2.n_rows);
          onesA=arma::ones(positionsBary.n_rows,1);
          Cmat=(z*onesB+onesA*arma::trans(x)-(2.0*positionsBary*arma::trans(positions2)));
          duals.col(head_i)=transport_network_dual_arma(a,b,Cmat);
          s(head_i,0)=arma::sum(duals.col(head_i))/m;
          duals.col(head_i)=duals.col(head_i)-(arma::ones(m,1)*s(head_i,0));
          alpha=duals.col(head_i);
          atilde=atilde%exp(((-1)*(beta)*alpha));
          atsum=accu(atilde);
          atilde=atilde/atsum;
          ahat=(1-(1/beta))*ahat+(atilde*(1/beta));
          head_iter++;
        }
      }
      head_iter=0;
      //Weights Iterations
      
      achange(0)=10000;
      weights_iter=0;
      while( (weights_iter<weights_maxIter) && (achange(0)>thresh)){
        beta=(weights_iter+2.0)/2.0;
        a=(1-(1/beta))*ahat+((1/beta)*atilde);
        a=a+0.0000000000000001;
        a=a/arma::accu(a);
        

          #pragma omp parallel for
          for (int i=0;i<(n2);i++){
            arma::mat positions2=Positions.at(i);
            //arma::mat positions2=positionsBary;
            arma::mat z=arma::diagvec((positionsBary*arma::trans(positionsBary)));
            arma::mat x=arma::diagvec(positions2*arma::trans(positions2));
            arma::mat onesB=arma::ones(1,positions2.n_rows);
            arma::mat onesA=arma::ones(positionsBary.n_rows,1);
            arma::mat Cmat=(z*onesB+onesA*arma::trans(x)-(2.0*positionsBary*arma::trans(positions2)));
            duals.col(i)=transport_network_dual_arma(a,Weights.at(i),Cmat);
            s(i,0)=arma::sum(duals.col(i))/m;
            duals.col(i)=duals.col(i)-(arma::ones(m,1)*s(i,0));
          }
        
        alpha=arma::sum(duals*frechet,1);
        atilde=atilde%exp(((-1)*beta*alpha));
        atsum=accu(atilde);
        atilde=atilde/atsum;
        ahat=(1-(1/beta))*ahat+(atilde*(1/beta));
        achange=arma::sum(abs((a-aold)));
        aold=a;
        weights_iter++;
      }
      
    }
    
    full_change=arma::accu(abs((a-aVeryOld)))+arma::accu(abs((positionsBary-positionsVeryOld)));
    aVeryOld=a;
    positionsVeryOld=positionsBary;
    iter++;
  }
  return(Rcpp::List::create(positionsBary,a,iter));
}
