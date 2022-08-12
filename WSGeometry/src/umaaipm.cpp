#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <ctime>
#ifdef _OPENMP
# include <omp.h>
#endif
using namespace Rcpp;
using namespace std;
using namespace arma;
//[[Rcpp::plugins(openmp)]]
// [[Rcpp::depends("RcppArmadillo")]]


arma::vec USLRM_intern(vector<arma::mat> A1,arma::mat Y,arma::vec g,arma::sp_mat U,int N,int m){
  arma::mat Asum = zeros(Y.n_rows,Y.n_cols) ;
#pragma omp parallel for
  for (int i=0;i<N;i++){
    A1.at(i)=inv(A1.at(i));
  }
  for (int i=0;i<N;i++){
    Asum+=A1.at(i);
  }
  arma::vec x1=zeros(g.n_elem);
  
  for (int i=0;i<N;i++){
    x1.subvec(((i)*(m-1)),((i+1)*(m-1)-1))=(A1.at(i)*g.subvec(((i)*(m-1)),((i+1)*(m-1)-1)));
  }  
  arma::vec x2=U.t()*x1;
  int L=U.n_cols;
  arma::vec x3=zeros(L);
  x3.subvec((L-m+1),L-1)=solve(inv(Y)+Asum,x2.subvec((L-m+1),L-1),solve_opts::likely_sympd);
  arma::vec x4=U*x3;
  for (int i=0;i<N;i++){
    x4.subvec(((i)*(m-1)),((i+1)*(m-1)-1))=(A1.at(i)*x4.subvec(((i)*(m-1)),((i+1)*(m-1)-1)));
  }  
  return(x1-x4);
}


arma::vec UDLRM_intern(vector<arma::sp_mat> B1,vector<arma::mat> B2,vector<arma::sp_mat> B3inv,arma::mat Y,const arma::vec sizes_csum,const arma::vec g,const int N,const int m,const arma::sp_mat U) {
  vector<arma::mat> Binv;
  arma::mat Asum = zeros(Y.n_rows,Y.n_cols) ;
  arma::mat B3sum = zeros(Y.n_rows,Y.n_cols) ;
  int diff=0;
  for (int i=0;i<N;i++){
    diff=sizes_csum(i+1)-sizes_csum(i);
    Binv.push_back(zeros(diff,diff));
  }
#pragma omp parallel for
  for (int i=0;i<N;i++){
    Binv.at(i)=(inv(B1.at(i)-(B2.at(i)*B3inv.at(i)*B2.at(i).t())));
  }
  
#pragma omp parallel for 
  for (int i=0;i<N;i++){
    Asum=Asum+(B3inv.at(i)*B2.at(i).t()*Binv.at(i)*B2.at(i)*B3inv.at(i));
    B3sum=B3sum+B3inv.at(i);
  }
  Asum+=B3sum;
  arma::vec x1=g;
  for (int i=0; i<N;i++){
    x1.subvec(((i)*(m-1)),((i+1)*(m-1)-1))=(B3inv.at(i)*(B2.at(i).t()*(Binv.at(i)*(B2.at(i)*(B3inv.at(i)*x1.subvec(((i)*(m-1)),((i+1)*(m-1)-1)))))))+(B3inv.at(i)* x1.subvec(((i)*(m-1)),((i+1)*(m-1)-1)));
  }
  arma::vec x2=U.t()*x1;
  int L=U.n_cols;
  arma::vec x3=zeros(L);
  x3.subvec((L-m+1),L-1)=solve(inv(Y)+Asum,x2.subvec((L-m+1),L-1));
  arma::vec x4=U*x3;
  for (int i=0; i<N;i++){
    x4.subvec(((i)*(m-1)),((i+1)*(m-1)-1))=(B3inv.at(i)*(B2.at(i).t()*(Binv.at(i)*(B2.at(i)*(B3inv.at(i)*x4.subvec(((i)*(m-1)),((i+1)*(m-1)-1)))))))+(B3inv.at(i)* x4.subvec(((i)*(m-1)),((i+1)*(m-1)-1)));
  }
  
  return(x1-x4);
}

arma::mat Udist_mat(const arma::mat X,const arma::mat Y,const double C, const double p){
  int n=X.n_cols;
  int m=Y.n_cols;
  arma::mat x(n,1,fill::zeros);
  arma::mat y(m,1,fill::zeros);
  for (int i=0;i<n;i++){
    x(i)=accu(X.col(i)%X.col(i));
  }
  for (int i=0;i<m;i++){
    y(i)=accu(Y.col(i)%Y.col(i));
  }
  arma::mat distMat=(x*ones(1,m))+(ones(n,1)*trans(y))-2.0*(trans(X)*Y);
  distMat.elem(find(distMat>pow(10,12))).fill(pow(C,p)/2);
  return(distMat);
}


double Usqnorm(arma::vec x){
  return(accu(x%x));
}




//[[Rcpp::export]]
Rcpp::List umaaipm_free_cpp(arma::vec p, arma::vec s, arma::vec x,const arma::vec b,arma::vec costvec,arma::mat support,arma::mat fullsupport, const arma::sp_mat constMat,const int N,const int m,const int M,const arma::vec sizes,const arma::vec sizescsum,const int nr, const int nc,const arma::sp_mat U,const double C,const int maxIter,const int outer_maxIter,const double thresh,const int threads){
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  int iter=0;
  double m2=m;
  const int Mm=M*m;
  double bc = 1+std::max(sqrt(Usqnorm(costvec)), sqrt(Usqnorm(b)));
  arma::vec Rc = x%s;
  double mu=mean(Rc);
  arma::vec Rd=constMat.t()*p+s-costvec;
  arma::vec Rp=constMat*x-b;
  double relResidual=(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
  double rel_gap = 1e15;
  const double maxDiag=5.e+14;
  const int data_dim=support.n_rows;
  arma::vec d=x/s;
  arma::vec t1;
  arma::vec t2;
  arma::mat dpile = zeros(m,M);
  arma::vec B1diag=zeros(M);
  vector<arma::sp_mat> B1;
  vector<arma::sp_mat> B1inv;
  vector<arma::mat> B2;
  vector<arma::sp_mat> B3;
  vector<arma::sp_mat> B3inv;
  vector<arma::mat> A1;
  vector<arma::mat> T;
  arma::vec y=zeros(m-1);
  arma::mat Y=zeros(m-1,m-1);
  arma::vec alpha;
  arma::vec xx;
  arma::vec dp=zeros(p.n_elem);
  arma::vec dx=zeros(x.n_elem);
  arma::vec ds=zeros(s.n_elem);
  int L=x.n_elem;
  double cc;
  double eta;
  double sigma;
  double muaff;
  double etaMin=0.95;
  double alphax;
  double alphas;
  arma::mat support2=zeros(data_dim,m);
  bool largesupp=FALSE;
  if ((m^2)>(4*(sum(sizes%sizes)))){
    largesupp=TRUE;
  }
  arma::mat X(m,M,fill::zeros);
  arma::vec Xrs(m,fill::zeros);
  int outer_iter=0;
  
  while (outer_iter<=outer_maxIter){
    iter=0;
    Rc = x%s;
    Rd=constMat.t()*p+s-costvec;
    Rp=constMat*x-b;
    mu=mean(Rc);
    rel_gap = 1e15;
    relResidual=(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
    for (int i=0;i<N;i++){
      x.subvec(m*sizescsum(i),m*sizescsum(i+1)-1)=((1/m2)*ones(m,1)*trans(b.subvec(sizescsum(i),sizescsum(i+1)-1))).as_col();
    }
    x.subvec(Mm,m*(M+1)-1)=ones(m)/m2;
    p.subvec(0,M-1)=ones(M)*(-1);
    p.subvec(M,M+(N*(m-1))-1)=zeros(N*(m-1));
    p(M+N*(m-1))=-1;
    s=trans(trans(costvec)-trans(p)*constMat);
    rel_gap=999;
    
    
    d=x/s;
    while((iter<maxIter) & ((mu>thresh)| (rel_gap>thresh)| (relResidual>thresh))){

      rel_gap=(accu(costvec%x) - accu(b%p))/(abs(accu(b%p))+abs(accu(costvec%x))+1);
      //optval=accu(costvec%x)/N2;
      
      d=x/s; 
      d.elem(find(d>=maxDiag)).fill(maxDiag); 
      t1=x%Rd-Rc; 
      t2=-(Rp+constMat*(t1/s)); 
      //start first solve
      dpile=reshape(d.subvec(0,(Mm-1)),m,M); 
      B1diag=(sum(dpile,0)).t(); 
      dpile.shed_row(0); 
      if (iter==0){
        for (int i=0;i<N;i++){
          B2.push_back((dpile.cols(sizescsum(i),sizescsum(i+1)-1)).t());
          B1.push_back(speye(sizes(i),sizes(i)));
          B1inv.push_back(speye(sizes(i),sizes(i)));
          B1.at(i).diag()=(B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          B1inv.at(i).diag()=(1.0/B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          T.push_back(B2.at(i).t()*B1inv.at(i)); 
          B3.push_back(speye(m-1,m-1)); 
          B3inv.push_back(speye(m-1,m-1)); 
          B3.at(i).diag()=sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          B3inv.at(i).diag()=1.0/sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          A1.push_back((B3.at(i)-(T.at(i)*B2.at(i))));
        }
      }
      else{
#pragma omp parallel for
        for (int i=0;i<N;i++){
          B2.at(i)=(dpile.cols(sizescsum(i),sizescsum(i+1)-1)).t();
          B1.at(i).diag()=(B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          B1inv.at(i).diag()=(1.0/B1diag.subvec(sizescsum(i),sizescsum(i+1)-1));
          T.at(i)=B2.at(i).t()*B1inv.at(i);
          B3.at(i).diag()=sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          B3inv.at(i).diag()=1.0/sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          A1.at(i)=(B3.at(i)-(T.at(i)*B2.at(i)));
        }
      }
      cc=accu(d.subvec(nc-m,nc-1));
      y=d.subvec(nc-m+1,nc-1);
      Y=diagmat(y)-(y*y.t())/cc;
      alpha=kron(-1*ones(N),y);
      xx=t2;
      for (int i=0;i<N;i++){
        xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))= xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))-T.at(i)*xx.subvec(sizescsum(i),sizescsum(i+1)-1);
      }  
      xx.subvec(M,nr-2)= xx.subvec(M,nr-2)-(alpha*xx(nr-1)/cc);
      xx.subvec(0,M-1)=xx.subvec(0,M-1)/B1diag;
      xx(nr-1)/=cc;
      if ((rel_gap>=(1e-3))&& (largesupp)){
        xx.subvec(M,nr-2)=UDLRM_intern(B1,B2,B3inv,Y,sizescsum, xx.subvec(M,nr-2),N,m,U);
      }
      else{   
        xx.subvec(M,nr-2)=USLRM_intern(A1,Y,xx.subvec(M,nr-2),U,N,m);
      }
      
      
      xx(nr-1)=xx(nr-1)-(accu(xx.subvec(M,nr-2)%alpha)/cc);
      for (int i=0;i<N;i++){
        xx.subvec(sizescsum(i),sizescsum(i+1)-1)= xx.subvec(sizescsum(i),sizescsum(i+1)-1)-(T.at(i).t()*xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1)));
      } 
      dp=xx;
      
      //end first solve
      
      
      
      dx=((constMat.t()*dp)%x+t1)/s;
      ds=-(s%dx+Rc)/x;
      eta=std::max(etaMin,1-mu);
      alphax=-1.0/std::min(arma::min(dx/x),-1.0);
      alphas=-1.0/std::min(arma::min(ds/s),-1.0);
      muaff=accu((x+alphax*dx)%(s+alphas*ds))/(nc*1.0);
      sigma=pow((muaff/mu),3);
      
      Rc=(Rc+(dx%ds))-(sigma*mu);
      t1=x%Rd-Rc;
      t2=-(Rp+constMat*(t1/s));
      
      
      //start second solve
      
      xx=t2;
      for (int i=0;i<N;i++){
        xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))= xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))-T.at(i)*xx.subvec(sizescsum(i),sizescsum(i+1)-1);
      }  
      xx.subvec(M,nr-2)= xx.subvec(M,nr-2)-(alpha*xx(nr-1)/cc);
      xx.subvec(0,M-1)=xx.subvec(0,M-1)/B1diag;
      xx(nr-1)/=cc;
      if ((rel_gap>=(1e-3))&& (largesupp)){
        xx.subvec(M,nr-2)=UDLRM_intern(B1,B2,B3inv,Y,sizescsum, xx.subvec(M,nr-2),N,m,U);
      }
      else{   
        xx.subvec(M,nr-2)=USLRM_intern(A1,Y,xx.subvec(M,nr-2),U,N,m);
      }
      xx(nr-1)=xx(nr-1)-(accu(xx.subvec(M,nr-2)%alpha)/cc);
      for (int i=0;i<N;i++){
        xx.subvec(sizescsum(i),sizescsum(i+1)-1)= xx.subvec(sizescsum(i),sizescsum(i+1)-1)-(T.at(i).t()*xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1)));
      } 
      dp=xx;
      
      //end second solve
      
      
      dx=((constMat.t()*dp)%x+t1)/s;
      ds=-(s%dx+Rc)/x;
      
      alphax=-1.0/std::min(arma::min(dx/x),-1.0);
      alphax=std::min(1.0,eta*alphax);
      alphas=-1.0/std::min(arma::min(ds/s),-1.0);
      alphas=std::min(1.0,eta*alphas);
      //alpha = std::min(alphax, alphas);
      
      x+=alphax*dx;
      s+=alphas*ds;
      p+=alphas*dp;
      
      
      Rd=constMat.t()*p+s-costvec;
      Rp=constMat*x-b;
      Rc=x%s;
      mu=mean(Rc);
      relResidual=sqrt(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
      
      
      iter++;
    }
    
    //update positions 
    X=reshape(x.subvec(0,(L-m)-1),m,M);
    X.row(m-1).fill(0.0);
    for (int i=0;i<N;i++){
      X.col(sizescsum(i+1)-1).fill(0.0);
    }
    arma::mat tmp(data_dim,m,fill::zeros);
    Xrs=(sum(X,1));
    for (int i=0;i<data_dim;i++){
      tmp.row(i)=trans(Xrs);
    }
    support2=(fullsupport*X.t());
    support2.elem(find(tmp>0))/=tmp.elem(find(tmp>0));
    support.elem(find(tmp>0))=support2.elem(find(tmp>0));
    costvec.subvec(0,M*m-1)=vectorise(Udist_mat(support,fullsupport,C,2.0));

    outer_iter++;
  }
  
  
  
  return(Rcpp::List::create(support,x,X,Xrs));
  
}





//[[Rcpp::export]]
Rcpp::List umaaipm_fixed_cpp(arma::vec p, arma::vec s, arma::vec x,const arma::vec b,arma::vec costvec,arma::mat support,arma::mat fullsupport, const arma::sp_mat constMat,const int N,const int m,const int M,const arma::vec sizes,const arma::vec sizescsum,const int nr, const int nc,const arma::sp_mat U,const double C,const int maxIter,const double thresh,const int threads){
#ifdef _OPENMP
  omp_set_num_threads(threads);
#endif
  int iter=0;
  double m2=m;
  const int Mm=M*m;
  double bc = 1+std::max(sqrt(Usqnorm(costvec)), sqrt(Usqnorm(b)));
  arma::vec Rc = x%s;
  double mu=mean(Rc);
  arma::vec Rd=constMat.t()*p+s-costvec;
  arma::vec Rp=constMat*x-b;
  double relResidual=(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
  double rel_gap = 1e15;
  const double maxDiag=5.e+14;
  const int data_dim=support.n_rows;
  arma::vec d=x/s;
  arma::vec t1;
  arma::vec t2;
  arma::mat dpile = zeros(m,M);
  arma::vec B1diag=zeros(M);
  vector<arma::sp_mat> B1;
  vector<arma::sp_mat> B1inv;
  vector<arma::mat> B2;
  vector<arma::sp_mat> B3;
  vector<arma::sp_mat> B3inv;
  vector<arma::mat> A1;
  vector<arma::mat> T;
  arma::vec y=zeros(m-1);
  arma::mat Y=zeros(m-1,m-1);
  arma::vec alpha;
  arma::vec xx;
  arma::vec dp=zeros(p.n_elem);
  arma::vec dx=zeros(x.n_elem);
  arma::vec ds=zeros(s.n_elem);
  double cc;
  double eta;
  double sigma;
  double muaff;
  double etaMin=0.95;
  double alphax;
  double alphas;
  arma::mat support2=zeros(data_dim,m);
  bool largesupp=FALSE;
  if ((m^2)>(4*(sum(sizes%sizes)))){
    largesupp=TRUE;
  }
  arma::mat X(m,M,fill::zeros);
  arma::vec Xrs(m,fill::zeros);
    iter=0;
    Rc = x%s;
    Rd=constMat.t()*p+s-costvec;
    Rp=constMat*x-b;
    mu=mean(Rc);
    rel_gap = 1e15;
    relResidual=(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
    for (int i=0;i<N;i++){
      x.subvec(m*sizescsum(i),m*sizescsum(i+1)-1)=((1/m2)*ones(m,1)*trans(b.subvec(sizescsum(i),sizescsum(i+1)-1))).as_col();
    }
    x.subvec(Mm,m*(M+1)-1)=ones(m)/m2;
    p.subvec(0,M-1)=ones(M)*(-1);
    p.subvec(M,M+(N*(m-1))-1)=zeros(N*(m-1));
    p(M+N*(m-1))=-1;
    s=trans(trans(costvec)-trans(p)*constMat);
    rel_gap=999;
    
    
    d=x/s;
    while((iter<maxIter) & ((mu>thresh)| (rel_gap>thresh)| (relResidual>thresh))){

      rel_gap=(accu(costvec%x) - accu(b%p))/(abs(accu(b%p))+abs(accu(costvec%x))+1);
      //optval=accu(costvec%x)/N2;
      
      d=x/s; 
      d.elem(find(d>=maxDiag)).fill(maxDiag); 
      t1=x%Rd-Rc; 
      t2=-(Rp+constMat*(t1/s)); 
      //start first solve
      dpile=reshape(d.subvec(0,(Mm-1)),m,M); 
      B1diag=(sum(dpile,0)).t(); 
      dpile.shed_row(0); 
      if (iter==0){
        for (int i=0;i<N;i++){
          B2.push_back((dpile.cols(sizescsum(i),sizescsum(i+1)-1)).t()); 
          B1.push_back(speye(sizes(i),sizes(i))); 
          B1inv.push_back(speye(sizes(i),sizes(i)));
          B1.at(i).diag()=(B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          B1inv.at(i).diag()=(1.0/B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          T.push_back(B2.at(i).t()*B1inv.at(i)); 
          B3.push_back(speye(m-1,m-1)); 
          B3inv.push_back(speye(m-1,m-1)); 
          B3.at(i).diag()=sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          B3inv.at(i).diag()=1.0/sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          A1.push_back((B3.at(i)-(T.at(i)*B2.at(i))));
        }
      }
      else{
#pragma omp parallel for
        for (int i=0;i<N;i++){
          B2.at(i)=(dpile.cols(sizescsum(i),sizescsum(i+1)-1)).t();
          B1.at(i).diag()=(B1diag.subvec(sizescsum(i),sizescsum(i+1)-1)); 
          B1inv.at(i).diag()=(1.0/B1diag.subvec(sizescsum(i),sizescsum(i+1)-1));
          T.at(i)=B2.at(i).t()*B1inv.at(i);
          B3.at(i).diag()=sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          B3inv.at(i).diag()=1.0/sum(dpile.cols(sizescsum(i),sizescsum(i+1)-1),1); 
          A1.at(i)=(B3.at(i)-(T.at(i)*B2.at(i)));
        }
      }
      cc=accu(d.subvec(nc-m,nc-1));
      y=d.subvec(nc-m+1,nc-1);
      Y=diagmat(y)-(y*y.t())/cc;
      alpha=kron(-1*ones(N),y);
      xx=t2;
      for (int i=0;i<N;i++){
        xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))= xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))-T.at(i)*xx.subvec(sizescsum(i),sizescsum(i+1)-1);
      }  
      xx.subvec(M,nr-2)= xx.subvec(M,nr-2)-(alpha*xx(nr-1)/cc);
      xx.subvec(0,M-1)=xx.subvec(0,M-1)/B1diag;
      xx(nr-1)/=cc;
      if ((rel_gap>=(1e-3))&& (largesupp)){
        xx.subvec(M,nr-2)=UDLRM_intern(B1,B2,B3inv,Y,sizescsum, xx.subvec(M,nr-2),N,m,U);
      }
      else{   
        xx.subvec(M,nr-2)=USLRM_intern(A1,Y,xx.subvec(M,nr-2),U,N,m);
      }
      
      
      xx(nr-1)=xx(nr-1)-(accu(xx.subvec(M,nr-2)%alpha)/cc);
      for (int i=0;i<N;i++){
        xx.subvec(sizescsum(i),sizescsum(i+1)-1)= xx.subvec(sizescsum(i),sizescsum(i+1)-1)-(T.at(i).t()*xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1)));
      } 
      dp=xx;
      
      //end first solve
      
      

      dx=((constMat.t()*dp)%x+t1)/s;
      ds=-(s%dx+Rc)/x;
      eta=std::max(etaMin,1-mu);
      alphax=-1.0/std::min(arma::min(dx/x),-1.0);
      alphas=-1.0/std::min(arma::min(ds/s),-1.0);
      muaff=accu((x+alphax*dx)%(s+alphas*ds))/(nc*1.0);
      sigma=pow((muaff/mu),3);
      
      Rc=(Rc+(dx%ds))-(sigma*mu);
      t1=x%Rd-Rc;
      t2=-(Rp+constMat*(t1/s));

      //start second solve
      
      xx=t2;
      for (int i=0;i<N;i++){
        xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))= xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1))-T.at(i)*xx.subvec(sizescsum(i),sizescsum(i+1)-1);
      }  
      xx.subvec(M,nr-2)= xx.subvec(M,nr-2)-(alpha*xx(nr-1)/cc);
      xx.subvec(0,M-1)=xx.subvec(0,M-1)/B1diag;
      xx(nr-1)/=cc;
      if ((rel_gap>=(1e-3))&& (largesupp)){
        xx.subvec(M,nr-2)=UDLRM_intern(B1,B2,B3inv,Y,sizescsum, xx.subvec(M,nr-2),N,m,U);
      }
      else{   
        xx.subvec(M,nr-2)=USLRM_intern(A1,Y,xx.subvec(M,nr-2),U,N,m);
      }
      xx(nr-1)=xx(nr-1)-(accu(xx.subvec(M,nr-2)%alpha)/cc);
      for (int i=0;i<N;i++){
        xx.subvec(sizescsum(i),sizescsum(i+1)-1)= xx.subvec(sizescsum(i),sizescsum(i+1)-1)-(T.at(i).t()*xx.subvec((M+(i)*(m-1)),M-1+(m-1)*(i+1)));
      } 
      dp=xx;
      
      //end second solve

      
      dx=((constMat.t()*dp)%x+t1)/s;
      ds=-(s%dx+Rc)/x;
      
      alphax=-1.0/std::min(arma::min(dx/x),-1.0);
      alphax=std::min(1.0,eta*alphax);
      alphas=-1.0/std::min(arma::min(ds/s),-1.0);
      alphas=std::min(1.0,eta*alphas);
      //alpha = std::min(alphax, alphas);
      
      x+=alphax*dx;
      s+=alphas*ds;
      p+=alphas*dp;
      
      
      Rd=constMat.t()*p+s-costvec;
      Rp=constMat*x-b;
      Rc=x%s;
      mu=mean(Rc);
      relResidual=sqrt(Usqnorm(Rd)+Usqnorm(Rp)+Usqnorm(Rc))/bc;
      
      
      iter++;
    }
    

  return(Rcpp::List::create(support,x,X,Xrs));
  
}
