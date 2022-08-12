maaipm_fixed_wrap<-function(pos.list,weights.list,support,maxIter=200,startvec=NULL,thresh=10^-5,threads=1){
  t<-Matrix::t
  d<-dim(pos.list[[1]])[2]
  N<-length(pos.list)
  M<-dim(support)[2]
  sizes<-unlist(lapply(weights.list,length))
  #constraints
  const<-build_const(M,sizes)
  gc()
  #build rhs
  b<-c(unlist(weights.list),rep(0,M*N),1)
  b<-b[-seq(sum(sizes)+1,length(b)-1,M)]
  #build_cost
  costvec<-c(newpgen_cost(support,Reduce("cbind",pos.list),10^10),rep(0,M))
  N<-length(sizes)
  sizes_csum<-c(0,cumsum(sizes))
  m<-M
  M<-sum(sizes)
  n_row<-M+N*(m-1)+1
  n_col<-(M+1)*m
  Mm<-M*m
  if(is.null(startvec)){
    x<-rep(0,(M+1)*m)
    for (i in 1:N){
      pi<-(1/m)*rep(1,m)%*%t(b[(sizes_csum[i]+1):(sizes_csum[i+1])])
      x[(m*sizes_csum[i]+1):(m*sizes_csum[i+1])]<-pi
    }
    x[(Mm+1):(m*(M+1))]<-rep(1/m,m)
  }
  else{
    x<-startvec
  }
  p<-c(rep(-1,M),rep(0,N*(m-1)),-1)
  s<-as.vector(t(costvec)-t(p)%*%const)
  gc() 
  x<-maaipm_fixed_cpp(p,s,x,b,costvec, const,N,m,M,sizes,sizes_csum,n_row, n_col,build_U(N,m),maxIter,thresh,threads)[[1]]
  L<-length(x)
 # return(x)
 return(x[(L-m+1):L])
}

