#' Exact computation of 2-Wasserstein barycenters in R^d using linear programming
#' @description This function solves the 2-Wasserstein barycenter problem of N finitely supported input
#' measures explicitly by solving the corresponding linear program.
#' @param pos.list A list of Mxd matrices, specifying the positions of the data measures.
#' @param weights.list A list of vectors with non-negative entries and identical total sum specifying
#' the weights of the data measures.
#' @param frechet.weights A vector of positive entries summing to one specifying the weights of each data measure
#' in the Fréchet functional.
#' @return A list with two entries. The first entry contains the positions of
#'  the computed barycenter and the second entry contains the corresponding weights.
#' @examples
#' pos.list<-vector("list",4)
#' weights.list<-vector("list",4)
#' wpp.list<-vector("list",4)
#' pos.list[[1]]<-matrix(c(0,0,1,1,1,0,0,1),nrow=4,ncol=2)/10
#' pos.list[[2]]<-matrix(c(9,9,10,10,10,9,9,10),nrow=4,ncol=2)/10
#' pos.list[[3]]<-matrix(c(9,9,10,10,1,0,0,1),nrow=4,ncol=2)/10
#' pos.list[[4]]<-matrix(c(0,0,1,1,10,9,9,10),nrow=4,ncol=2)/10
#' plot(0, 0, xlab = "", ylab = "", type = "n", xlim = c(0, 1), ylim = c(0, 1))
#' for(i in 1:4)
#'   points(pos.list[[i]][,1], pos.list[[i]][,2], col = i)
#' weights.list[[1]]<-rep(1/4,4)
#' weights.list[[2]]<-rep(1/4,4)
#' weights.list[[3]]<-rep(1/4,4)
#' weights.list[[4]]<-rep(1/4,4)
#' bary<-barycenter_lp(pos.list,weights.list)
#' wpp.list[[1]]<-wpp(pos.list[[1]],weights.list[[1]])
#' wpp.list[[2]]<-wpp(pos.list[[2]],weights.list[[2]])
#' wpp.list[[3]]<-wpp(pos.list[[3]],weights.list[[3]])
#' wpp.list[[4]]<-wpp(pos.list[[4]],weights.list[[4]])
#' bary.wpp<-wpp(bary$positions,bary$weights)
#' wpp.list[[5]]<-bary.wpp
#' multiwppplot(wpp.list,pch=c(21,21,21,21,24),bg=c("blue","red","darkgreen","purple","orange")
#' ,col=rep("black",5),cex=3,clean=TRUE)
#' @references E Anderes, S Borgwardt, and J Miller (2016). Discrete Wasserstein barycenters: 
#' Optimal transport for discrete data. Mathematical Methods of Operations Research, 84(2):389-409.\cr
#' S Borgwardt and S Patterson (2020). Improved linear programs for discrete barycenters.
#' Informs Journal on Optimization 2(1):14-33.
#' @export

barycenter_lp<-function(pos.list,weights.list,frechet.weights=NULL){
  if (!requireNamespace("Rsymphony", quietly = TRUE)) {
    warning("Package Rsymphony not detected: Please install Rsymphony for optimal performance 
            if you are not using a Mac. Alternatively, install lpSolve.")
    Rsym<-FALSE
    if (!requireNamespace("lpSolve", quietly = TRUE)) {
      stop("Neither Rsymphony nor lpSolve detected. Please install at least one of the two to utilise this method.")
    }
  }
  else{
    Rsym<-TRUE
  }
  d<-dim(pos.list[[1]])[2]
  sizes<-unlist(lapply(weights.list,length))
  const<-(gen_constraints_multi(sizes))
  bol_const<-const>0
  rhs<-unlist(weights.list)
  N<-length(pos.list)
  if (is.null(frechet.weights)){
    frechet.weights<-rep(1/N,N)
  }
  pos.full<-matrix(0,0,d)
  for (i in 1:N){
    pos.full<-rbind(pos.full,pos.list[[i]])
  }
  M<-prod(sizes)
  cost<-NULL
  for (k in 1:M){
    mat<-pos.full[bol_const[,k],]
    cost<-c(cost,euclid_bary_func_w(mat,frechet.weights))
  }
  if (Rsym){
    out<-Rsymphony::Rsymphony_solve_LP(obj=cost,mat=const,dir=rep("==",length(rhs)),rhs=rhs,max=FALSE)
  }
  else{
    out<-lpSolve::lp("min",cost,const,rep("==",length(rhs)),rhs)
  }

  inds<-out$solution>0
  A.sub<-const[,inds]
  out.pos<-matrix(0,0,d)
  for (k in 1:(dim(A.sub)[2])){
    pos<-A.sub[,k]>0
    pos.sub<-diag(frechet.weights)%*%(pos.full[pos,])
    out.pos<-rbind(out.pos,colSums(pos.sub))
  }
  weights.out<-out$sol[inds]
  return(list(positions=out.pos,weights=weights.out))
}

