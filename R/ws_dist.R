#' Compute the p-Wasserstein distance between two measures
#' @description This is essentially a wrapper function of \link[transport]{transport}. It has the advantage of allowing 
#' more general input objects, such as images or matrices, without the user having to manually convert these objects. 
#' @param A One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}. 
#' @param B One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}.
#' @param p A positive real number specifying the power of the Wasserstein distance.
#' @param sampling A boolean specifying whether a stochastic approximation (Sommerfeld et al., 2019) should be used to approximate the distance. 
#' @param S A positive integer specifying the number of samples drawn in the stochastic approximation. 
#' @param R The number of repetitions averaged over in the stochastic approximation. 
#' @return A number specifying the computed p-Wasserstein distance.
#' @examples 
#' \dontrun{
#' P1<-transport::random32a$mass
#' P2<-transport::random32b$mass
#' P1<-P1/sum(P1)
#' P2<-P2/sum(P2)
#' ws_dist(P1,P2)
#' }
#' @references M Sommerfeld, J Schrieber, Y Zemel, and A Munk (2019).
#'Optimal transport: Fast probabilistic approximations with exact solvers.  Journal of Machine Learning Research 20(105):1--23.
#' @export
ws_dist<-function(A,B,p=2,sampling=FALSE,S=NULL,R=NULL){
  Atype<-type_check(A)
  Btype<-type_check(B)
  A<-process_data(A,Atype)
  A<-transport::wpp(A$positions,A$weights)
  B<-process_data(B,Btype)
  B<-transport::wpp(B$positions,B$weights)
  if (sampling){
    res<-0
    for (r in 1:R){
      A.sub<-data_sample(list(A$coordinates),list(A$mass),S)
      A.sub<-transport::wpp(A.sub$positions[[1]],A.sub$weights[[1]])
      B.sub<-data_sample(list(B$coordinates),list(B$mass),S)
      B.sub<-transport::wpp(B.sub$positions[[1]],B.sub$weights[[1]])
      res<-res+(transport::transport(A.sub,B.sub,p,method="networkflow",fullreturn=TRUE)$cost)^(1/p)
    }
    res<-res/R
  }
  else{
    res<-(transport::transport(A,B,p,method="networkflow",fullreturn=TRUE)$cost)^(1/p)
  }
  return(res)
}