#' Compute Wasserstein geodesics
#' @description Computes the geodesic between two measures P1 and P2 in arbitrary dimensions at given timepoints.
#' @param P1 One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}.
#' @param P2 One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}.
#' @param p A real number >=1 specifying the exponent of the Wasserstein distance.
#' @param steps A vector of numbers in [0,1] describing the time points at which the geodesic should be evaluated.
#' @return A list of the same length as steps where each element is a \link[transport]{wpp-object} describing the value of the geodesic at the 
#' corresponding times in steps. 
#' @examples 
#' \dontrun{
#' U<-runif(20)
#' U<-U/sum(U)
#' pos<-matrix(runif(2*20),nrow=20,ncol=2)
#' P1<-transport::wpp(pos,U)
#' U<-runif(20)
#' U<-U/sum(U)
#' pos<-matrix(runif(2*20),nrow=20,ncol=2)
#' P2<-transport::wpp(pos,U)
#' geodesic<-geodesic_pos(P1,P2,p=2,seq(0,1,0.1))
#' plotGeodesic(geodesic,File="GeodesicR2")
#' }
#' @export
geodesic_pos<-function(P1,P2,p=2,steps){
  P1type<-type_check(P1)
  P2type<-type_check(P2)
  P1<-process_data(P1,P1type,return_type="wpp")
  P2<-process_data(P2,P2type,return_type="wpp")
  plan<-transport::transport(P1,P2,p)
  L<-length(plan$from)
  d<-P1$dimension
  start<-matrix(0,L,d)
  V<-matrix(0,L,d)
  W<-rep(0,L)
  for (k in 1:L){
    f<-plan$from[k]
    t<-plan$to[k]
    start[k,]<-P1$coordinates[f,]
    V[k,]<-P2$coordinates[t,]-P1$coordinates[f,]
    W[k]<-plan$mass[k]
  }
  N<-length(steps)
  out.list<-vector("list",N)
  for (i in 1:N){
    out.list[[i]]<-transport::wpp(start+(V*steps[i]),W)
  }
  return(out.list)
}
