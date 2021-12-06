#' Compute the p-Kantorovich-Rubinstein distance between two measures of possibly unequal total mass.
#' @description This function constructs the corresponding problem and solves it suing the \link[transport]{transport}-function.
#' @param A One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}. 
#' @param B One of the following: A matrix, representing an image;
#'  A file name containing an image; A \link[transport]{wpp-object}.
#' @param p A positive real number specifying the order of the Kantorovich-Rubinstein distance.
#' @param C A positive real number specifying the cost parameter of the Kantorovich-Rubinstein distance.
#' @return A list containing an entry "distance" (specifying the KR distance between the two measures) and an 
#' entry "plan" containing an optimal plan for the unbalanced optimal transport problem.
#' @examples 
#' M<-1000
#' W1<-runif(M)
#' W2<-runif(M)
#' pos1<-matrix(runif(M*2),M,2)
#' pos2<-matrix(runif(M*2),M,2)
#' wpp1<-transport::wpp(pos1,W1)
#' wpp2<-transport::wpp(pos2,W2)
#' system.time(res<-WSGeometry:::kr_dist(wpp1,wpp2,2,2))
#' @references Kantorovich-Rubinstein distance and barycenter for finitely supported measures:  Foundations and Algorithms; Florian Heinemann, Marcel Klatt, Axel Munk
#' @export
kr_dist<-function(A,B,p=2,C){
  Atype<-type_check(A)
  Btype<-type_check(B)
  A<-process_data(A,Atype)
  A<-transport::wpp(A$positions,A$weights)
  B<-process_data(B,Btype)
  B<-transport::wpp(B$positions,B$weights)
  costM<-gen_cost(B$coordinates,A$coordinates)^(p/2)
  costM<-rbind(costM,C^p/2)
  costM<-cbind(costM,C^p/2)
  costM[dim(costM)[1],dim(costM)[2]]<-0
  a<-c(A$mass,sum(B$mass))
  b<-c(B$mass,sum(A$mass))
  res<-(transport::transport(a,b,costM,method="networkflow",fullreturn=TRUE))
  plan<-res$primal
  plan<-plan[1:(length(a)-1),1:(length(b)-1)]
  return(list(distance=res$cost^(1/p),plan=plan))
}