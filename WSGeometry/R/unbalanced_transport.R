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
#' M<-2000
#' W1<-runif(M)
#' W2<-runif(M)
#' pos1<-matrix(runif(M*2),M,2)
#' pos2<-matrix(runif(M*2),M,2)
#' wpp1<-wpp(pos1,W1)
#' wpp2<-wpp(pos2,W2)
#' system.time(res<-kr_dist(wpp1,wpp2,2,2))
#' @references Kantorovich-Rubinstein distance and barycenter for finitely supported measures:  Foundations and Algorithms; Florian Heinemann, Marcel Klatt, Axel Munk; https://arxiv.org/pdf/2112.03581.pdf.
#' @export
kr_dist<-function(A,B,p=2,C){
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("The package transport is required for this method. Please install it to use this function.")
  }
  Atype<-type_check(A)
  Btype<-type_check(B)
  A<-process_data_unb(A,Atype)
  A<-wpp(A$positions,A$weights)
  B<-process_data_unb(B,Btype)
  B<-wpp(B$positions,B$weights)
  if ((sum(A$mass)==0)&(sum(B$mass)==0)){
    return(list(distance=0,plan=matrix(0,A$N,B$N)))
  }
  if ((sum(A$mass)>0)&(sum(B$mass)==0)){
    return(list(distance=0.5*C^p*sum(A$mass),plan=matrix(0,A$N,B$N)))
  }
  if ((sum(A$mass)==0)&(sum(B$mass)>0)){
    return(list(distance=0.5*C^p*sum(B$mass),plan=matrix(0,A$N,B$N)))
  }
  
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