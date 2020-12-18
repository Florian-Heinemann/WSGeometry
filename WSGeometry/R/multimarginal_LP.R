#' Solve the multimarginal optimal transport problem by linear programming
#' @description Solves the N-fold multimarginal optimal transport problem between N specified measures and a specified cost.  This is essentially a convenient wrapper function that builds and solves the corresponding linear program.
#' @param weights A list of vectors specifying the weights of the marginal distributions. These vectors do not need to be of the same size.
#' @param costA An array where the entry (i1,i2,...,iN) specifies the value of the cost functional for the point i1 in the first measure, i2 in the second measure
#' and so on.
#' @return A list with two entries. The first entry contains the optimal multicoupling in array form, and the second entry contains 
#' the cost of the optimal solution. 
#' @examples
#' W<-list(rep(1,10),rep(1,10),rep(1,10))
#' C<-array(runif(10^3),c(10,10,10))
#' MM<-multi_marginal(W,C)
#' @export

multi_marginal<-function(weights,costA){
  D<-dim(costA)
  const<-gen_constraints_multi(D)
  costVec<-build_MM_costvec(costA,const)
  rhs<-unlist(weights)
  out<-Rsymphony::Rsymphony_solve_LP(obj=costVec,mat=const,dir=rep("==",sum(D)),rhs=rep(1,sum(D)),max=FALSE)
  optMMCoupling<-build_MMCoupling(out,const,D)
  return(list(MMCoupling=optMMCoupling,cost=out$objval))
}

