#' Bin data onto a grid.
#' @description Bin data onto a equidistant grid in [0,1]^2. 
#' @param data.pos A Mx2 matrix specifying the positions of the data measure.
#' @param data.weights A list of vectors of the same size as the number of rows in data.pos.
#' All entries in the vector must be non-negative and the entries in the vector must sum to one.
#' @param gridsize A vector of two integers specifying the dimensions of the grid, which the data should be binned to.
#' @param turn A boolean specifying whether the output should be rotated to keep the previous orientation when the matrix
#' is plotted with the image function.
#' @return A matrix containing the weights of the measure in each bin.
#' @export
bin2d<-function(data.pos,data.weights,gridsize,turn=FALSE){
  M.out<-matrix(0,gridsize[1],gridsize[2])
  data.pos<-t(t(data.pos)*gridsize)
  for (i in 1:gridsize[1]){
    for (j in 1:gridsize[2]){
      lb1<-i-0.5
      lb2<-j-0.5
      ub1<-i+0.5
      ub2<-j+0.5
      bool<-data.pos[,1]>=lb1 & data.pos[,1]<ub1 & data.pos[,2]>=lb2 & data.pos[,2]<ub2
      w<-sum(data.weights[bool])
      M.out[i,j]<-w
    }
  }
  if (turn==TRUE){
    M.out<- apply(t(M.out),2,rev) 
  }
  return(M.out)
}
