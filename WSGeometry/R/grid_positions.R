#' Generate a 2d grid in [0,1]^2 of given size.
#' @description Generates a matrix containing the positions of the points of an equidistant 2d grid in [0,1]^2.
#' @param n Integer giving one dimension of the grid.
#' @param m integer giving the other dimension of the grid.
#' @return A (nm)x2 grid containing the positions of the desired grid. 
#' @examples 
#' grid<-grid_positions(32,32)
#' grid.wpp<-wpp(grid,rep(1,1024))
#' wppplot(grid.wpp)
#' @export
grid_positions<-function(n,m){
  G1<-expand.grid(n:1,1:m)
  G2<-cbind(G1[,2],G1[,1])
  G2<-G2-0.5
  G2[,1]<-G2[,1]/m
  G2[,2]<-G2[,2]/n
  return(G2)
}