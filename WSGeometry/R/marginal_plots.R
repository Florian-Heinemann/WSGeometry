#' Plot transport between 1D marginals.
#' @description Generates a plot displaying two one-dimensional marginal distribution and a specified transportplan between them.
#' @param A A list with the entries "positions" and "weights" specifying a one-dimensional measure.
#' @param B A list with the entries "positions" and "weights" specifying a one-dimensional measure.
#' @param plan A matrix representing a transport plan between A and B. 
#' @param A.type A control parameter specifying the type of marginal plot. The option "line" yields a lineplot and the option
#' "bars" yields a bar plot.
#' @param B.type A control parameter specifying the type of marginal plot. The option "line" yields a lineplot and the option
#' "bars" yields a bar plot.
#' @param A.col A vector one two colours specifying the colours uses for the marginal plot of A. For "line" the first entry 
#' specifies the colour below the line and the second one the colour of the line. For "bars" the first entry specifies the 
#' colour of the bar and the second one the colour of the dot.
#' @param B.col A vector one two colours specifying the colours uses for the marginal plot of A. For "line" the first entry 
#' specifies the colour below the line and the second one the colour of the line. For "bars" the first entry specifies the 
#' colour of the bar and the second one the colour of the dot.
#' @param A.size For "line" a real number specifying the width of the line. For "bars" a vector of two real numbers, where 
#' the first specifies the size of the dot and the second one the width of the bar.
#' @param B.size For "line" a real number specifying the width of the line. For "bars" a vector of two real numbers, where 
#' the first specifies the size of the dot and the second one the width of the bar.
#' @param plan.type A control parameter specifying the type of plot for the transport plan. "line" provides a lineplot.
#' "pointsCol" displays the plan in terms of dots where the colour decodes the mass at a given location. "pointsSize" displays
#' the plan in terms of dots where the size of a point decodes its mass.
#' @param plan.size For "line" a real number specifying the width of the line. For "pointsCol" a real number specifying the 
#' size of the dots. For "pointsSize" a vector of two real numbers specifying the minimal and maximal size of the dots.
#' @param plan.col For "line" a colour specifying the colour of the line. For "pointsSize" a colour specifying the colour of 
#' the points. For "pointsCol" a vector of two colours, where the first corresponds to the colour of the minimal mass and the 
#' second to the colour of the maximum of the masses in the plan.
#' @param grid A vector of two booleans specifiyng whether a grid should be drawn. The first entry specifies whether the lines
#' starting from the left marginal are drawn and the second whether the ones from the one at the top are drawn.
#' @param grid.type A control parameter specifying the type of line used for the grid. See \link[ggplot2]{linetype} for details.       
#' @param grid.col A colour specifying the colour of the lines of the grid.
#' @param grid.width A real number specifying the width of the lines of the grid.
#' @param axis.width A real number specifying the width of the axis.
#' @param axis.length A number between 0 and 1 specifying the length of the axis drawn after their intersection point.                    
#' @return This function does not provide any return value, but is only called for the plot it generates.
#' @examples 
#' #Example Discrete
#' x.vec<-matrix(seq(-8,8,1))
#' y.vec<-matrix(seq(-8,8,2))
#' W1<-dnorm(x.vec,3,1)+dnorm(x.vec,0,4)+dnorm(x.vec,0,1)
#' W1<-W1/sum(W1)
#' W2<-dnorm(y.vec,-3,0.4)+dnorm(y.vec,-1,2)+dnorm(y.vec,2,1.5)
#' W2<-W2/sum(W2)
#' G1<-list(positions=x.vec,weights=W1)
#' G2<-list(positions=y.vec,weights=W2)
#' wpp1<-WSGeometry:::wpp(x.vec,W1)
#' wpp2<-WSGeometry:::wpp(y.vec,W2)
#' cmat<-t((WSGeometry:::gen_cost(x.vec,y.vec)))
#' plan<-transport::transport(W1,W2,cmat,method="networkflow",fullreturn=TRUE)$primal
#' marginal_plot(G1,G2,plan,A.type ="bars",B.type="bars",A.size=c(2,4),B.size=c(2,4),
#' plan.size = c(0,6),plan.col=c("purple"),axis.length = 0.5,plan.type = "pointsSize",
#' grid=c(TRUE,TRUE),grid.type = "dashed",grid.col="lightgrey")
#' 
#' 
#' 
#' 
#' 
#' #Example Semi-Discrete
#' x.vec<-matrix(seq(-8,8,0.01))
#' y.vec<-matrix(seq(-8,8,2))
#' W1<-dnorm(x.vec,3,1)+dnorm(x.vec,0,4)+dnorm(x.vec,0,1)
#' W1<-W1/sum(W1)
#' W2<-dnorm(y.vec,-3,0.4)+dnorm(y.vec,-1,2)+dnorm(y.vec,2,1.5)
#' W2<-W2/sum(W2)
#' G1<-list(positions=x.vec,weights=W1)
#' G2<-list(positions=y.vec,weights=W2)
#' wpp1<-WSGeometry:::wpp(x.vec,W1)
#' wpp2<-WSGeometry:::wpp(y.vec,W2)
#' cmat<-t((WSGeometry:::gen_cost(x.vec,y.vec)))
#' plan<-transport::transport(W1,W2,cmat,method="networkflow",fullreturn=TRUE)$primal
#' marginal_plot(G1,G2,plan,A.type ="line",B.type="bars",A.size=c(2,4),B.size=c(2,4),
#' plan.size = c(1),plan.col=c("yellow","red"),axis.length = 0.5,
#' plan.type = "pointsCol",grid=c(FALSE,FALSE))
#' 
#' 
#' 
#' 
#' 
#' #Example Continuous
#' x.vec<-matrix(seq(-8,8,0.01))
#' y.vec<-matrix(seq(-8,8,0.01))
#' W1<-dnorm(x.vec,3,1)+dnorm(x.vec,0,4)+dnorm(x.vec,0,1)
#' W1<-W1/sum(W1)
#' W2<-dnorm(y.vec,-3,0.4)+dnorm(y.vec,-1,2)+dnorm(y.vec,2,1.5)
#' W2<-W2/sum(W2)
#' G1<-list(positions=x.vec,weights=W1)
#' G2<-list(positions=y.vec,weights=W2)
#' wpp1<-WSGeometry:::wpp(x.vec,W1)
#' wpp2<-WSGeometry:::wpp(y.vec,W2)
#' cmat<-t((WSGeometry:::gen_cost(x.vec,y.vec)))
#' plan<-transport::transport(W1,W2,cmat,method="networkflow",fullreturn=TRUE)$primal
#' marginal_plot(G1,G2,plan,A.type ="line",B.type="line",A.size=c(2,4),
#' B.size=c(2,4),plan.size = c(1),plan.col=c("black","black"),
#' axis.length = 0.5,plan.type = "line",grid=c(FALSE,FALSE))
#' @export
marginal_plot<-function(A,B,plan,A.type="line",B.type="line",A.col=c("red","darkred"),B.col=c("blue","darkblue"),A.size=c(2,4),B.size=c(2,4),plan.type="pointsCol",plan.size=c(1),plan.col=c("white","purple"),grid=c(FALSE,FALSE),grid.type="dashed",grid.col="lightgrey",grid.width=0.5,axis.width=5,axis.length=0.6){
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("This method requires ggplot2 to run. Please install it to use this function,")
  }
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("This method requires gridExtra to run. Please install it to use this function,")
  }
  A$positions<-as.matrix(A$positions)
  B$positions<-as.matrix(B$positions)
  if (A.type=="line"){
    P2<-mplot(A,linecol = A.col[2],fillcol=A.col[1],width=A.size[1]) 
  }
  if (A.type=="bars"){
    P2<-mplot_bar(A,"top",A.col[1],A.col[2],A.size[1],A.size[2])
  }
  if (B.type=="line"){
    P1<-mplot(B,"left",linecol = B.col[2],fillcol=B.col[1],width=B.size[1])
  }
  if (B.type=="bars"){
    P1<-mplot_bar(B,"left",B.col[1],B.col[2],B.size[1],B.size[2])
  }

  P3<-Tplot(A,B,plan,size=plan.size,plan.col,plan.type,grid[1],grid[2],grid.type,grid.col,grid.width)
  layout<-matrix(c(NA,1,1,1,2,3,3,3,2,3,3,3,2,3,3,3),4,4)
  P<-gridExtra::grid.arrange(P1,P2,P3,layout_matrix=layout)
  axis.shift<-0.25*(1-axis.length)
  P<-P +  ggplot2::annotation_custom(  grid::grid.polygon(c(axis.shift, 1/4, 1, 1/4, 1/4, 1/4), 
          c(3/4, 3/4, 3/4,0, 3/4, 1-axis.shift), id = c(1,1,1,2,2,2),   gp = grid::gpar(lwd = axis.width)))
}

