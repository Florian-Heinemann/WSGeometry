#' Plot wpp-objects and transport plans between them.
#' @description This function plots one or two wpp-objects and transport plans between them.
#' @param A A \link[transport]{wpp-object}.
#' @param B A \link[transport]{wpp-object}.
#' @param plan A matrix representing the OT plan between A and B. 
#' @param p An integer specifying the power of the Euclidean distance used for computing the OT plan 
#' if none is provided. 
#' @param C The penalty C of the Kantorovich-Rubinstein distance used for computing the UOT plan if no
#' plan is provided. 
#' @param pch A graphical parameter specifying the shape of the locations in the plots. The first entry determined this
#' for A and the second for B.
#' @param bg A graphical parameter determining the colour of the locations for certain choices of pch. The first entry determined this
#' for A and the second for B.
#' @param col A graphical parameter specifying the colour of the locations or the colour of the border of 
#' the locations depending on pch. The first entry determined this
#' for A and the second for B.
#' @param linecol A graphical parameter specifying the colour of the transport plan.
#' @param  cex A graphical parameter specifying the size of the plotted points.
#' @param bwd A graphical parameter specifying the size of the border of the plotted points when applicable.
#' @param lwd A graphical parameter specifying the width of the OT plans. 
#' @param xlim A vector with two entries specifying the range of the plot on the x-axis.
#' @param ylim A vector with two entries specifying the range of the plot on the y-axis.
#' @param clean A graphical parameter specifying whether the axis are plotted. 
#' @param tstart A graphical parameter specifying the visuals of the transport plans. The available options are "center",
#' "border" and "arrow". For "center" the lines start and end in the center of the points. For "border" they start and end 
#' at the border. For "arrow" the transport plan is represented by arrows and the start and end point are controlled by the 
#' parameter arrowScale.
#' @param arrowScale A real number specifying the starting and end positions of the transport plan when using the "arrow" 
#' parameter for tstart.
#' @return This function has no return value and is only called for generating plots.
#' @examples 
#' M<-8
#' set.seed(76597)
#' pos1<-cbind(runif(M,0.1,0.8),runif(M,0.1,0.9))
#' pos2<-cbind(runif(M,0.1,0.9),runif(M,0.2,0.9))
#' W1<-rep(1,M)
#' W2<-rep(1,M)
#' wpp1<-wpp(pos1,W1)
#' wpp2<-wpp(pos2,W2)
#' wppplot(wpp1,wpp2,plan=transport::transport(wpp1,wpp2,p=2),cex=3,bwd=2,lwd=2,clean = TRUE
#' ,tstart = "border",pch=c(21,21),col=c("black","black"))
#' wppplot(wpp1,wpp2,plan=transport::transport(wpp1,wpp2,p=2),cex=3,bwd=2,lwd=2,clean = TRUE
#' ,tstart = "border",pch=c(22,25),col=c("black","black"),bg=c("green","orange"))
#' wppplot(wpp1,wpp2,plan=transport::transport(wpp1,wpp2,p=2),cex=3,bwd=2,lwd=2,clean = TRUE
#' ,tstart = "center",pch=c(22,25),col=c("blue","red"),bg=c("green","orange"))
#' wppplot(wpp1,wpp2,plan=transport::transport(wpp1,wpp2,p=2),cex=3,bwd=2,lwd=2,clean = TRUE
#' ,tstart = "arrow",pch=c(20,20),col=c("purple","brown"),bg=c("brown","purple"),arrowScale = 0.02)
#' @export
wppplot<-function(A,B=NULL,plan=NULL,p=NULL,C=NULL,pch=c(16,16),bg=c("blue","red"),col=c("blue","red"),linecol="black",
                  cex=1,bwd=2,lwd=2,xlim=c(0,1),ylim=c(0,1),clean=FALSE,tstart="center",arrowScale=0.03){
  oldpars <- par(xpd = TRUE, xaxs = "i", yaxs = "i")
  if (!is.null(plan)){
    plan<-plan_convert(plan,"frame") 
  }
  if (clean){
    plot(A$coordinates, type = "n", xlab = "", 
         ylab = "",axes=FALSE,frame.plot=TRUE,xlim=xlim,ylim=ylim)
  }
  else{
    plot(A$coordinates, type = "n", xlab = "", 
         ylab = "",xlim=xlim,ylim=ylim)
  }
  
  if (is.null(plan)){
    if (!is.null(p)){
      if (!requireNamespace("transport", quietly = TRUE)) {
        stop("If no transport plan is provided, but the parameter p is set, then transport is required to 
             compute an optimal transport plan. Please install it to use this option or provide a transport
             plan of your own.")
      }
      if (is.null(C)){
        if ((abs(A$totmass-B$totmass))>10^(-14)){
          warning("Measures did not have equal total intensity, but balanced OT was chosen. Normalising both measures
                to be probability measures.")
          A.tmp<-wpp(A$coordinates,A$mass/A$totmass)
          B.tmp<-wpp(B$coordinates,B$mass/B$totmass)
        }
        else{
          A.tmp<-A
          B.tmp<-B
        }
        plan<-transport::transport(A.tmp,B.tmp,p=p)
      }
      else{
        uotplan<-kr_dist(A,B,p=p,C)
        plan<-plan_convert(uotplan$plan,"frame")
      }
    } 
  }
  if (!is.null(plan)){
    if (tstart=="border"){
      segments(A$coordinates[plan$from, 1], A$coordinates[plan$from,2], 
               B$coordinates[plan$to, 1], B$coordinates[plan$to,2], col = linecol, lwd =lwd*plan$mass,lty=1)
    }
  }
  points(A$coordinates, col = col[1], 
         pch = pch[1], cex = cex*sqrt(A$mass), lwd = bwd,bg=bg[1],xlim=xlim,ylim=ylim)
  if (!is.null(B)){
    points(B$coordinates, col = col[2], 
           pch = pch[2], cex = cex*sqrt(B$mass), lwd = bwd,bg=bg[2],xlim=xlim,ylim=ylim)
  }
  if (!is.null(plan)){
    if (tstart=="center"){
      segments(A$coordinates[plan$from, 1], A$coordinates[plan$from,2], 
               B$coordinates[plan$to, 1], B$coordinates[plan$to,2], col = linecol, lwd =lwd,lty=1)
    }
    if (tstart=="arrow"){
      V<-B$coordinates[plan$to,]-A$coordinates[plan$from,]
      V<-t(apply(V,1,normalise))
      Vfac<-arrowScale
      arrows(A$coordinates[plan$from, 1]+Vfac*V[,1], A$coordinates[plan$from,2]+Vfac*V[,2], 
             B$coordinates[plan$to, 1]-Vfac*V[,1], B$coordinates[plan$to,2]-Vfac*V[,2], col = linecol, lwd =lwd,lty=1)
    }
  }
}