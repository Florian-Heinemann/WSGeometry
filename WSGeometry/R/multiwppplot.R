#' Plot wpp-objects and transport plans between them.
#' @description This function plots one or two wpp-objects and transport plans between them.
#' @param A.list A list of \link[transport]{wpp-object}s.
#' @param plan.list A list of matrices representing transport plans between the elements of A.list.
#' @param plan.which A list of the same length as plan.list. Each element of the list is a vector of length two
#' specifying the indices of the measures in A.list between the plan is transporting.
#' @param pch A graphical parameter specifying the shape of the locations in the plots. 
#' @param bg A graphical parameter determining the colour of the locations for certain choices of pch.
#' @param col A graphical parameter specigying the colour of the locations or the colour of the border of
#' the locations depending on pch.
#' @param linecol A graphical parameter specifying the colour of the transport plan.
#' @param  cex A graphical parameter specifying the size of the plotted points.
#' @param bwd A graphical parameter specifying the size of the border of the plotted points when applicable.
#' @param lwd A graphical parameter specifying the width of the OT plans. 
#' @param xlim A vector with two entries speciying the range of the plot on the x-axis.
#' @param ylim A vector with two entries speciying the range of the plot on the y-axis.
#' @param clean A graphical parameter specifying whether the axis are plotted. 
#' @param tstart A graphical parameter specifying the visuals of the transport plans. The available options are "center",
#' "border" and "arrow". For "center" the lines start and end in the center of the points. For "border" they start and end 
#' at the border. For "arrow" the transport plan is represented by arrows and the start and end point are controlled by the 
#' parameter arrowScale.
#' @param arrowScale A real number specifying the starting and end positions of the transport plan when using the "arrow" 
#' parameter for tstart.
#' @return This function has no return value and is only called for generating plots.
#' @examples 
#' pos1<-rbind(c(1,1),c(5,9),c(9,1)) #blue
#' pos2<-rbind(c(2,8),c(5,4),c(8,9)) #green
#' pos3<-rbind(c(3,4),c(6,7),c(7.5,3)) #red
#' W1<-rep(1,3)
#' W2<-rep(1,3)
#' W3<-rep(1,3)
#' wpp1<-wpp(pos1,W1)
#' wpp2<-wpp(pos2,W2)
#' wpp3<-wpp(pos3,W3)
#' centroid<-matrix(0,27,2)
#' count<-1
#' for (i in 1:3){
#'   for (j in 1:3){
#'     for (k in 1:3){
#'       centroid[count,]<-(pos1[i,]+pos2[j,]+pos3[k,])/3
#'       count<-count+1
#'     }
#'   }
#' }
#' W.cent<-rep(1,27)
#' wpp.cent<-wpp(centroid,W.cent)
#' wpp.full<-wpp(rbind(pos1,pos2,pos3),c(W1,W2,W3))
#' multiwppplot(list(wpp1,wpp2,wpp3,wpp.cent),cex=3,pch=c(21,21,21,24),
#' bg=c("blue","darkgreen","red","purple"),col=rep("black",4),xlim=c(0,10),
#' ylim=c(0,10),clean = TRUE,bwd=2)
#' plan1<-data.frame(from=c(1,2,3),to=c(1,14,27),mass=c(1,1,1))
#' plan2<-data.frame(from=c(1,2,3),to=c(1,14,27),mass=c(1,1,1))
#' plan3<-data.frame(from=c(1,2,3),to=c(1,14,27),mass=c(1,1,1))
#' plan.which<-list(c(1,4),c(2,4),c(3,4))
#' plan.list<-list(plan1,plan2,plan3)
#' multiwppplot(list(wpp1,wpp2,wpp3,wpp.cent),plan.list=plan.list,plan.which=plan.which
#' ,cex=3,pch=c(21,21,21,24),bg=c("blue","darkgreen","red","purple"),col=rep("black",4)
#' ,xlim=c(0,10),ylim=c(0,10),clean = TRUE,bwd=2,lwd=5,linecol="orange",tstart="border")
#' 
#' 
#' @export
multiwppplot<-function(A.list,plan.list=NULL,plan.which=NULL,pch=c(16,16),bg=c("blue","red"),col=c("blue","red"),linecol="black",
                  cex=1,bwd=2,lwd=2,xlim=c(0,1),ylim=c(0,1),clean=FALSE,tstart="center",arrowScale=0.03){
  oldpars <- par(xpd = TRUE, xaxs = "i", yaxs = "i")
  if (!is.null(plan.list)){
    plan.list<-lapply(plan.list,plan_convert,target="frame")
  }
  if (clean){
    plot(A.list[[1]]$coordinates, type = "n", xlab = "", 
         ylab = "",axes=FALSE,frame.plot=TRUE,xlim=xlim,ylim=ylim)
  }
  else{
    plot(A.list[[1]]$coordinates, type = "n", xlab = "", 
         ylab = "",xlim=xlim,ylim=ylim)
  }
  if (tstart=="border"){
    if (!is.null(plan.list)){
      P<-length(plan.list)
      count<-1
      for (plan in plan.list){
        inds<-plan.which[[count]]
        A<-A.list[[inds[1]]]
        B<-A.list[[inds[2]]]
        segments(A$coordinates[plan$from, 1], A$coordinates[plan$from,2], 
                 B$coordinates[plan$to, 1], B$coordinates[plan$to,2], col = linecol, lwd =lwd,lty=1)
        count<-count+1
      }
    }
  }

  N<-length(A.list)
  for(i in 1:N){
    A<-A.list[[i]]
    points(A$coordinates, col = col[i] , 
           pch = pch[i], cex = cex*sqrt(A$mass), lwd = bwd,bg=bg[i],xlim=xlim,ylim=ylim)
  }
  
  
  if (tstart=="center"){
    if (!is.null(plan.list)){
      P<-length(plan.list)
      count<-1
      for (plan in plan.list){
        inds<-plan.which[[count]]
        A<-A.list[[inds[1]]]
        B<-A.list[[inds[2]]]
        segments(A$coordinates[plan$from, 1], A$coordinates[plan$from,2], 
                 B$coordinates[plan$to, 1], B$coordinates[plan$to,2], col = linecol, lwd =lwd,lty=1)
        count<-count+1
      }
    }
  }
  if (tstart=="arrow"){
    if (!is.null(plan.list)){
      P<-length(plan.list)
      count<-1
      for (plan in plan.list){
        inds<-plan.which[[count]]
        A<-A.list[[inds[1]]]
        B<-A.list[[inds[2]]]
        V<-B$coordinates[plan$to,]-A$coordinates[plan$from,]
        V<-t(apply(V,1,normalise))
        Vfac<-arrowScale
        arrows(A$coordinates[plan$from, 1]+Vfac*V[,1], A$coordinates[plan$from,2]+Vfac*V[,2], 
               B$coordinates[plan$to, 1]-Vfac*V[,1], B$coordinates[plan$to,2]-Vfac*V[,2], col = linecol, lwd =lwd)
        count<-count+1
      }
    }
  }
}


