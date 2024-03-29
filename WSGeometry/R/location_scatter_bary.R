#' Computes the 2-Wasserstein barycenter of location-scatter families of measures
#' @description This function solves the 2-Wasserstein barycenter problem of N measures from a location-scatter 
#' family, where each data distribution is given by a mean vector and a covariance matrix. In particular, this can be
#' used to compute the barycenter of Gaussian distributions.
#' @param means A list of mean vectors of the elements of the location-scatter family.
#' @param cov A list of semipositive-definite covariance matrices of the elements of the 
#' location-scatter family.
#' @param thresh A real number specifying the threshold for terminating the iterative algorithm.
#' @param maxiter An integer specifying after how many iterations the algorithm should be terminated
#' even if the specified threshold has not been reached.
#' @param showIter A boolean specifying whether the number of performed iterations should be shown at the end.
#' @return A list of two elements. The first element "mean" gives the mean 
#' vector of the barycenter measure. The second elemenent "cov" gives the covariance matrix of the 
#' barycenter measure. 
#' @examples
#' #One dimensional example
#' mean.list<-list(5,15)
#' var.list<-list(1,4)
#' res<-location_scatter_bary(mean.list,var.list)
#' x<-seq(0,22,10^-4)
#' y1<-dnorm(x,mean=5,sd=1)
#' y2<-dnorm(x,mean=15,sd=2)
#' y3<-dnorm(x,mean=res$mean,sd=sqrt(res$cov))
#' plot(x,y1,type="l",main = "Barycenter of two 1-d Gaussian distributions",
#' ylab = "density",col="blue")
#' lines(x,y2,col="green")
#' lines(x,y3,col="red")
#' legend(15,0.4, legend=c("N(5,1)", "N(15,4)","Barycenter"),
#'        col=c("blue","green","red"),lty=1,cex=0.9)
#'        
#'        
#' #two dimensional example
#' # packages graphics and mvtnorm are required to run this example
#'set.seed(2898581)
#'mean.list <- list(c(0,0), c(0,0), c(0,0))
#'COV <- 0.3 + rWishart(3, df = 2, Sigma = diag(2))
#'cov.list <- list(COV[,, 1], COV[,, 2], COV[,, 3])
#'res<-location_scatter_bary(mean.list, cov.list)
#'
#'x <- y <- seq(-3, 3, .1)
#'z <- array(0.0, dim = c(length(x), length(y), 4))
#'for(i in seq_along(x))
#'  for(j in seq_along(y))
#'  {
#'    for(n in 1:3)
#'      z[i, j, n] <- mvtnorm::dmvnorm(c(x[i], y[j]), sigma = COV[, , n])
#'    
#'    z[i, j, 4] <- mvtnorm::dmvnorm(c(x[i], y[j]), sigma = res$cov)
#'  }
#'
#'op <- par(mfrow = c(2, 2),  mai = c(0, 0, 0, 0))
#'for(n in 1:3)
#'{
#'  graphics::persp(x, y, z[, , n], theta = 30, phi = 30, expand = 0.5, col = "lightblue",
#'   zlab = "", ticktype = "detailed", shade = .75, lphi = 45, ltheta = 135)
#'  text(x = 0, y = 0.2, labels = paste("COV[,, ", n, "]", sep = ""))
#'}
#'
#'graphics::persp(x, y, z[, , 4], theta = 30, phi = 30, expand = 0.5, col = "red", zlab = "",
#' ticktype = "detailed", shade = .75, lphi = 45, ltheta = 135)
#'text(x = 0, y = 0.2, labels = "barycenter")
#'par(op)
#' @references PC Álvarez-Esteban, E del Barrio, JA Cuesta-Albertos, and C Matrán (2016).
#'A fixed-point approach to barycenters in Wasserstein space. J. Math. Anal.
#'Appl., 441(2):744–762.\cr
#' Y Zemel and VM Panaretos (2019). Fréchet Means and Procrustes Analysis in Wasserstein Space. Bernoulli 25(2):932-976.
#' @export
location_scatter_bary<-function(means,cov,thresh=10^(-5),maxiter=100,showIter=FALSE){
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("This method requires expm to run. Please install it to use this function,")
  }
  S.list<-cov
  N<-length(S.list)
  d<-dim(S.list[[1]])
  if (is.null(d)){
    d<-rep(0,2)
    d[1]<-1
    d[2]<-length(S.list[[1]])
  }
  out.mean<-Reduce("+",means)/N
  S<-diag(1,d[1],d[2])
  S.old<-diag(0,d[1],d[2])
  err<-Inf
  iter<-0
  while(err>thresh){
    tic<-proc.time()
    S.old<-S
    S.sqrt<-expm::sqrtm(S)
    S.sqrt.m<-expm::sqrtm(solve(S))
    S.sum<-lapply(S.list,sqrtm_transform,S.sqrt)
    S.sum<-Reduce("+",S.sum)/N
    S.sum<-S.sum%*%S.sum
    S<-S.sqrt.m%*%S.sum%*%S.sqrt.m
    
    
    S.sqrt<-expm::sqrtm(S)
    S.sqrt.m<-expm::sqrtm(solve(S))
    S.sum<-lapply(S.list,sqrtm_transform,S.sqrt)
    S.sum<-Reduce("+",S.sum)/N
    
    
    err.mat<-S-S.sum
    err<-sum(abs(err.mat))
    iter<-1
  }
  if (showIter){
    print(iter)
  }
  return(list(mean=out.mean,cov=S))
}

