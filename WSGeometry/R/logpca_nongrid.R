#'Computes Wasserstein principal components
#'@description Computes principal components in the 2-Wasserstein Space for a dataset of weighted point measures in R^2.
#'@details This function computes the principal components of a dataset consisting of weighted point measures in R^2.
#'To do this it first maps the data to the tangent space at the barycenter and then performs standard Euclidean PCA on this space.
#'Afterwards the resulting components are mapped back to the 2-Wasserstein space.
#'@param data.list A list of objects of which the principal components should be computed. Each element should be one of the following:
#' A matrix, representing an image; A file name containing an image; 
#' A \link[transport]{wpp-object}; 
#' A \link[transport]{pp-object};
#' A list containing an entry `positions` with the support of the measure and an entry `weights` containing the weights of the support points;
#' A list containing en entry `positions` specifying the support of a measure with uniform weights.  
#' @param barycenter A barycenter of the dataset. See data.list for possible object types. 
#'@param pca.count An integer specifying the number of principal components to be computed.
#'@param steps_number An integer specifying the number of discretisation steps for the output of the principal components.
#'@return A list with two entries. The first contains a list (of length pca.count), where each entry is a list of \link[transport]{wpp-object}s specifying one point on the corresponding principal
#'component. The second entry is a vector containing the eigenvalues of the computed principal components. The output can be plotted by applying
#'\link[WSGeometry]{plotGeodesic} to each component of the list.
#'@export
#'@examples 
#' set.seed(2020)
#' N<-20
#' supp.size<-10^2
#' L<-sqrt(supp.size)
#' d<-2
#' data.list<-vector("list",N)
#' image.list<-vector("list",N)
#' for (i in 1:N){
#'   t.vec<-seq(0,2*pi,length.out=supp.size)
#'   pos<-cbind(cos(t.vec)*runif(1,0.2,1),sin(t.vec)*runif(1,0.2,1))
#'   theta<-runif(1,0,2*pi)
#'   rotation<-matrix(c(cos(theta),sin(theta),-1*sin(theta),cos(theta)),2,2)
#'   pos<-pos%*%rotation
#'   pos<-pos+1
#'   pos<-pos/2
#'   W<-rep(1/supp.size,supp.size)
#'   data.list[[i]]<-wpp(pos,W)
#' }
#' 
#' res1<-wasserstein_bary(data.list,return_type = "wpp",method="alternating",
#' supp.size = supp.size,warmstartlength = 2,pos_maxIter = 50,
#' weights_maxIter = 0,maxIter = 1,stepsize=1)
#' pcomps<-ws_logpca(data.list,res1,3)
#' ## Set the image and/or gif flags to TRUE to run the example. 
#' ## CRAN policy prevents examples from generating files in the working directory,
#' ## so this had to be disabled.
#' plotGeodesic(pcomps$components[[1]],File="PCA1",images=FALSE,gif=FALSE)
#' plotGeodesic(pcomps$components[[2]],File="PCA2",images=FALSE,gif=FALSE)
#' plotGeodesic(pcomps$components[[3]],File="PCA3",images=FALSE,gif=FALSE)
#'
#' @references E Cazelles, V Seguy, J Bigot, M Cuturi, and N Papadakis (2017); Log-PCA versus Geodesic PCA of histograms 
#' in the Wasserstein space. SIAM Journal on Scientific Computing 40(2):B429–B456. \cr
#' W Wei, D Slepcev, S Basu, JA Ozolek, and GK Rohde (2013). A Linear Optimal Transportation Framework for Quantifying and Visualizing Variations in Sets of Images. International Journal of Computer Vision 101(2):254-269.
ws_logpca<-function(data.list,barycenter,pca.count,steps_number=21){
  if (!requireNamespace("transport", quietly = TRUE)) {
    stop("This method requires transport to run. Please install it to use this function,")
  }
  if (!requireNamespace("RSpectra", quietly = TRUE)) {
    stop("This method requires RSpectra to run. Please install it to use this function,")
  }
  N<-length(data.list)
  data.types<-lapply(data.list,type_check)
  data.list<-mapply(process_data,data.list,data.types,SIMPLIFY = FALSE)
  data.wpp<-mapply(wpp,lapply(data.list,"[[",1),lapply(data.list,"[[",2),SIMPLIFY = FALSE)
  type<-type_check(barycenter)
  bary<-process_data(barycenter,type)
  bary<-wpp(bary$positions,bary$weights)
  P<-lapply(data.wpp,transport::transport,a=bary,method="networkflow",fullreturn=TRUE)
  P<-lapply(P,"[[",2)
  V<-mapply(bary_proj,lapply(data.list,"[[",1),P,MoreArgs = list(bary=bary$mass),SIMPLIFY=FALSE)
  V<-lapply(V,as.vector)
  V.mat<-t(matrix(unlist(V),2*length(bary$mass),N))
  V.bar<-Reduce('+',V)/N
  tmp<-lapply(V,cov_mat,mean=V.bar)
  S<-Reduce('+',tmp)/N
  decomp<-RSpectra::eigs_sym(S,pca.count)
  EigV<-decomp$vectors/sqrt(bary$mass)
  tV<-t(t(V.mat)-V.bar)%*%diag(c(bary$mass,bary$mass))%*%EigV
  PC<-vector("list",pca.count)
  for (c in 1:pca.count){
    TV<-EigV[,c]*sqrt(decomp$values[c])
    PC[[c]]<-matrix(TV,length(bary$mass),2)
  }
  ###turn velocity fields in geodesics
  steps<-seq(-1,1,length.out=steps_number)
  out.list<-vector("list",pca.count)
  for (k in 1:pca.count){
    sub.list<-vector("list",steps_number)
    V<-PC[[k]]
    for (i in 1:steps_number){
      sub.list[[i]]<-wpp(bary$coordinates+(V*steps[i]),bary$mass)
    }
    out.list[[k]]<-sub.list
  }
  return(list(components=out.list,evalues=decomp$values))
}

