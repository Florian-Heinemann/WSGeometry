#' Compute Wasserstein barycenters
#' @description This function computes the Wasserstein barycenter of a list of suitable objects and returns the barycenter in a prespecified form. 
#' @param data.list A list of objects of which the barycenter should be computed. Each element should be one of the following:
#' A matrix, representing an image; A path to a file containing an image; 
#' A \link[transport]{wpp-object}; 
#' A \link[transport]{pp-object};
#' A list containing an entry named `positions` with the support of the measure and an entry named `weights` containing the weights of the support points;
#' A list containing en entry named `positions`` specifying the support of a measure with uniform weights.  
#' @param frechet.weights A real vector summing to 1, specifying the weights in the Frechet functional.  Should be of the same length as data.list. 
#' @param method A string specifiying the method to be used. This also determines which of the other parameters are active/used in this function call. See details
#' for the specific methods currently available.
#' @param return_type A string specifying the format of the output. The currently available options are "default" (which gives list with entries `positions` and `weights`);
#' "wpp"- which gives a \link[transport]{wpp-object}; and "image_mat" for a matrix of the same dimensions as the input matrices (only for the regular method).
#' @param supp.size A positive integer specifying the size of the support used to approximate the barycenter in the "alternating" method.
#' @param output.supp An Mxd matrix specifying the support set on which the optimal weights of the barycenter should be approximated when method = "fixed_support". Each row of the matrix represents one support point in R^d.
#' @param shared A boolean flag specifying whether all measures have the same support set and the weights of the barycenter should be optimised over this set as well.
#' @param sample.size A positive integer specifying the number of samples drawn in the stochastic approximation of the barycenter for method "sampling". 
#' @param maxIter A positive integer specifyng the maximum number of "outer" iterations. The full number of iteration steps performed is 
#' maxIter * (weights_maxIter+pos_maxIter).
#' @param weights_maxIter A positive integer specifying the maximum number of iterations on the weights of the barycenter.
#' @param pos_maxIter A positive integer specifying the maximum number of iterations on the support of the barycenter. 
#' @param stepsize A positive number specifying the stepsize in the position iterations. 
#' @param thresh A positive number specifying the minimal amount of change between iterations, which does not cause the algorithm to terminate.
#' @param regular A positive number specifying the regularisation parameter in the "regular" method. 
#' @param warmstart A boolean specifying whether the algorithm should use a warmstart based on a stochastic subgradient descent. 
#' @param warmstartlength A positive integer specifying the length of the warmstart. The number of steps in the SGD in the warmstart is `length(data.list)*warmstartlength`.
#' @param showIter A boolean specifying whether the number of "outer" iterations performed should be shown at the end.
#' @param threads A positive integer specifying the number of threads used for parallel computing. 
#' @details This is the main function of this package.  It computes/approximates 2-Wasserstein barycenters using the different methods outlined in the following. \cr
#' "lp". Here the barycenter problem can be posed as a linear program. This method builds and solves this linear program. While this gives exact
#' solutions to the problem, this method is highly run-time extensive and should only be used on small datasets. (See Anderes et al. (2016) and Borgwardt & Patterson (2020) for details). \cr
#' "regular". This method solves the entropy-regularised fixed support barycenter problem. Here, a penalisation term is introduced to the problem, which 
#' yields a strictly convex problem that can be solved with Sinkhorn's algorithm (for details see Benamou et al. (2015)). Additionally, it is 
#' assumed that all the measures have the same support set, and instead of an exact (regularised) barycenter, the methods finds the best solution
#' having the same support as the data. This is quite reasonable when the dataset consists of images and the barycenter should be an image as well.  The choice of the
#' regularisation parameter "regular" is a delicate issue. Large values reduce the runtime, but yield "blurry" barycenters.  Small values
#' yield sharper results, but have longer run-time and may cause numerical instabilities. The choice of this parameter depends on the dataset at hand,
#' and will typically require tuning. \cr
#' "fixed_support". This method computes the best approximation of the barycenter, which is supported on a pre-specified support set (as supplied by the parameter "output.support"). Contrary to the "regular"
#' method, here this set does not need to coincide with any of the support sets of the data measures.  See Cuturi & Doucet (2014) for details. \cr
#' "alternating". This method computes the best approximation of the barycenter with a certain support size. It alternates between finding the best positions
#' for given weights and then finding the best weights for these positions. See Cuturi and Doucet (2014) for details. \cr
#' "sampling". This method uses the SUA method of Heinemann et al. (2020) to generate a stochastic approximation of the barycenter.  It 
#' replaces the original measures by empirical measures obtained from samples of size 'sample.size' from each data measure.\cr
#' The unregularised optimal transport problems, which need to be solved for the iterative methods, without regularisation, in each iteration step, are solved using
#' a fast network simplex implementation (which is a modification of the LEMON Library by Nicolas Bonneel).
#' @examples 
#' ##Basic Examples
#' #build list
#' K<-1
#' N<-4*K
#' M<-9
#' d<-2
#' data.list<-vector("list",N)
#' 
#' ###image_mat
#' for (i in 1:K){
#'   U<-runif(M)
#'   U<-U/sum(U)
#'   data.list[[i]]<-matrix(U,sqrt(M))
#' }
#' 
#' #wpp
#' for (i in (K+1):(2*K)){
#'   U<-runif(M)
#'   U<-U/sum(U)
#'   pos<-matrix(runif(d*M),M,d)
#'   data.list[[i]]<-transport::wpp(pos,U)
#' }
#' 
#' #point pattern
#' for (i in (2*K+1):(3*K)){
#'   pos<-matrix(runif(d*M),M,d)
#'   data.list[[i]]<-list(positions=pos)
#' }
#' 
#' #weighted point pattern
#' 
#' for (i in (3*K+1):(4*K)){
#'   U<-runif(M)
#'   U<-U/sum(U)
#'   pos<-matrix(runif(d*M),M,d)
#'   data.list[[i]]<-list(positions=pos,weights=U)
#' }
#' 
#' system.time(res1<-wasserstein_bary(data.list,return_type = "wpp",method="lp"))
#' frechet_func(res1,data.list)
#' \donttest{
#' system.time(res2<-wasserstein_bary(data.list,return_type = "wpp",method="alternating",
#' supp.size = M*N-N+1,warmstartlength = 3,pos_maxIter = 100,weights_maxIter = 100))
#' frechet_func(res2,data.list)
#' }
#' system.time(res3<-wasserstein_bary(data.list,return_type = "wpp",
#' method="fixed_support",warmstartlength = 3,weights_maxIter = 100,output.supp = res1$coordinates))
#' frechet_func(res3,data.list)
#' system.time(res4<-wasserstein_bary(data.list,return_type = "wpp",
#' method="sampling",sample.size=8,warmstartlength = 3,pos_maxIter = 100))
#' frechet_func(res4,data.list)
#' \donttest{
#' ##Visual Examples
#' ###ellipses
#' set.seed(420)
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
#'   data.list[[i]]<-transport::wpp(pos,W)
#'   I<-bin2d(data.list[[i]]$coordinates,data.list[[i]]$mass,c(L*2,L*2))
#'   I<-smear(I,1,1)
#'   I<-I/sum(I)
#'   image.list[[i]]<-I
#' }
#' 
#' system.time(res1<-wasserstein_bary(data.list,return_type = "wpp",method="alternating"
#' ,supp.size = supp.size,warmstartlength = 0,pos_maxIter = 10,weights_maxIter = 10,maxIter = 10))
#' plot(res1)
#' system.time(res2<-wasserstein_bary(data.list,return_type = "wpp",
#' method="fixed_support",warmstartlength = 0,weights_maxIter = 50,maxIter=1,
#' output.supp = grid_positions(2*L,2*L)))
#' plot(res2)
#' system.time(res3<-wasserstein_bary(data.list,return_type = "wpp",method="sampling",
#' sample.size=400,warmstartlength = 0,pos_maxIter = 100,stepsize = 1,maxIter=1))
#' plot(res3)
#' 
#' system.time(res4<-wasserstein_bary(image.list,return_type = "wpp",
#' method="regular",stepsize = 1,weights_maxIter = 50))
#' plot(res4)
#' system.time(res5<-wasserstein_bary(image.list,return_type = "wpp",
#' method="fixed_support",shared=TRUE,warmstartlength = 0,weights_maxIter = 50,maxIter=1,
#' output.supp = grid_positions(2*L,2*L)))
#' plot(res5)
#' }
#'@references E Anderes, S Borgwardt, and J Miller (2016). Discrete Wasserstein barycenters: 
#' Optimal transport for discrete data. Mathematical Methods of Operations Research, 84(2):389-409. \cr
#' S Borgwardt and S Patterson (2020). Improved linear programs for discrete barycenters.
#' Informs Journal on Optimization 2(1):14-33.\cr
#' J-D Benamou, G Carlier, M Cuturi, L Nenna, and G Peyré (2015). 
#' Iterative Bregman projections for regularized transportation problems. SIAM Journal on Scientific Computing 37(2):A1111-A1138. \cr
#' M Cuturi and A Doucet (2014). Fast Computation of Wasserstein Barycenters. Proceedings of the 31st International Conference on Machine Learning, PMLR 32(2):685-693. \cr
#' F Heinemann, A Munk, and Y Zemel (2020). Randomised Wasserstein barycenter computation: Resampling with statistical guarantees. arXiv preprint.\cr
#' N. Bonneel (2018). Fast Network Simplex for Optimal Transport. \cr Github repository, nbonneel/network_simplex. \cr
#' N. Bonneel, M. van de Panne, S. Paris and W. Heidrich (2011). Displacement interpolation using Lagrangian mass transport. ACM Transactions on Graphics (SIGGRAPH ASIA 2011) 30(6).
#'@export
wasserstein_bary<-function(data.list,frechet.weights=NULL,method="alternating",return_type="wpp",supp.size=NULL,
                        output.supp=NULL,shared=FALSE,sample.size=NULL,maxIter=10,weights_maxIter=100,
                        pos_maxIter=100,stepsize=0.1,thresh=10^(-12),regular=10^-3,warmstart=TRUE,warmstartlength=2,showIter=FALSE, threads=1){
  
  

  
  N<-length(data.list)
  if (!is.list(data.list)){
    stop("Input data is not a list")
  }
  types<-lapply(data.list,type_check)
  if (!all(is.element(types,c("image_mat","image_file","wpp","weighted point pattern", "point pattern")))){
    stop("One or multiple entries of the input list are of an unsupported type!")
  }
  if (is.null(frechet.weights)){
    frechet.weights<-rep(1/N,N)
  }
  if (length(frechet.weights)!=N){
    stop("The length of frechet.weights should be equal to the length of data.list")
  }
  if (!identical(sum(frechet.weights),1)){
    warning("The weights of the Frechet functional do not sum to 1 and have been rescaled.")
    frechet.weights<-frechet.weights/sum(frechet.weights)
  }  
  if (!is.null(output.supp)){
    if (output.supp[1]=="grid"){
      G<-sqrt(supp.size)
      output.supp<-grid_positions(G,G)
    }
  }

  
  if ((method=="lp")&(return_type=="image")){
    stop("The use of the LP formulation for image data has been disabled due to run-time constraints.")
  }
  if ((method=="lp")){
    data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
    bary<-barycenter_lp(lapply(data.list,"[[",1),lapply(data.list,"[[",2),frechet.weights)
    if ((return_type=="wpp")){
      return(transport::wpp(bary$positions,bary$weights))
    }
    if (return_type=="default"){
      return(list(positions=bary$positions,weights=bary$weights))
    }
    else{
      warnings("Unknown return_type: Using default output instead.")
      return(list(positions=bary$positions,weights=bary$weights))
    }
  }
  if ((method=="regular")&(shared==TRUE)){
    if (!all(is.element(types,c("image_mat","image_file")))){
      stop("This method is only available for data with shared support or image data.")
    }
    data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
    support<-data.list[[1]]$positions
    C<-gen_cost((support),(support))
    weight.mat<-matrix(unlist(data.list),length(data.list[[1]]),N)
    bary.ret<-bary_sinkhorn_arma(weight.mat,matrix(frechet.weights),weights_maxIter,regular,C,thresh,threads)
    bary.est<-bary.ret$bary
    iter.run<-bary.ret$iterations
    if (showIter){
      print(iter.run)
    }
    if ((return_type=="wpp")){
      return(transport::wpp(support,bary.est))
    }
    if (return_type=="default"){
      return(list(positions=support,weights=bary.est))
    }
    else{
      warnings("Unknown or unsupported return_type: Using default output instead.")
      return(list(positions=support,weights=bary.est))
    }
  }
  if ((method=="regular")){
    if (!all(is.element(types,c("image_mat","image_file")))){
      stop("This method is only available for data with shared support or image data.")
    }
    data.list<-mapply(process_data,data.list,types,MoreArgs = list(return_type="image_mat"),SIMPLIFY = FALSE)
    return_grid<-dim(data.list[[1]])
    support<-grid_positions(return_grid[1],return_grid[2])
    C<-gen_cost((support),(support))
    weight.mat<-matrix(unlist(data.list),length(data.list[[1]]),N)
    bary.ret<-bary_sinkhorn_arma(weight.mat,matrix(frechet.weights),weights_maxIter,regular,C,thresh,threads)
    bary.est<-bary.ret$bary
    iter.run<-bary.ret$iterations
    if (showIter){
      print(iter.run)
    }
    if (return_type=="image_mat"){
      return(matrix(bary.est,return_grid[1],return_grid[2]))
    }
    if ((return_type=="wpp")){
      return(transport::wpp(support,bary.est))
    }
    if (return_type=="default"){
      return(list(positions=support,weights=bary.est))
    }
    else{
      warnings("Unknown or unsupported return_type: Using default output instead.")
      return(list(positions=support,weights=bary.est))
    }
  }
  if (method=="fixed_support"){
    if (shared){
      if (all(is.element(types,c("image_mat","image_file")))){
        data.list<-mapply(process_data,data.list,types,MoreArgs = list(return_type="image_mat"),SIMPLIFY = FALSE)
        return_grid<-dim(data.list[[1]])
        support<-grid_positions(return_grid[1],return_grid[2])
      }
      else{
        data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
        support<-data.list[[1]]$positions
      }
      for (i in 1:N){
        data.list[[i]]<-list(positions=support,weights=c(data.list[[i]]))
      }
      bary<-wsbary_cxx_armaP(lapply(lapply(data.list,"[[",2),matrix),support,
                             lapply(data.list,"[[",1),matrix(frechet.weights),TRUE,maxIter,weights_maxIter,0,stepsize,thresh,warmstart,warmstartlength,threads)
    }
    else{
      support<-output.supp
      types<-lapply(data.list,type_check)
      data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
      bary<-wsbary_cxx_armaP(lapply(lapply(data.list,"[[",2),matrix),support,
                             lapply(data.list,"[[",1),matrix(frechet.weights),FALSE,maxIter,weights_maxIter,0,stepsize,thresh,warmstart,warmstartlength,threads)
    }
    iter.run<-bary[[3]]
    if (showIter){
      print(iter.run)
    }
    if ((return_type=="wpp")){
      return(transport::wpp(support,bary[[2]]))
    }
    if (return_type=="default"){
      return(list(positions=support,weights=bary[[2]]))
    }
    else{
      warnings("Unknown or unsupported return_type: Using default output instead.")
      return(list(positions=support,weights=bary[[2]]))
    }
  }
  
  if (method=="alternating"){
    if (is.null(supp.size)){
      warning("supp.size has not been set and will be set automatically")
      types<-lapply(data.list,type_check)
      data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
      supp.size<-max(unlist(lapply(lapply(data.list,"[[",2),length)))
    }
    types<-lapply(data.list,type_check)
    data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
    start_support<-matrix(runif(supp.size*(dim(data.list[[1]]$positions)[2])),supp.size,(dim(data.list[[1]]$positions)[2]))
    bary<-wsbary_cxx_armaP(lapply(lapply(data.list,"[[",2),matrix),start_support,
                           lapply(data.list,"[[",1),matrix(frechet.weights),FALSE,maxIter,weights_maxIter,pos_maxIter,stepsize,thresh,warmstart,warmstartlength,threads)
    iter.run<-bary[[3]]
    if (showIter){
      print(iter.run)
    }
    if ((return_type=="wpp")){
      return(transport::wpp(bary[[1]],bary[[2]]))
    }
    if (return_type=="default"){
      return(list(positions=bary[[1]],weights=bary[[2]]))
    }
    else{
      warnings("Unknown or unsupported return_type: Using default output instead.")
      return(list(positions=bary[[1]],weights=bary[[2]]))
    }
  }
  
  if (method=="sampling"){
    if (is.null(sample.size)){
      stop("Sample size has not been specified!")
    }
    types<-lapply(data.list,type_check)
    data.list<-mapply(process_data,data.list,types,SIMPLIFY = FALSE)
    data.list<-data_sample(lapply(data.list,"[[",1),lapply(data.list,"[[",2),sample.size)
    supp.size<-sample.size
    start_support<-matrix(runif(supp.size*(dim(data.list$positions[[1]])[2])),supp.size,(dim(data.list$positions[[1]])[2]))
    bary<-wsbary_cxx_armaP(lapply(data.list$weights,matrix),start_support,
                           data.list$positions,matrix(frechet.weights),FALSE,maxIter,0,pos_maxIter,stepsize,thresh,warmstart,warmstartlength,threads)
    iter.run<-bary[[3]]
    if (showIter){
      print(iter.run)
    }
    if ((return_type=="wpp")){
      return(transport::wpp(bary[[1]],bary[[2]]))
    }
    if (return_type=="default"){
      return(list(positions=bary[[1]],weights=bary[[2]]))
    }
    else{
      warnings("Unknown or unsupported return_type: Using default output instead.")
      return(list(positions=bary[[1]],weights=bary[[2]]))
    }
  }
  stop("Unknown method input.")
}
