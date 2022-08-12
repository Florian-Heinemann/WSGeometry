#This file contains a collection of small helpful functions which are called by a wide range of other functions
#within this package. They are not exported to the package's namespace for this reason. 

#Compute the value of the euclidean barycenter functional
euclid_bary_func_w<-function(x,weights){
  K<-dim(x)[1]
  d<-dim(x)[2]
  x<-matrix(unlist(x),K,d)
  x.bar<-matrix(colSums(diag(weights)%*%x),1,d)
  out.sum<-0
  for (k in 1:K){
    out.sum<-out.sum+weights[k]*gen_cost(matrix(x[k,],1,d),x.bar)
  }
  return(out.sum)
}

#Generate a specific column of the constraint matrix of a multi-marginal optimal transport problem.
gen_constraints_multi_col<-function(sizes,index){
  N<-length(sizes)
  index<-index-1
  A<-rep(0,sum(sizes))
  j<-floor(index/prod(sizes[2:N]))
  A[(j+1)]<-1
  if (N>2){
    for (i in 2:(N-1)){
      j<-sum(sizes[1:(i-1)])+floor((index%%prod(sizes[i:N]))/(prod(sizes[(i+1):N])))
      A[(j+1)]<-1
    }
  }
  j<-sum(sizes[1:(N-1)])+(index%%(sizes[N]))
  A[(j+1)]<-1
  return(A)
}

#Generate the constraint matrix of the multi-marginal optimal transport problem.
gen_constraints_multi<-function(sizes){
  A<-matrix(0,sum(sizes),prod(sizes))
  for (i in 1:prod(sizes)){
    A[,i]<-gen_constraints_multi_col(sizes,i)
  }
  return(A)
}

# #Map weighted point pattern to a grid by deploying a 2d Gaussian kernel density estimator.
# smooth2d<-function(data.pos,data.weights,gridsize,H=matrix(c(1,0,0,1)/2000,2,2),turn=FALSE,flip=FALSE){
#   G<-grid_positions(gridsize[1],gridsize[2])
#   data.bin<-ks::kde(data.pos,w=data.weights,eval.points=G,H = H,xmin=c(0,0),xmax=c(1,1))
#   data.mat<-matrix(data.bin$estimate, gridsize[1],gridsize[2])
#   data.mat<-data.mat/sum(data.mat)
#   if (turn==TRUE){
#     data.mat<- apply(t(data.mat),2,rev) 
#   }
#   if (flip==TRUE){
#     data.mat<-data.mat[,gridsize[1]:1]
#   }
#   return(data.mat)
# }

#Generate a empirical data set from an input data set by drawing S independent samples from each measure, respectively.
data_sample<-function(data.pos,data.weights,S){
  N<-length(data.pos)
  sample.pos<-vector("list",N)
  sample.weights<-vector("list",N)
  for (i in 1:N){
    M<-length(data.weights[[i]])
    samp.weights<-rep(0,M)
    samps<-sample(1:M,S,replace=TRUE,prob=data.weights[[i]])
    for (s in 1:S){
      samp.weights[samps[s]]<-samp.weights[samps[s]]+(1/S)
    }
    sample.pos[[i]]<-data.pos[[i]][samp.weights>0,]
    sample.weights[[i]]<-samp.weights[samp.weights>0]
  }
  return(list(positions=sample.pos,weights=sample.weights))
}

# ####Bin data onto a grid.
# 
# bin2d<-function(data.pos,data.weights,gridsize,turn=FALSE){
#   M.out<-matrix(0,gridsize[1],gridsize[2])
#   data.pos<-t(t(data.pos)*gridsize)
#   for (i in 1:gridsize[1]){
#     for (j in 1:gridsize[2]){
#       lb1<-i-0.5
#       lb2<-j-0.5
#       ub1<-i+0.5
#       ub2<-j+0.5
#       bool<-data.pos[,1]>=lb1 & data.pos[,1]<ub1 & data.pos[,2]>=lb2 & data.pos[,2]<ub2
#       w<-sum(data.weights[bool])
#       M.out[i,j]<-w
#     }
#   }
#   if (turn==TRUE){
#     M.out<- apply(t(M.out),2,rev) 
#   }
#   return(M.out)
# }

###########Generate Centroid Superset
generate_superset<-function(positions){
  N<-length(positions)
  d<-dim(positions[[1]])[2]
  sizes<-unlist(lapply(positions,length))/d
  support.size<-prod(sizes)
  vec<-list()
  for (k in 1:(N)){
    sup<-positions[[k]]
    V<-matrix(0,0,d)
    count<-0
    current.pos<-1
    sub.count<-0
    if (k<N){
      sub.max<-prod(sizes[(k+1):N])
    }
    else{
      sub.max<-1
    }
    while (count <support.size){
      tmp<-sup[current.pos,]
      V<-rbind(V,tmp)
      count<-count+1
      sub.count<-sub.count+1
      if (sub.count>=sub.max){
        sub.count<-0
        current.pos<-current.pos+1
        if (current.pos>sizes[k]){
          current.pos<-1
        }
      }
    }
    vec[[k]]<-V
  }
  out<-Reduce('+',vec)/N
  return(out)
}





#######Check type of data object
type_check<-function(object){
  if (!(("weights" %in% names(object))&("positions" %in% names(object)))){
    if (!("positions" %in% names(object))){
      if (!((class(object)[1])=="wpp")){
        if (!is.matrix(object)){
          if (!((is.character(object))&(length(object)=1))){
            return("Input type not supported!")
          }
          else{
            return("image_file")
          }
        }
        else{
          return("image_mat")
        }
      }
      else{
        return("wpp")
      }
    }
    else{
      return("point pattern")
    }
  }
  else{
    return("weighted point pattern")
  }
}
####Convert types of data 
process_data<-function(object,type,return_type="weighted point pattern"){
  if (type=="image_file"){
    if (!requireNamespace("imager", quietly = TRUE)) {
      stop("The package imager is required to use images as input measures. Please install it to use this functionality.")
    }
    I<-imager::load.image(object)
    I<-imager::grayscale(I)
    d<-(attributes(I)$dim)[1:2]
    IM<-I[1:d[1],1:d[2],1,1]
    IM<-IM/sum(IM)
    object<-IM
    type<-"image_mat"
  }
  if ((!(return_type==type))&(type=="image_mat")){
    d<-dim(object)
    G<-grid_positions(d[1],d[2])
    W<-as.vector(object)
    G<-G[W>0,]
    W<-W[W>0]
    object<-list(positions=G,weights=W)
    type<-"weighted point pattern"
  }
  if ((!(return_type==type))&(type=="wpp")){
    object<-list(positions=object$coordinates,weights=as.vector(object$mass))
    type<-"weighted point pattern"
  }
  if ((!(return_type==type))&(type=="point pattern")){
    M<-dim(object[[1]])[1]
    object<-list(positions=object[[1]],weights=rep(1/M,M))
    type<-"weighted point pattern"
  }
  if (type=="weighted point pattern"){
    object<-list(positions=object$positions,weights=object$weights)
    return(object)
  }
  if (type==return_type){
    return(object)
  }
  return(list(positions=NaN,weights=NaN))
}




########Helper Functions for the Multi-Marginal Wrapper Function
select_costInds<-function(A,col,sizes){
  a<-A[,col]
  b<-NULL
  for (k in 1:length(sizes)){
    b<-c(b,1:sizes[k])
  }
  inds<-a*b
  inds<-inds[inds>0]
  return(inds)
}


build_MM_costvec<-function(costA,A){
  D<-dim(costA)
  costvec<-rep(0,prod(D))
  for (r in 1:prod(D)){
    costInds<-select_costInds(A,r,D)
    costvec[r]<-costA[t(costInds)]
  }
  return(costvec)
}

build_MMCoupling<-function(optSol,A,D){
  indBool<-optSol$solution>0
  x.vec<-optSol$solution[indBool]
  MMCoupling<-array(0,D)
  inds<-1:prod(D)
  inds<-inds[indBool]
  count<-1
  for (r in inds){
    costInds<-select_costInds(A,r,D)
    MMCoupling[t(costInds)]<-x.vec[count]
    count<-count+1
  }
  return(MMCoupling)
}



#Generate numbering ordered accordingly for gif-generation.
numbering_gen<-function(N){
  pot<-0
  check<-FALSE
  while(!check){
    pot<-pot+1
    pow<-10^pot
    tmp<-N/pow
    if (tmp<1){
      check<-TRUE
    }
  }
  pot<-pot-1
  strings<-NULL
  zero.count<-0
  for (k in N:1){
    zeros<-gen_zero_string(zero.count)
    if (k==10^pot){
      pot<-pot-1
      zero.count<-zero.count+1
    }
    tmp<-paste(zeros,k,sep="")
    strings<-c(strings,tmp)
  }
  
  return(strings[N:1])
}
#generate string of zeros of given length
gen_zero_string<-function(n){
  if (n==0){
    return("")
  }
  out.string<-""
  for (k in 1:n){
    out.string<-paste(out.string,"0",sep="")
  }
  return(out.string)
}
#Rotate a 3D point pattern along a certain axis for a given angle.
rotate3D<-function(pos,axis,angle){
  R<-matrix(0,3,3)
  R[1]<-cos(angle)+axis[1]^2*(1-cos(angle))
  R[2]<-axis[1]*axis[2]*(1-cos(angle))+axis[3]*sin(angle)
  R[3]<-axis[3]*axis[1]*(1-cos(angle))-axis[2]*sin(angle)
  R[4]<-axis[1]*axis[2]*(1-cos(angle))-axis[3]*sin(angle)
  R[5]<-cos(angle)+axis[2]^2*(1-cos(angle))
  R[6]<-axis[2]*axis[3]*(1-cos(angle))+axis[1]*sin(angle)
  R[7]<-axis[1]*axis[3]*(1-cos(angle))+axis[2]*sin(angle)
  R[8]<-axis[2]*axis[2]*(1-cos(angle))-axis[1]*sin(angle)
  R[9]<-cos(angle)+axis[3]^2*(1-cos(angle))
  return(t(diag(c(2,3,1))%*%(R%*%t(pos))))
}
#Generate a Torus in 3D.
gen_torus<-function(M,R,r){
  theta<-seq(0,2*pi,length.out=M)
  phi<-seq(0,2*pi,length.out=M)
  G<-expand.grid(theta,phi)
  x<-(R+r*cos(G[,1]))*cos(G[,2])
  y<-(R+r*cos(G[,1]))*sin(G[,2])
  z<-r*sin(G[,1])
  return(cbind(x,y,z))
}
#Generate a pentagonal prism in 3D.
draw3Dpentagon<-function(M,G){
  c1<-cos(2*pi/5)
  c2<-cos(pi/5)
  s1<-sin(2*pi/5)
  s2<-sin(pi/5)
  P1<-c(0,1)
  P2<-c(-s1,c1)
  P3<-c(-s2,-c2)
  P4<-c(s2,-c2)
  P5<-c(s1,c1)
  V1<-P2-P1
  V2<-P3-P2
  V3<-P4-P3
  V4<-P5-P4
  V5<-P1-P5
  t.vec<-seq(0,1,length.out = M)
  penta<-NULL
  penta<-rbind(penta,cbind(P1[1]+V1[1]*t.vec,P1[2]+V1[2]*t.vec))
  penta<-rbind(penta,cbind(P2[1]+V2[1]*t.vec,P2[2]+V2[2]*t.vec))
  penta<-rbind(penta,cbind(P3[1]+V3[1]*t.vec,P3[2]+V3[2]*t.vec))
  penta<-rbind(penta,cbind(P4[1]+V4[1]*t.vec,P4[2]+V4[2]*t.vec))
  penta<-rbind(penta,cbind(P5[1]+V5[1]*t.vec,P5[2]+V5[2]*t.vec))
  R<-sum(V1^2)/2/(M+1)
  grid<-expand.grid(-G[1]:G[1],-G[2]:G[2],-G[3]:G[3])*R
  penta3<-NULL
  L<-length(penta[,1])
  for (k in 1:L){
    penta3<-rbind(penta3,cbind(penta[k,1]+grid[,1],penta[k,2]+grid[,2],grid[,3]))
  }
  return(penta3)
}
#Euclidean Norm of a vector.
sq_norm<-function(v){
  return(sqrt(sum(v^2)))
}
#Set a vector to unit norm (Euclidean)
normalize<-function(v){
  return(v/(sq_norm(v)))
}
#Perform barycentric coordinate projection.
bary_proj<-function(pos,bary,P){
  D<-sqrt(diag(bary^(-1)))
  D[D==Inf]<-0
  V<-D%*%(P)%*%pos
  return(V)
}
#Compute Covariance Matrix for PCA
cov_mat<-function(data,mean){
  out<-(data-mean)%*%t(data-mean)
  return(out)
}
#Matrix squareroot transform for PCA.
sqrtm_transform<-function(A,B){
  if (!requireNamespace("expm", quietly = TRUE)) {
    stop("The package expm is required for this method. Please install it to use this function.")
  }
  tmp<-expm::sqrtm(B%*%A%*%B)
  return(tmp)
}

#Add Dummy point for augmented problem.
extend_positions<-function(pos){
  d<-dim(pos)[1]
  return(cbind(pos,rep(10^14,d)))
}
#Compute cost vor augmented problem.
newpgen_cost<-function(mat1,mat2,C){
  costM<-t(gen_cost(t(mat1),t(mat2)))
  costM[costM>10^14]<-C^2/2
  return(costM)
}
#Compute Euclidean Norm of a Vector
sqnorm<-function(x){
  return(sqrt(sum(x^2)))
}

#Build Matrix for MAAIPM
build_U<-function(N,m){
  rows<-c(1:(N-1),1:N)
  cols<-c(1:(N-1),rep(N,N))
  val<-rep(1,2*N-1)
  tmp1<-Matrix::sparseMatrix(rows,cols,x=val)
  tmp2<-Matrix::sparseMatrix(1:(m-1),1:(m-1),x=rep(1,m-1))
  return(kronecker(tmp1,tmp2))
}

#Build Constraint Matrix
build_const<-function(m,sizes){
  N<-length(sizes)
  sizes_csum<-c(0,cumsum(sizes))
  M<-sum(sizes)
  n_row<-M+N*(m-1)+1
  n_col<-(M+1)*m
  Mm<-M*m
  
  row1<-kronecker(1:M,rep(1,m))
  row2<-rep(0,M*(m-1))
  for (i in 1:N){
    index<-(sizes_csum[i]*(m-1)+1):(sizes_csum[i+1]*(m-1))
    row2[index]<-kronecker(rep(1,sizes[i]),(M+(i-1)*(m-1)+1):(M+i*(m-1)))
  }
  row3<-n_row*rep(1,m)
  row<-c(row1,row2,row3)
  
  #build col inds
  col1<-1:Mm
  col2<-(1:Mm)-kronecker((m*(1:M)+1-m),c(1,rep(0,m-1)))
  col2<-col2[col2>0]
  col3<-(n_col-m+1):n_col
  col<-c(col1,col2,col3)
  val<-rep(1,Mm+M*(m-1)+m)
  const<-Matrix::sparseMatrix(row,col,x=val)
  const[(M+1): (M+N*(m-1)), (n_col -m+2) : n_col]<-kronecker(rep(1,N),Matrix::sparseMatrix(1:(m-1),1:(m-1),x=rep(-1,m-1)))
  return(const)
}

#Normalise measure
probm<-function(x){
  return(x/sum(x))
}

#Refine grid for multiscale approach
refine_grid<-function(current.grid){
  current.size<-dim(current.grid)[1]
  ind.vec<-rep(1:current.size,rep(2,current.size))
  new.grid<-current.grid[,ind.vec]
  new.grid<-new.grid[ind.vec,]/4
  return(new.grid)
}
#convert data types for measures which are not probability measures.
process_data_unb<-function(object,type,return_type="weighted point pattern"){
  if (type=="image_file"){
    if (!requireNamespace("imager", quietly = TRUE)) {
      stop("The package imager is required to use images as input measures. Please install it to use this functionality.")
    }
    I<-imager::load.image(object)
    I<-imager::grayscale(I)
    d<-(attributes(I)$dim)[1:2]
    IM<-I[1:d[1],1:d[2],1,1]
    object<-IM
    type<-"image_mat"
  }
  if ((!(return_type==type))&(type=="image_mat")){
    d<-dim(object)
    G<-grid_positions(d[1],d[2])
    W<-as.vector(object)
    G<-G[W>0,]
    W<-W[W>0]
    object<-list(positions=G,weights=W)
    type<-"weighted point pattern"
  }
  if ((!(return_type==type))&(type=="wpp")){
    object<-list(positions=object$coordinates,weights=as.vector(object$mass))
    if (length(object$weights)==0){
      object<-list(positions=matrix(0,1,2),weights=10^-14)
    }
    type<-"weighted point pattern"
  }
  if ((!(return_type==type))&(type=="point pattern")){
    M<-dim(object[[1]])[1]
    object<-list(positions=object[[1]],weights=rep(1,M))
    type<-"weighted point pattern"
  }
  if (type=="weighted point pattern"){
    object<-list(positions=object$positions,weights=object$weights)
    return(object)
  }
  if (type==return_type){
    return(object)
  }
  return(list(positions=NaN,weights=NaN))
}

# altwpp<-function (coordinates, mass,dir="row") {
#   if (dir=="row"){
#     d<-dim(coordinates)[2]
#   }
#   if (dir=="col"){
#     d<-dim(coordinates)[1]
#   }
#   massnonzero<-mass>0
#   mass <- mass[massnonzero]
#   N<-length(mass)
#   if (N>1){
#     if (dir=="row"){
#       coordinates <- coordinates[massnonzero, ]
#     }
#     if (dir=="col"){
#       coordinates <- coordinates[ ,massnonzero]
#     }
#   }
#   if (N==1){
#     if (dir=="row"){
#       coordinates <- matrix(coordinates[massnonzero, ],1,d)
#     }
#     if (dir=="col"){
#       coordinates <- matrix(coordinates[massnonzero, ],d,1)
#     }
#   }
#   if (N==0){
#     if (is.null(d)){
#       d<-2
#     }
#     if (dir=="row"){
#       coordinates <- matrix(0,1,d)
#     }
#     if (dir=="col"){
#       coordinates <- matrix(0,d,1)
#     }
#   }
#   totmass<-sum(mass)
#   res <- list(dimension = d, N = N, coordinates = coordinates, 
#               mass = mass, totmass = totmass,dir=dir)
#   class(res) <- "wpp"
#   return(res)
# }

# print.wpp <- function(x, ...) {
#   stopifnot(class(x) == "wpp")
#   cat("Pattern of ",x$N," points in ",x$dimension, " dimensions with total mass ", x$totmass,".\n", sep="")
#   if (x$dir=="row"){
#     cat("Minimal coordinates:", apply(x$coordinates,2,min), "\n")
#     cat("Maximal coordinates:", apply(x$coordinates,2,max), "\n")
#   }
#   if (x$dir=="col"){
#     cat("Minimal coordinates:", apply(x$coordinates,1,min), "\n")
#     cat("Maximal coordinates:", apply(x$coordinates,1,max), "\n")
#   }
#   
#   invisible(x)
# }

#Sampling method for randomised UOT computations.
multinomial_sampling<-function(pos,weights,S){
  total_mass<-sum(weights)
  weights_normed<-weights/sum(weights)
  M<-length(weights)
  inds<-rmultinom(1,S,weights_normed)
  w.samp<-rep(total_mass/S,M)
  w.samp<-w.samp*inds
  pos.samp<-pos[,w.samp>0]
  wpp.out<-wpp(t(pos),w.samp)
  return(wpp.out)
}


#Normalise a vector with respect to the Euclidean norm
normalise<-function(v){
  return(v/(sq_norm(v)))
}

#Convert OT plans from data frame form to matrix form and vice versa.
plan_convert<-function(plan,target=NULL){
  if (target=="matrix"){
    if (is.data.frame(plan)){
      Tplan<-as.matrix(Matrix::sparseMatrix(i=plan$from,j=plan$to,x=plan$mass))
      return(Tplan)
    }
    else{
      return(plan)
    }
  }
  if (target=="frame"){
    if (is.data.frame(plan)){
      return(plan)
    }
    else{
      d<-dim(plan)
      rowInds<-matrix(0,d[1],d[2])
      colInds<-matrix(0,d[1],d[2])
      for (i in 1:(d[1])){
        rowInds[i,]<-rep(i,d[2])
      }
      for (i in 1:(d[2])){
        colInds[,i]<-rep(i,d[1])
      }
      from<-rowInds[plan>0]
      to<-colInds[plan>0]
      mass<-plan[plan>0]
      return(data.frame(from=from,to=to,mass=mass))
    }
  }
}

mplot_bar<-function(A,pos="top",col,fill,width,size){
  x<-NULL
  y<-NULL
  m<-NULL
  mdf<-data.frame(x=A$positions,m=A$weights)
  P<-ggplot2::ggplot(mdf)
  L<-length(A$weights)
  for (l in 1:L){
    P<-P+ggplot2::geom_line(data = data.frame(x = c(A$positions[l], A$positions[l]), y = c(0, A$weights[l])), ggplot2::aes(x = x, y = y),size=width,col=col) 
  }
  P<-P+ggplot2::geom_point(ggplot2::aes(x=x,y=m),size=size,show.legend = FALSE,col=fill)+ggplot2::theme_void()
  if (pos=="left"){
    P<-P+  ggplot2::scale_y_reverse(position = "right")+ggplot2::coord_flip()
  }
  return(P)
}


#plot a single marginal as a lineplot
mplot<-function(A,pos="top",linecol="black",fillcol="blue",width=3){
  x<-NULL
  y<-NULL
  m<-NULL
  mdf<-data.frame(x=A$positions,m=A$weights)
  P<-ggplot2::ggplot(mdf)+ggplot2::geom_area(size=width,color=linecol,fill=fillcol,ggplot2::aes(x=x,y=m))+ggplot2::theme_void()
  if (pos=="left"){
    P<-P+  ggplot2::scale_y_reverse(position = "right")+ggplot2::coord_flip()
  }
  return(P)
}
#Plot transport plan for marginal plots
Tplot<-function(A,B,plan,size,colgrad,type="points",lgrid=FALSE,tgrid=FALSE,grid.type="dashed",grid.col="lightgrey",grid.width=1){
  plan<-plan_convert(plan,target="frame")
  tdf<-data.frame(x=A$positions[plan$from,],y=B$positions[plan$to,],m=plan$mass)
  #Make the CRAN CHECK happy. 
  x<-NULL
  y<-NULL
  m<-NULL
  if (type=="pointsCol"){
    tplot<-ggplot2::ggplot(tdf,ggplot2::aes(x=x,y=y,col=m))
  }
  if (type=="pointsSize"){
    tplot<-ggplot2::ggplot(tdf,ggplot2::aes(x=x,y=y,size=m))
  }
  if (type=="line"){
    tplot<-ggplot2::ggplot(tdf,ggplot2::aes(x=x,y=y,col=m))
  }
  
  
  if (lgrid){
    y.vec<-B$positions
    for (y in y.vec){
      tplot<-tplot+ggplot2::geom_hline(yintercept=y,linetype=grid.type,col=grid.col,size=grid.width)
    }
  }
  if (tgrid){
    x.vec<-A$positions
    for (x in x.vec){
      tplot<-tplot+ggplot2::geom_vline(xintercept=x,linetype=grid.type,col=grid.col,size=grid.width)
    }
  }
  if (type=="pointsCol"){
    tplot<-tplot+ggplot2::geom_point(size=size[1],show.legend = FALSE)+
      ggplot2::scale_colour_continuous(low=colgrad[1],high=colgrad[2])+ggplot2::theme_void()
  }
  if (type=="pointsSize"){
    tplot<-tplot+ggplot2::geom_point(show.legend = FALSE,col=colgrad[1])+ggplot2::theme_void()+ 
      ggplot2::scale_size_continuous(range = size)
  }
  if (type=="line"){
    tplot<-tplot+ggplot2::geom_line(size=size[1],show.legend = FALSE)+
      ggplot2::scale_colour_continuous(low=colgrad[1],high=colgrad[2])+ggplot2::theme_void()
  }
  return(tplot)
}