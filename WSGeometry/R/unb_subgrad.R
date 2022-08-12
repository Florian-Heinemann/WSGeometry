unb_subgrad_wrap<-function(pos.list,weights.list,support,C,frechet.weights=NULL,maxIter=100,startvec=NULL,stepsize=10^-1,thresh=10^-5,warmstart=TRUE,warmstartlength=100,threads=1){
d<-dim(pos.list[[1]])[2]
N<-length(pos.list)
total_masses<-(unlist(lapply(weights.list,sum)))
mass_cap<-max(total_masses)+2*mean(total_masses)
for (i in 1:N){
  weights.list[[i]]<-c(weights.list[[i]],mass_cap-total_masses[i])
}
weights.list<-lapply(weights.list,probm)
weights.list<-lapply(weights.list,matrix)
pos.list<-lapply(pos.list,extend_positions)
M<-dim(support)[2]
support<-extend_positions(support)
sizes<-unlist(lapply(weights.list,length))
costMat.list<-vector("list",N)
costMat.list<-lapply(pos.list,newpgen_cost,mat1=support,C=C)
if (is.null(frechet.weights)){
  frechet.weights=matrix(rep(1/N,N))
}
res<-krbary_subgrad_cxx(weights.list,costMat.list,frechet.weights,maxIter,stepsize,thresh,warmstart, warmstartlength,threads)
bary.w<-res[[1]][1:M]
bary.w[bary.w<10^(-14)]<-0
bary.w<-bary.w*mass_cap
return(bary.w)
#return(wpp(t(support[,1:M]),mass_cap*res[[1]][1:M]))

}