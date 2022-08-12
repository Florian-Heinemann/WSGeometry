#' Compute the Frechet functional/The objective value of the KR/(p,C)-barycenter problem
#' @description This function computes the objective value of the KR/(p,C)-barycenter problem for a given measure and a given dataset of measures. 
#' @param bary An object representing a measure, for which the Frechet value should be computed. Should be one of the following: 
#' A matrix, representing an image; A path to a file containing an image; 
#' A \link[transport]{wpp-object}; 
#' A \link[transport]{pp-object};
#' A list containing an entry named `positions` with the support of the measure and an entry named `weights` containing the weights of the support points;
#' A list containing en entry named `positions`` specifying the support of a measure with uniform weights.  
#' @param data A list of objects which should be compared to bary. Each element should be one of the following:
#' A matrix, representing an image; A path to a file containing an image; 
#' A \link[transport]{wpp-object}; 
#' A \link[transport]{pp-object};
#' A list containing an entry named `positions` with the support of the measure and an entry named `weights` containing the weights of the support points;
#' A list containing en entry named `positions`` specifying the support of a measure with uniform weights.  
#' @param C A double specifying the penalty for the (p,C)-Kantorovich-Rubinstein distance.
#' @return A real number specifying the Frechet value of the input object for the given dataset.
#' @export
krfrechet_func<-function(bary,data,C){
  type<-type_check(bary)
  bary<-process_data_unb(bary,type)
  bary<-wpp(bary$positions,bary$weights)
  data.types<-lapply(data,type_check)
  data<-mapply(process_data_unb,data,data.types,SIMPLIFY = FALSE)
  data<-mapply(wpp,lapply(data,"[[",1),lapply(data,"[[",2),SIMPLIFY = FALSE)
  N<-length(data)
  val<-0
  for (k in 1:N){
    val<-val+(kr_dist(bary,data[[k]],p=2,C)$distance^2)
  }
  return(val/N)
}