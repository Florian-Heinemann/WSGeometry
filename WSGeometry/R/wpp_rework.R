#' Constructor function for \link[transport]{wpp-object}s.
#' @description Generates a \link[transport]{wpp-object} from given coordinates and weights.  
#' @param coordinates A Mxd matrix where each row represents one location of a support point of the point pattern.
#' @param mass A vector of length M of non-negative values specifying the mass at the given locations.
#' @return A \link[transport]{wpp-object} corresponding to the given mass and coordinates.
#' @export
wpp<-function (coordinates, mass) {
  d<-dim(coordinates)[2] 
  massnonzero<-mass>0
  mass <- mass[massnonzero]
  N<-length(mass)
  if (N>1){
    zeropoints <- coordinates[!massnonzero, ]
    coordinates <- coordinates[massnonzero, ]
  }
  if (N==1){
    coordinates <- matrix(coordinates[massnonzero, ],1,d)
  }
  if (N==0){
    if (is.null(d)){
      d<-2
    }
    coordinates <- matrix(0,1,d)
  }
  totmass<-sum(mass)
  res <- list(dimension = d, N = N, coordinates = coordinates, 
              mass = mass, totmass = totmass,dir=dir)
  class(res) <- "wpp"
  if (N>1){
    attr(res, "zeropoints") <- zeropoints 
  }
  return(res)
}

