#' Split the values of entries in a matrix between a specified area round it. 
#' @description Takes a matrix M and splits the value of the matrix at a given coordinate (i,j) with a rectangle of 
#' positions around it given by r1 and r2. The position (i,j) will get its previous value divided by (r1 x r2) and the 
#' surrounding positions (r1 in horizontal and r2 in vertical direction) will have their entries increased by the same value.
#' @param M A matrix with real numbers as entries.
#' @param r1 Integer specifying the range of the split in the horizontal direction.
#' @param r2 Integer specifying the range of the split in the vertical direction.
#' @return A matrix of the same dimensions as M, which had the mass split applied to all entries simultaneously.
#' @export
smear<-function(M,r1,r2){
  d<-dim(M)
  M.out<-M
  M.out[M.out>0]<-0
  for (i in 1:d[1]){
    for (j in 1:d[2]){
      if (M[i,j]>0){
        M.out[seq(max(1,i-r1),min(i+r1,d[1])),seq(max(1,j-r2),min(j+r2,d[2]))]<-M.out[seq(max(1,i-r1),min(i+r1,d[1])),seq(max(1,j-r2),min(j+r2,d[2]))]+M[i,j]
      }
    }
  }
  return(M.out)
}
