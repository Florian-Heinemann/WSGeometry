
#' Cost Matrix
#' @description Calculates the cost matrix between points.
#'
#' @param x A numeric matrix. Each row corresponds to the coordinates of one point in the first point cloud.
#' @param y A numeric matrix. Each row corresponds to the coordinates of one point in the second point cloud
#' @param exp exponent to apply to the distance
#' @param wfr computing the wasserstein fisher rao distance
#' @param p parameter for the minkowski metric. standard p = 2 give the minkowski metric.
#' @return The distance matrix between the points. The rows correspond to the points in x, the columns to the
#' points in y
#' @noRd
cost_matrix <- function(x, y, exp = 1, p = 2, wfr = FALSE) {
    if (is.null(ncol(x))) {
        x <- matrix(x, ncol = 1)
    }
    
    if (is.null(ncol(y))) {
        y <- matrix(y, ncol = 1)
    }
    
    if (ncol(x) != ncol(y)) {
        print("Unequal dimensions.")
    }
    
    
    cMatrix <- matrix(rep(0, nrow(x) * nrow(y)), ncol = nrow(y))
    
    if (p == Inf) { 
        for (i in 1:nrow(x)) {
            cMatrix[i,] <- max(abs(t(t(y) - x[i,])))^exp
        }
    } else {
        for (i in 1:nrow(x)) {
            cMatrix[i,] <- ((rowSums((abs(t(t(y) - x[i,])))^p))^(1/p))^exp
        }
    }
    
    if (wfr) {
        cMatrix <- -2 * log(cospi(pmin(cMatrix / pi, 1 / 2)))
    }
    
    return(cMatrix)
}
