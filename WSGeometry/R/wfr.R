
#' @title Wasserstein-Fisher-Rao Distance
#' 
#' @description Approximate the Wasserstein-Fisher-Rao distance between two measures using entropic transport.
#' 
#' @details This function approximates the Wasserstein-Fisher-Rao distance between two measures using 
#' the Sinkhorn algorithm based on Corollary 5.9 in [1]. The theoretical exact result corresponds to the limit of \code{epsi} going to 0.
#'
#' @references 
#' [1] Lenaic Chizat, Gabriel Peyré, Bernhard Schmitzer, François-Xavier Vialard (2019). “Unbalanced Optimal Transport: Dynamic and Kantorovich Formulation” arXiv: 1508.05216v3.
#' 
#'
#' @param supplyList A list containing the supply measure and the support of the supply distribution (see example).
#' @param demandList A list containing the demand measure and the support of the demand distribution (see example).
#' @param delta Transport range parameter.
#' @param epsi (optional) Final regularization strength used for approximation (default: \code{1e-8}).
#' @param maxIteration (optional) Maximal number of iterations to perform (default: \code{5000}).
#' @param tol (optional) convergence tolerance for Sinkhorn algorithm (default: \code{1e-3}).
#' 
#' @return A single number: the entropic WFR distance.
#' 
#' @examples
#' supplyPoints <- seq(0, 1, length.out=100)
#' demandPoints <- seq(0, 1, length.out=100)
#' 
#' supply <- rep(0, 100)
#' supply[2:20] = 2
#' supply[91:100] = c(seq(0, 2, length.out=5), seq(2, 0, length.out=5))
#' 
#' demand <- rep(0, 100)
#' demand[21:40] = seq(0, 2, length.out=20)
#' demand[51:90] = 1.3 * sqrt(1 - seq(-1, 1, length.out=40)^2)
#' 
#' supplyList <- list(supply, supplyPoints)
#' demandList <- list(demand, demandPoints)
#' 
#' epsi <- 1e-3
#' delta <- 1
#' 
#' wfr(supplyList, demandList, delta, epsi)
#' 
#' @export
wfr <- function(supplyList, demandList, delta, epsi = 1e-8, maxIteration = 5000, tol = 1e-3) {
  
  costMatrix <- cost_matrix(supplyList[[2]], demandList[[2]], 1, 2, TRUE, delta)
  
  epsVector <- exp(seq(0, log(epsi), length.out = max(1, ceiling(-log(epsi)))))
  div <- list("KL", 1)
  
  res <- regularized_transport(supplyList, demandList, supplyDivList = div, demandDivList = div, 
                               epsVector = epsVector, maxIteration = maxIteration, scalingTol = tol, finalTol = tol,
                               costMatrix = costMatrix, algorithm = "sinkhorn")
  
  if (res$primalDualGap > 1e-2)
    print(sprintf("Warning: large primal dual gap: %f", res$primalDualGap))
  
  cost <- 0.5 * (res$primalCost + res$dualCost)
  
  return(sqrt(2 * delta**2 * cost))
}


