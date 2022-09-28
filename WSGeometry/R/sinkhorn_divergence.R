

#' @title Sinkhorn Divergence
#' 
#' @description Calculate the entropic Sinkhorn divergence for unbalanced transport
#' 
#' @details This function calculates the entropic Sinkhorn divergence for unbalanced transport defined as 
#' \eqn{OT_{\varepsilon}(\alpha, \beta) - \frac{1}{2} OT_{\varepsilon}(\alpha, \alpha) - \frac{1}{2} OT_{\varepsilon}(\beta, \beta) + \frac{\varepsilon}{2}(m(\alpha) - m(\beta))^2}{OT(a, b) - 0.5 OT(a, a) - 0.5 OT(b, b) - 0.5 (m(a) - m(b))^2}
#' where \eqn{m(\alpha)}{m(a)} is the total mass of \eqn{\alpha}{a}.
#'
#' @references 
#' Séjourné T, Feydy J, Vialard F, Trouvé A, Peyré G (2021). “Sinkhorn Divergences for Unbalanced Optimal Transport.” arXiv: 1910.12958. 1910.12958v2.
#' 
#' @param supplyList A list containing the supply measure and the support of the supply distribution (see \link[WSGeometry]{regularized_transport}).
#' @param demandList A list containing the demand measure and the support of the demand distribution (see \link[WSGeometry]{regularized_transport}).
#' @param supplyDivList A list describing the supply divergence (see \link[WSGeometry]{regularized_transport}).
#' @param demandDivList A list describing the demand divergence (see \link[WSGeometry]{regularized_transport}).
#' @param epsi (optional) Final regularization strength used for approximation (default: \code{1e-8}).
#' @param ssCostMatrix (optional) Cost matrix for transport between supply and itself, if \code{NULL} it is calculated from \code{supplyList[[2]]} (default: \code{NULL}).
#' @param sdCostMatrix (optional) Cost matrix for transport between supply and demand, if \code{NULL} it is calculated from \code{supplyList[[2]]} and \code{demandList[[2]]} (default: \code{NULL}).
#' @param ddCostMatrix (optional) Cost matrix for transport between demand and itself, if \code{NULL} it is calculated from \code{demandList[[2]]} (default: \code{NULL}).
#' @param p (optional) Norm to use for cost matrix calculation (default: 2).
#' @param exp (optional) Power to apply to calculated cost matrix (default: 1).
#' @param maxIteration (optional) Maximal number of iterations to perform (default: \code{5000}).
#' @param tol (optional) convergence tolerance for Sinkhorn algorithm (default: \code{1e-3}).
#' 
#' @return A single number: the calculated Sinkhorn divergence
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
#' supplyDiv <- list("TV", 1, 1)
#' demandDiv <- list("TV", 1, 1)
#' 
#' epsi <- 1e-3
#' 
#' sinkhorn_divergence(supplyList, demandList, supplyDiv, demandDiv, 
#'                     epsi, maxIteration=10000, tol=1e-4)
#' 
#' @export
sinkhorn_divergence <- function(supplyList, demandList, supplyDivList, demandDivList, epsi = 1e-8, 
                                p = 2, exp = 1, ssCostMatrix = NULL, sdCostMatrix = NULL, ddCostMatrix = NULL,
                                maxIteration = 5000, tol = 1e-3) {
  
  epsVector <- exp(seq(0, log(epsi), length.out = max(1, ceiling(-log(epsi)))))
  
  if (is.null(ssCostMatrix))
    ssCostMatrix <- cost_matrix(supplyList[[2]], supplyList[[2]], p = p, exp = exp)
  if (is.null(sdCostMatrix))
    sdCostMatrix <- cost_matrix(supplyList[[2]], demandList[[2]], p = p, exp = exp)
  if (is.null(ddCostMatrix))
    ddCostMatrix <- cost_matrix(demandList[[2]], demandList[[2]], p = p, exp = exp)
  
  res_sd <- regularized_transport(supplyList, demandList, supplyDivList, demandDivList, 
                                  epsVector = epsVector, costMatrix = sdCostMatrix,
                                  maxIteration = maxIteration, scalingTol = tol, finalTol = tol)
  res_ss <- regularized_transport(supplyList, supplyList, supplyDivList, supplyDivList, 
                                  epsVector = epsVector, costMatrix = ssCostMatrix,
                                  maxIteration = maxIteration, scalingTol = tol, finalTol = tol)
  res_dd <- regularized_transport(demandList, demandList, demandDivList, demandDivList, 
                                  epsVector = epsVector, costMatrix = ddCostMatrix,
                                  maxIteration = maxIteration, scalingTol = tol, finalTol = tol)
  
  if (res_sd$primalDualGap > 1e-2)
    print(sprintf("Warning: large primal-dual gap: %f", res_sd$primalDualGap))
  if (res_ss$primalDualGap > 1e-2)
    print(sprintf("Warning: large primal-dual gap: %f", res_ss$primalDualGap))
  if (res_dd$primalDualGap > 1e-2)
    print(sprintf("Warning: large primal-dual gap: %f", res_dd$primalDualGap))
  
  cost_sd <- 0.5 * (res_sd$primalCost + res_sd$dualCost)
  cost_ss <- 0.5 * (res_ss$primalCost + res_ss$dualCost)
  cost_dd <- 0.5 * (res_dd$primalCost + res_dd$dualCost)
  
  return(cost_sd - 0.5 * cost_ss - 0.5 * cost_dd + 0.5 * epsi * (sum(supplyList[[1]]) - sum(demandList[[1]]))**2)
}