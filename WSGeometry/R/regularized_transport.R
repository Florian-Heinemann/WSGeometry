




#' @title Regularized Unbalanced Optimal Transport
#' 
#' @description  Solving regularized unbalanced optimal transport problems.
#' 
#' 
#' @details  This function solves the regularized unbalanced optimal transport problem using either the Scaling 
#' or Sinkhorn algorithm described in the papers given below.
#' The algorithms minimize the regularized unbalanced optimal transport problem given by
#' 
#' \eqn{\min_{b \in \!R^{X \times Y}_+} F_1(P_X b) + F_2(P_Y b) + \varepsilon KL(b|K)}{min_b F_1(P_X b) +  F_2(P_Y b) + epsilon * KL(b|K)}
#' 
#' where \eqn{F_1}{F_1} and \eqn{F_2}{F_2} are divergence functions, \eqn{P_X b} and \eqn{P_Y b} the marginals of \eqn{b}, 
#' \eqn{KL(\cdot | \cdot)}{KL( | )} the Kullback Leibner divergence and \eqn{K = \exp(-c(x,y)/eps)}{K = exp(-c(x,y)/eps)} the kernel associated with the cost matrix.
#' 
#' 
#' 
#' The following divergence functions are currently implemented: 
#' \itemize{
#'   \item Kullback-Leibner divergence ("KL")
#'   \item Total variation ("TV")
#'   \item The divergence associated with the range constraint ("RG")
#'   \item The Power divergence ("Power") and its special cases the Hellinger ("Hellinger") and Berg ("Berg") divergence. (Only available for the Sinkhorn algorithm.)
#' }
#' 
#' Instead of providing a cost matrix to the function, it is possible to provide the support of the supply and demand measures and 
#' have the function compute the cost matrix itself. 
#' 
#' 
#' @references
#' 
#' Chizat L, Peyré G, Schmitzer B, Vialard F (2018). “Scaling algorithms for unbalanced optimal transport problems.” Mathematics of Computation, 87(314), 2563–2609.\cr
#' Schmitzer B (2019). “Stabilized sparse scaling algorithms for entropy regularized transport problems.” SIAM Journal on Scientific Computing, 41(3), A1443–A1481.\cr
#' Séjourné T, Feydy J, Vialard F, Trouvé A, Peyré G (2021). “Sinkhorn Divergences for Unbalanced Optimal Transport.” arXiv: 1910.12958. 1910.12958v2.
#' 
#' 
#'
#' @param supplyList A list containing the supply measure and, if the cost matrix is not provided, the support of the supply distribution (see example).
#' @param demandList A list containing the demand measure and, if the cost matrix is not provided, the support of the demand distribution (see example).
#' 
#' @param supplyDivList (optional) A list containing the information about the divergence function for the supply measure. The first element is the
#' abbreviation of the divergence function. The following values are available: 
#' 
#' \itemize{
#'   \item \code{"KL"} for the Kullback-Leibner divergence \eqn{F = \lambda \cdot KL()}{F1 = lambda * KL()}
#'   \item \code{"TV"} for total variation divergence \eqn{F = \lambda \cdot TV()}{F1 = lambda * TV()}
#'   \item \code{"RG"} for the range constraint 
#'   \item \code{"Power"} for the Power divergence with conjugate exponent r: \eqn{F_1 = \lambda \cdot Power_{r}()}{F1 = lambda * Power_r()}
#'   \item \code{"Hellinger"} for the Power divergence with conjugate exponent r = -1: \eqn{F_1 = \lambda \cdot Power_{-1}()}{F1 = lambda * Power_-1()}
#'   \item \code{"Berg"} for the Power divergence with conjugate exponent r = 0: \eqn{F_1 = \lambda \cdot Power_0()}{F1 = lambda * Power_0()}
#'   \item \code{"LSR"} for a range constraint with slopes at the boundaries rather than hard limits
#' }
#' 
#' The Power divergence and its special cases are only implemented for the Sinkhorn algorithm. 
#' 
#' The next elements in supplyDivList are the additional parameters for the divergence functions. For \code{"KL"}, \code{"TV"}, \code{"Hellinger"}
#' and \code{"Berg"} the only required parameter is the regularization parameter \eqn{\lambda}{lambda}. The \code{"RG"} divergence function requires two parameters
#' alpha and beta with \eqn{0 \leq \alpha \leq \beta}{0 <= alpha <= beta} that define the lower and upper bounds. Lastly, \code{"Power"} requires
#' two parameters as well. The first is the regularization parameter \eqn{\lambda}{lambda}, followed by the conjugate exponent \eqn{r}.  
#' The default value is \code{list("KL", 0.5)}.
#' 
#' @param demandDivList (optional) A list containing the information about the divergence
#' function for the demand measure in the same form as the supplyDivList. The default value is \code{list("KL", 0.5)}.
#' @param epsVector (optional) A numeric value or vector holding the regularization parameter. If a vector is given, epsilon scaling will be performed. 
#' This may be required for convergence if epsilon is small. The default value is \code{1e-3}.
#' @param maxIteration (optional) The maximum number of algorithm iterations. The default value is \code{5000}.
#' @param tol (optional) A numeric value. If the change of the dual variables from one step to the next is smaller than this value the algorithm is 
#' considered to be converged. If the primal-dual-gap returned is very large even though the maximal number of iterations is not reached, this value may be too large. 
#' The default value is \code{1e-7}.
#' @param exp (optional) The exponent applied to the cost function. Default value is \code{1}.
#' @param p (optional) The parameter for calculating the \eqn{L_p}{L_p} cost. Can be a positive real number or Inf. Default value is 2 calculating the euclidean distance.
#' @param costMatrix (optinal) A cost matrix for transport between the supply and demand points. 
#' Default value is \code{NULL} for which the cost matrix is calculated based on \code{supplyList[[2]]}, \code{demandList[[2]]}, \code{exp}, \code{p} and \code{wfr}.
#' @param algorithm (optinal) Set to \code{"sinkhorn"} to use the Sinkhorn algorithm. Default value is \code{"scaling"}. Only 
#' these two algorithms are available. Usually \code{"scaling"} is faster but may be less numerically stable, especially for small epsilon.
#' Also note that \code{"sinkhorn"} will not create / destroy mass where the supply / demand is zero.
#' @param indicatorSlack (optional) amount of violation to allow in indicator functions before returning infinity. 
#' This is needed for numeric stability and due to finite number of iterations (default is \code{1e-6}).
#' @param threadCount (optional) number of threads to use for Sinkhorn algorithm (default: 1 if problem is small, maximal number available if problem is large); 
#' the scaling algorithm uses the matrix-vector multiplication builtin in armadillo, which can be parallel depending on the installation.
#' @param maxLinesPerWork (optional) how many elements to assign to one thread invocation in Sinkhorn, tweaking this can improve performance (default: \code{5}).
#' 
#' @return A list of named values:
#' \itemize{
#'   \item \code{"transportPlan"} the primal solution
#'   \item \code{"primalCost"} the primal cost
#'   \item \code{"dualCost"} the dual cost
#'   \item \code{"primalDualGap"} the difference between primal and dual costs
#'   \item \code{"dualSupply"} the dual variable at the supply marginal
#'   \item \code{"dualDemand"} the dual variable at the demand marginal
#'   \item \code{"finalStepSize"} the final step size
#'   \item \code{"iterations"} the total numer of steps performed
#'   \item \code{"epsIterations"} the numer of iterations performed per value of epsilon in \code{epsVector}
#' }
#'   
#' @examples 
#' I <- 1000
#' J <- 1000
#' X <- seq(0, 1, length.out = I)
#' Y <- seq(0, 1, length.out = J)
#' 
#' p <- rep(0, I)
#' p[2:200] = 2
#' p[901:1000] = c(seq(0, 2, length.out = 50), seq(2, 0, length.out = 50))
#' 
#' q <- rep(0, J)
#' q[201:400] = seq(0, 2, length.out = 200)
#' q[501:900] = 1.3 * sqrt(1 - seq(-1, 1, length.out = 400)^2)
#' 
#' supply <- list(p, X)
#' demand <- list(q, Y)
#' 
#' maxIter <- 100
#' eps <- 1e-3
#' 
#' supplyDiv <- list("KL", 0.04)
#' demandDiv <- list("KL", 0.04)
#' 
#' res <- regularized_transport(supply, demand, supplyDiv, demandDiv,
#'                             maxIteration = maxIter, epsVector = eps, exp = 2)
#' plot_1D_transport(res$transportPlan, supply, demand)
#'
#' @export
regularized_transport <- function(supplyList, demandList, supplyDivList = list("KL", 0.5), demandDivList = list("KL", 0.5), 
                                  epsVector = 1e-3, maxIteration = 5000, tol = 1e-7, exp = 1, p = 2,
                                  costMatrix = NULL, algorithm = "scaling", indicatorSlack = 1e-6,
                                  threadCount = -1, maxLinesPerWork = 5) {
    
    supplyDivType <- NULL
    supplyDivParameters <- NULL
    demandDivType <- NULL
    demandDivParameters <- NULL
    
    if (any(epsVector <= 0)) 
        stop("epsilon has to be positive")
    
    
    if (supplyDivList[[1]] == "KL") {
        supplyDivType <- 1
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "TV") {
        supplyDivType <- 2
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "RG") {
        supplyDivType <- 3
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "Berg") {
        supplyDivType <- 4
        supplyDivList[3] = 0
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "Hellinger") {
        supplyDivType <- 4
        supplyDivList[3] = 0.5
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "Power") {
        supplyDivType <- 4
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else if (supplyDivList[[1]] == "LSR") {
        supplyDivType <- 5
        supplyDivParameters <- unlist(supplyDivList[-1])
    } else {
        stop("Please supply a valid divergence")
    }
    
    if (demandDivList[[1]] == "KL") {
        demandDivType <- 1
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "TV") {
        demandDivType <- 2
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "RG") {
        demandDivType <- 3
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "Berg") {
        demandDivType <- 4
        demandDivList[3] = 0
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "Hellinger") {
        demandDivType <- 4
        demandDivList[3] = 0.5
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "Power") {
        demandDivType <- 4
        demandDivParameters <- unlist(demandDivList[-1])
    } else if (demandDivList[[1]] == "LSR") {
        demandDivType <- 5
        demandDivParameters <- unlist(demandDivList[-1])
    } else {
        stop("Please demand a valid divergence")
    }
    
    
    supply <- supplyList[[1]]
    demand <- demandList[[1]]
    
    
    if (is.null(costMatrix)) {
        costMatrix <- cost_matrix(supplyList[[2]], demandList[[2]], exp, p)
    }  
    
    
    if (algorithm == "sinkhorn") {
        
        res <- Sinkhorn_Rcpp(cost_matrix = costMatrix, supply = supply, demand = demand,
                             supply_div_type = supplyDivType, demand_div_type = demandDivType, 
                             supply_div_parameters = supplyDivParameters, demand_div_parameters = demandDivParameters,
                             iter_max = maxIteration, epsvec = epsVector, tol = tol, indicator_slack = indicatorSlack,
                             thread_cnt = threadCount, max_lines_per_work = maxLinesPerWork)
        
    } else if (algorithm == "scaling") {
        
        res <- StabilizedScaling_Rcpp(cost_matrix = costMatrix, supply = supply, demand = demand,
                                      supply_div_type = supplyDivType, demand_div_type = demandDivType, 
                                      supply_div_parameters = supplyDivParameters, demand_div_parameters = demandDivParameters,
                                      iter_max = maxIteration, epsvec = epsVector, tol = tol, indicator_slack = indicatorSlack)
        
    }
    
    returnList <- list("transportPlan" = res$transportPlan,
                       "primalCost" = res$primalCost,
                       "dualCost" = res$dualCost,
                       "primalDualGap" = abs(res$primalCost - res$dualCost),
                       "dualSupply" = res$supplyDual[1,1:length(res$supplyDual)],
                       "dualDemand" = res$demandDual[1,1:length(res$demandDual)],
                       "finalStepSize" = res$finalStepSize,
                       "iterations" = res$iterations,
                       "epsiIterations" = res$epsiIterations[1,1:length(res$epsiIterations)])
    return(returnList)
}