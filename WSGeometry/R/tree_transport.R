#' Transport on Trees
#' 
#' Solving unbalanced optimal transport problems on trees.
#' 
#' This function can be used to solve unbalanced optimal transport problems that use a tree metric as cost function.
#' If the cost matrix is derived from a tree metric, the unbalanced optimal transport problem  \eqn{\min_r <C,r> + \sum_i p_i(\alpha_i - \sum_j r_{ij}) + \sum_j q_j(\beta-\sum_i r_{ij})}{
#' min_r <C,r> + sum_i(p_i (a-sum_j r_ij)) + sum_j(q_j (b-sum_i r_ij))} with supply and demand measure \eqn{\alpha}{a} and \eqn{\beta}{b}
#' , transport plan \eqn{r}, construction and destruction costs \eqn{p} and \eqn{q} and cost matrix \eqn{C} can be solved efficiently in
#' \eqn{O(n \log^2(n))}{O(n log²(n))} time. The basic C++ implementation of the algorithm used in this package can be found at \url{https://github.com/joisino/treegkr}. It was extended in order to make it accessible 
#' from R and compute the transport map. 
#'
#' @references
#' Sato R, Yamada M, Kashima H (2021). “Fast Unbalanced Optimal Transport on a Tree.” arXiv: 2006.02703.
#'
#'
#' @param tree A tree structure in list format:
#' The first element is the index of the root node.
#' The other elements are vectors defining the edges of the tree. Each of these vectors has to be 
#'  of the form \code{c(parentNodeIndex, childNodeIndex, edgeWeight)}.
#' @param supplyList A list containing two vectors: the amount of supply and the cost of destruction / export per node.
#' @param demandList A list containing two vectors: the amount of demand and the cost of construction / import per node.
#' @param output Determines the format of the output (Default: \code{"cost"}):
#' \itemize{
#' \item \code{"transportPlan"} give a standard transport plan matrix where the value at \eqn{i,j} gives the amount of mass transported from 
#' note \eqn{i} to \eqn{j}. Note that this requires \eqn{O(n^2)}{O(n²)} space and time.
#' \item \code{"list"} gives a list of vectors. Each vector holds three entries: The first gives the supply node, the second the target node of the transport
#' and the third the amount of mass transported between those nodes.
#' \item \code{"cost"} returns only the transport cost.
#' }
#'
#' @return A list containing the transport cost, the mass import and export vectors, and the transport list or plan if specified in \code{output}.
#'
#'
#' @examples
#'tree <- list(1, c(1,2,1), c(2,3,1), c(2,4,1), c(1,5,1), c(5,6,1),
#'                c(5,7,1), c(1,8,10), c(8,9,1), c(8,10,1))
#'
#'constructionCost <- rep(2, 10) 
#'destructionCost <- rep(2, 10)
#'
#'supply <- c(0,0,2,0,0,4,0,0,2,0)
#'demand <- c(0,0,0,0,5,0,0,0,0,1)
#'
#'supplyList = list(supply, destructionCost)
#'demandList = list(demand, constructionCost)
#'
#'transport <- tree_transport(tree, supplyList, demandList,
#'                           output = "list")
#'
#'plot_tree(tree, tList = transport$transportList,
#'          supply = supply, demand = demand)
#'
#'
#'tree <- list(1, c(1,2,1), c(2,3,1), c(3,4,1), c(3,5,1), c(2,6,1),
#'                c(1,7,1), c(7,8,1), c(7,9,1), c(9,10,1), c(9,11,1),
#'                c(11,12,1), c(11,13,1))
#'
#'
#'
#' constructionCost <- rep(2, 13)
#' destructionCost <- rep(1, 13)
#'
#'
#' supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
#' demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)
#' 
#' supplyList = list(supply, destructionCost)
#' demandList = list(demand, constructionCost)
#'
#' transport <- tree_transport(tree, supplyList, demandList,
#'                            output = "list")
#'
#' plot_tree(tree, tList = transport$transportList,
#'           supply = supply, demand = demand)
#'
#'
#' @export
tree_transport <- function(tree, supplyList, demandList, output = "cost") {

    supply = supplyList[[1]]
    demand = demandList[[1]]
    
    creationCost = supplyList[[2]]
    destructionCost = demandList[[2]]
    

    if (output != "list" & output != "transportPlan"  & output != "cost") {

        print("Unknown output setting. Output was set to 'cost'")
        output <- "cost"

    }

    treegkrOut <- treegkr_Rcpp(tree[-1], supply, demand, creationCost, destructionCost)

    importVec = treegkrOut$import
    transportList = treegkrOut$transportList
    
    # handle mass transport in place if supply
    for (i in 1:length(supply)) {
        m = min(supply[i], demand[i])
        if (m > 0) {
            transportList$to[length(transportList$to) + 1] = i
            transportList$from[length(transportList$from) + 1] = i
            transportList$mass[length(transportList$mass) + 1] = m
        }
    }
    
    # importVec hat positive and negative entries. The positive entries indicate
    # imported mass, the negative entries exported mass.
    # Splitting the vector in import and export
    import <- pmax(importVec, 0)
    export <- pmax(-importVec, 0)
    
    if (output == "cost") {
        # return only cost, import and export
        result <- list(treegkrOut$cost, import, export)
        names(result) <- c("cost", "import", "export")
        
        return(result)
    } else {
        res <- transportList
        
        # Computing the transport plan in list form. Each entry consists of the node the mass
        # comes from, the node the mass is transported to and the amount of mass.
        transportList <- list()
        if (length(res$mass) > 0) {
            for (i in 1:length(res$mass)) {
                if (res$mass[i] != 0) {
                    transportList[[length(transportList)+1]] <- c(res$from[i], res$to[i], res$mass[i])
                }
            }
        }
        
        # return a list according to the "output" variable, calculating matrix if necessary
        if (output == "list") {
            result <- list(treegkrOut$cost, transportList, import, export)
            names(result) <- c("cost", "transportList", "import", "export")
            return(result)
        } else if (output == "transportPlan") {
            tPlan <- matrix(0, length(supply), length(demand))
            if (length(res) > 0) {
                tPlan[cbind(res$from, res$to)] <- res$mass
            } 
            result <- list(treegkrOut$cost, tPlan, import, export)
            names(result) <- c("cost", "transportPlan", "import", "export")
            
            return(result)
        }
    }
}



