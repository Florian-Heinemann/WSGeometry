#' Plotting dense 1D transport
#' 
#' A function to plot optimal transport between discretizations of continuous supply and demand measures in 1D. 
#' 
#' 
#' @param transportPlan A non negative numeric matrix giving the mass transport. The value at \eqn{(i, j)}
#' gives the mass transported from supply point \eqn{i} to demand point \eqn{j}.
#' @param supplyList A list containing the non negative supply measure and the underlying discretization as vectors.
#' @param demandList A list containing the non negative demand measure and the underlying discretization as vectors.
#' @param stacked A boolean declaring whether to draw the transported intervals vertically stacked (Default: \code{FALSE}).
#' @param mirrored A boolean declaring whether to draw supply and demand in opposite directions (Default: \code{FALSE}).
#' @param palette The color palette to use for transport visualization, may be "rainbow" or any accepted by \link[grDevices]{hcl.colors} (Default: "spectral").
#'
#' \if{html}{\figure{1D.png}}
#' \if{latex}{\figure{1D.png}{options: width=0.5in}}
#'
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
#' maxIter <- 200
#' eps <- 1e-3
#' 
#' suppyDiv <- list("KL", 0.04)
#' demandDiv <- list("KL", 0.04)
#' res <- regularized_transport(supply, demand, suppyDiv, demandDiv,
#'                              maxIteration = maxIter, epsVector = eps, p = 2)
#' plot_1D_transport(res$transportPlan, supply, demand)
#' 
#'
#' @export
plot_1D_transport <- function(transportPlan, supplyList, demandList, stacked = TRUE, mirrored = FALSE, palette = "spectral") {
    
    # Number of color intervals
    numIntervals <- 50
    
    X <- supplyList[[2]]
    
    x1Measure <- rep(1, length(supplyList[[1]]))
    y1Measure <- rep(1, length(demandList[[1]]))
    
    if (palette == "rainbow") {
        colors <- rainbow(numIntervals, s = 1, v = 1, start = 0, 
                          end = max(1, numIntervals - 1) / numIntervals, alpha = 1, rev = FALSE)   
    }
    else {
        colors <- hcl.colors(numIntervals, palette)
    }
    

    
    ymax <- max(supplyList[[1]], demandList[[1]], t(transportPlan) %*% x1Measure, transportPlan %*% y1Measure)
    
    supplyMax = max(supplyList[[1]], transportPlan %*% y1Measure)
    demandMax = max(demandList[[1]], t(transportPlan) %*% x1Measure)
    
    dir = if (mirrored) -1 else 1
    
    # Initialize plot
    if (mirrored)
        plot(supplyList[[2]], supplyList[[1]], type = "n", ylim= c(-demandMax, supplyMax), xlab = "Position", ylab = "Mass")
    else
        plot(supplyList[[2]], supplyList[[1]], type = "n", ylim= c(0, ymax), xlab = "Position", ylab = "Mass")
    
    firstSupp <- supplyList[[2]][1]
    firstDem <- demandList[[2]][1]
    
    if (stacked) {
        supBaseLevel <- rep(0, length(X))
        demBaseLevel <- rep(0, length(X))
    }
    
    # Plot the intervals for different colors
    for (i in 1:numIntervals) {
        # Calculating the interval of X
        colSuppInter <- rep(1, length(X))
        colSuppInter[(i - 1) / numIntervals > X | i / numIntervals < X] <- 0
        
        # Restricting the transport map on the interval
        subK <- transportPlan * colSuppInter
        
        # Adding the color intervals
        if (stacked) {
            sup = rowSums(subK)
            dem = colSums(subK)
            
            polygon(c(demandList[[2]], rev(demandList[[2]])), dir * c(demBaseLevel, rev(demBaseLevel + dem)), col = colors[i])
            polygon(c(supplyList[[2]], rev(supplyList[[2]])), c(supBaseLevel, rev(supBaseLevel + sup)), col = colors[i])
            
            supBaseLevel <- supBaseLevel + sup
            demBaseLevel <- demBaseLevel + dem
        } else {
            polygon(c(firstDem, demandList[[2]]), dir * c(0, colSums(subK)), col = colors[i])
            polygon(c(firstSupp, supplyList[[2]]), c(0, rowSums(subK)) , col = colors[i])
        }
    }
    
    # plot outlines of transported measures
    lines(demandList[[2]], dir * t(transportPlan) %*% x1Measure, type = "l", col = "green")
    lines(supplyList[[2]], transportPlan %*% y1Measure, type = "l", col = "blue")
    
    # plot outlines of input measures
    lines(supplyList[[2]], supplyList[[1]], type = "l", lty = 3 , col = "blue")
    lines(demandList[[2]], dir * demandList[[1]], type = "l", lty = 3, col = "green")
    
    # plot ground line
    lines(supplyList[[2]], rep(0, length(supplyList[[2]])), type = "l", col = "black")
}


#' Next Layer
#'
#' This recursively calculates the coordinates descendents of \code{node}
#'
#' @param treeDF A tree in data.frame format
#' @param coordinates The coordinates of the tree nodes as data.frame
#' @param node The current node index
#' @param layer The current layer
#'
#' @return data.frame with coordinates and information for all tree nodes
#'
#' @noRd
next_layer <- function(treeDF, coordinates, node, layer) {
    
    children <- treeDF[treeDF$parent == node,]$child
    numChildren <- length(children)
    
    if (numChildren == 0) {
        return(coordinates)
    }
    # maximum and minimum x coordinate for the next layer
    maxX <- coordinates[coordinates$node == node,]$maxX
    minX <- coordinates[coordinates$node == node,]$minX
    
    # distance between each child nodes
    distance = (maxX - minX) / numChildren
    
    # computing the x-cooridantes of the child nodes.
    for (i in 0:(numChildren-1)) {
        coordinates[nrow(coordinates) + 1,] <- c(children[i + 1], 
                                                 maxX - (i + 0.5) * distance,
                                                 layer, 
                                                 maxX - (i * distance), 
                                                 maxX - (i + 1) * distance, 
                                                 node)
        # Compute the next layer for each child.
        coordinates = next_layer(treeDF, coordinates, children[i + 1], layer - 1)
    }
    
    return(coordinates)
}


#' find_path
#'
#' Finding the path between a node and a descendent of it.
#'
#' @param from node at which the path starts
#' @param to node at which the path ends. Must be in the subtree of \code{from}
#' @param treeDF tree in data.frame format
#'
#' @return A list of nodes from \code{from} to \code{to}
#' @noRd
find_path <- function(from, to, treeDF) {
    
    parent = treeDF[treeDF$child == to,]$parent
    
    if (from == to) {
        return(from)
    } else if (identical(parent, numeric(0))) {
        return(NULL)
    } else if (parent == from) {
        return(c(from, to))
    } else {
        path = find_path(from, parent, treeDF)
        
        if (is.null(path)) {
            return(NULL)
        } else {
            return(c(path, to))
        }
    }
}


#' Plotting transport on trees
#'
#'
#' This function visualizes the transport of mass on a tree. 
#'
#'
#' @param tree A tree structure in list format:
#' The first element is the index of the root node.
#' The other elements are vectors defining the edges of the tree. Each of these vectors has to be 
#'  of the form \code{(parent_node_index, child_node_index, edge_weight)}.
#' @param tList (optional) The mass transport as list. Each element is a vector of the form
#'  \code{(source_node_index, target_node_index, mass)}.
#' @param supply (optional) A non negative numeric vector giving the mass supply at each tree node. The value at 
#' position \eqn{i} gives the supply at node \eqn{i}.
#' @param demand (optional) A non negative numeric vector giving the mass demand at each tree node. The value at 
#' position \eqn{i} gives the supply at node \eqn{i}.
#' @param supply.col Which color to use for supply points (Default: \code{"chartreuse3"})
#' @param demand.col Which color to use for demand points (Default: \code{"dodgerblue3"})
#' @param arr.pos Where to draw the arrow head along the transport arrow (Default: \code{0.5})
#' @param arr.adj Where to center the arrow head (Default: \code{0.5})
#' @param arr.type The type of arrow head (Default: \code{"triangle"})
#' @param arr.curve The curving strength to use for the arrow (Default: \code{0.1})
#' @param arr.col The color to use for arrows (Default: \code{"red"})
#' @param arr.length The length of the arrow head, before scaling with arrow line width (Default: \code{0.1})
#' @param arr.width The width of the arrow head, before scaling with arrow line width (Default: \code{0.1})
#' @param text.draw Whether to draw text labels of nodes (Default: \code{TRUE})
#' @param text.adj One or two values defining the justification of node labels (Default: \code{NULL})
#' @param text.pos A position specifier for node labels, see \link[graphics]{text} (Default: \code{4})
#' @param pch Symbol to use for vertices, see \link[graphics]{points} (Default: \code{19}).
#' @param bg Background color for \code{pch = 21:25} (Default : \code{NULL})
#'
#' \if{html}{\figure{tree.png}}
#' \if{latex}{\figure{tree.png}{options: width=0.5in}}
#'
#'
#' @examples
#' tree <- list(1, c(1, 2, 1), c(2, 3, 1), c(3, 4, 1), c(3, 5, 1),
#'                 c(2, 6, 1), c(1, 7, 1), c(7, 8, 1), c(7, 9, 1),
#'                 c(9, 10, 1), c(9, 11, 1), c(11, 12, 1), c(11, 13, 1))
#'
#' plot_tree(tree)
#'
#'
#' supply <- c(0,0,1,2,0,0,0,0,0,0,0,0,0)
#' demand <- c(0,0,0,0,0,1,0,1,0,0,0,1,1)
#'
#' plot_tree(tree, supply = supply, demand = demand)
#'
#'
#' tList = list(c(3, 6, 1), c(4, 8, 1))
#' plot_tree(tree, tList, supply, demand)
#'
#' @export
plot_tree <- function(tree, tList = NULL , supply = NULL, demand = NULL,
                     supply.col = "chartreuse3", demand.col = "dodgerblue3",
                     arr.pos = 0.5, arr.adj = 0.5, arr.type = "triangle", arr.col = "red",
                     arr.curve = 0.1, arr.length = 0.1, arr.width = 0.1,
                     text.draw = TRUE, text.adj = NULL, text.pos = 4,
                     pch = 19, bg = NULL) {
    
    if (length(tList) == 0) {
        tList <- NULL
    }
    
    # create a tree data frame
    treeDF <- as.data.frame(do.call(rbind, tree[-1]))
    colnames(treeDF) <- c("parent", "child", "weight")
    treeDF <- treeDF[order(treeDF$parent), ]
    
    rootNode <- tree[1]
    
    # initiate a data frame to hold the coordianates for the tree nodes
    coordinates <- data.frame(c(rootNode, 0, 0, 100, -100, -1))
    colnames(coordinates) <- c("node", "x", "layer", "maxX", "minX", "parent")
    
    # compute all coordinates
    coordinates <- next_layer(treeDF, coordinates, rootNode, -1)
    
    maxLayer <- min(coordinates$layer)
    
    coordinates$layer <- -100 * coordinates$layer / maxLayer
    maxLayer <- min(coordinates$layer) # == -100
    
    
    # add the supply to the coordiantes
    if (!is.null(supply) & !is.null(demand)) {
        supDem <- supply - demand
        coordinates$supply <- supDem[coordinates$node]
    }
    
    # create an empty plot
    plot(1, type = "n", xlab = "", ylab = "", xlim = c(-110, 110), ylim = c(maxLayer-1,1), axes = FALSE)
    
    # plot the edges
    for (i in 1:nrow(treeDF)) {
        segments(coordinates[coordinates$node == treeDF[i,]$parent,]$x,
                 coordinates[coordinates$node == treeDF[i,]$parent,]$layer,
                 coordinates[coordinates$node == treeDF[i,]$child,]$x,
                 coordinates[coordinates$node == treeDF[i,]$child,]$layer)
    }
    
    # If a transport plan is given, calculate transported mass per edge
    if (!is.null(tList)) {
        treeDF$tMass <- rep(0, nrow(treeDF))
        
        for (i in 1:length(tList)) {
            # compute the paths from the trees root node to each of the two nodes.
            pathTo <- unlist(find_path(rootNode, tList[[i]][2], treeDF))
            pathFrom <- unlist(find_path(rootNode, tList[[i]][1], treeDF))
            
            # remove common prefix of paths to get path between the nodes
            while (length(pathFrom) > 1 & length(pathTo) > 1 & pathFrom[2] == pathTo[2]) {
                pathTo <- pathTo[-1]
                pathFrom <- pathFrom[-1]
            }
            
            # add up transported mass along paths(s)
            if (length(pathTo) > 1) {
                for (j in 1:(length(pathTo) - 1)) {
                    treeDF[(treeDF$parent == pathTo[j] & treeDF$child == pathTo[j + 1]), ]$tMass =
                        treeDF[(treeDF$parent == pathTo[j] & treeDF$child == pathTo[j + 1]), ]$tMass + tList[[i]][3]
                }
            }
            if (length(pathFrom) > 1) {
                for (j in 1:(length(pathFrom) - 1)) {
                    treeDF[(treeDF$parent == pathFrom[j] & treeDF$child == pathFrom[j + 1]), ]$tMass =
                        treeDF[(treeDF$parent == pathFrom[j] & treeDF$child == pathFrom[j + 1]), ]$tMass - tList[[i]][3]
                }
            }
        }
        
        arrowsDF <- treeDF[treeDF$tMass != 0, ]
        # Plotting the arrows.
        # The line width indicates the amount of mass moves along that edge.
        for (i in 1:nrow(arrowsDF)) {
            
            fromX <- coordinates[coordinates$node == arrowsDF[i,]$parent,]$x
            fromY <- coordinates[coordinates$node == arrowsDF[i,]$parent,]$layer
            toX <- coordinates[coordinates$node == arrowsDF[i,]$child,]$x
            toY <- coordinates[coordinates$node == arrowsDF[i,]$child,]$layer
            
            w <- abs(arrowsDF[i,]$tMass)
            
            if (arrowsDF[i,]$tMass > 0) {
                diagram::curvedarrow(c(fromX, fromY), c(toX, toY), 
                            arr.adj = arr.adj, arr.pos = arr.pos, arr.type = arr.type, curve = arr.curve,
                            lwd = w * arr.width / 0.1, lcol = arr.col, arr.col = arr.col, 
                            arr.length = arr.length * w, arr.width = arr.width * w)
            } else {
                diagram::curvedarrow(c(toX, toY), c(fromX, fromY), 
                            arr.adj = arr.adj, arr.pos = arr.pos, arr.type = arr.type, curve = arr.curve,
                            lwd = w * arr.width / 0.1, lcol = arr.col, arr.col = arr.col,
                            arr.length = arr.length * w, arr.width = arr.width * w)
            }
        }
    }
    
    
    if (is.null(supply) | is.null(demand)) {
        # If the supply and demand are not given, plot all nodes in black.
        points(coordinates$x, coordinates$layer, pch = pch, bg = bg)
    } else {
        # Otherwise plot supply and demand nodes in their respective colors
        # The radius of the node indicates the amount of supply / demand: The bigger the node
        # the more mass is supplied or demanded.
        # Nodes without supply or demand are plotted in black.
        points(coordinates[coordinates$supply == 0, ]$x, coordinates[coordinates$supply == 0, ]$layer,  
               pch = pch, bg = bg)
        
        points(coordinates[coordinates$supply > 0, ]$x, coordinates[coordinates$supply > 0, ]$layer,
               pch = pch, bg = bg, cex = abs(coordinates[coordinates$supply > 0, ]$supply),  col = supply.col)
        
        points(coordinates[coordinates$supply < 0, ]$x, coordinates[coordinates$supply < 0, ]$layer,
               pch = pch, bg = bg, cex = abs(coordinates[coordinates$supply  < 0, ]$supply),  col = demand.col)
    }
    
    # Adding the keys to the plot.
    if (text.draw)
        text(coordinates$x, coordinates$layer, labels = coordinates$node, adj = text.adj, pos = text.pos)
    
}


