#' Title: Generates a dissimilarity matrix from weighted Euclidean matrix
#'
#' @param marker either a marker name from the rownames of X or an index
#' @param x matrix of the expression of the cells from a patient sample
#' @param y matrix of the expression of the cells taken from all the healthy samples.
#'          x and y should have the same dimension and ordering, as in WeightedDissimilarityMatrix
#' @return weight for specific gene marker calculated by how much it is conserved in the trajectory between x and y
#'
#' @import proxy
#' @import cellAlign
#'
#' @export
#' @examples
weightCalc <- function (marker, x, y, const, class = 'nonSig'){
  dist.method = "Euclidean"
  lm = proxy::dist(x[marker,],y[marker,],method=dist.method)
  cost = cellAlign:::globalCostMatrix(lm)
  traj_dist = cost$costMatrix[nrow(cost$costMatrix),ncol(cost$costMatrix)]
  if(!(class %in% c('Sig','nonSig'))){
    message('unrecognized class')
    stop()
  }
  if(class == 'nonSig'){w = 1/(const + traj_dist)}
  if(class == 'Sig'){w = 1/(1+exp(traj_dist - const))}
  return(w)
}
