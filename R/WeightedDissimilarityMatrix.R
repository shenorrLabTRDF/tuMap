#------------------------------------------------------------------------
# Function WeightedDissimilarityMatrix: generate dissimilarity matrix from weighted Euclidean matrix
#------------------------------------------------------------------------
#' Title
#'
#' @param x matrix of the expression of the cells from a patient sample: markers along rows and cells along columns
#' @param y matrix of the expression of the cells taken from all the healthy samples: markers along rows and cells along columns
#'          x and y should have the same dimension (number of rows)
#'
#' @return dissimilarity_struct -  a list with 2 fields:
#'         dissimilarity_matrix - the dissimilarity matrix between x and y computed using weighted Euclidean distance
#'         weights - the weights that were used for the distance calculation (per gene)
#' @export
#' @import proxy
#'
#' @examples
WeightedDissimilarityMatrix <- function(x,y, const = 1, class = 'nonSig'){
  markers = 1:nrow(x)
  # check the cost for each of the markers individually
  weights_dist = unlist(lapply(markers,function(marker){return(weightCalc(marker, x, y, const, class))}))
  # normalize the weights to act like a probability function
  weights_dist = weights_dist/sum(weights_dist)
  # calculate a dissimilarity matrix using the weights:
  dissimilarity_matrix = proxy::dist(t(x*sqrt(weights_dist)), t(y*sqrt(weights_dist)))

  # output the parameters
  dissimilarity_struct = list()
  dissimilarity_struct[["dissimilarity_matrix"]] = dissimilarity_matrix
  dissimilarity_struct[["weights"]] = weights_dist
  names(dissimilarity_struct[["weights"]]) = rownames(x)
  dissimilarity_struct
}
