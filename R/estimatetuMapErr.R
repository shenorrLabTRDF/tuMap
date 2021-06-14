#' Title Estimates alignment error resulting from tuMap alignment
#'
#' @param wDisMatPt the result of the patient's alignment with the healthy that includes markers' weights
#' @param hExp a list of interpolated expression matrices, each corresponds to a healthy individual with the relevant markers
#' Each element in hExp should include the following attributes: trajInd: the individual trajectory; trajAvH: the averaged healthy trajectory; interScaled: the interpolated, scaled expression
#' @param interScaledH interpolated scaled expression of markers along the averaged healthy
#'
#' @return a vector of errors per location along the trajectory
#' @export
#' @import cellAlign
#' @import proxy
#'
#' @examples
estimatetuMapErr <- function(wDisMatPt, hExp, interScaledH){
  #run on healthy individuals
  errorAlongTraj = do.call('rbind', lapply(hExp, function(hIntExpInd){
    #apply the markers' weighting scheme to calculate the weighted dissimilarity matrix:
    dissimilarity_matrix = proxy::dist(t(hIntExpInd$interScaled$scaledData*sqrt(wDisMatPt$weights)),
                                       t(interScaledH$scaledData*sqrt(wDisMatPt$weights)))
    #calculate the alignment using the weighted dissimilarity matrix:
    gAlignWweights = cellAlign::globalAlign(dissimilarity_matrix, step.pattern = symmetric1)
    #calculate tuMap pseudotime:
    tuMapPt = pseudotimeScaling(ptCancer = c(hIntExpInd$trajInd), interScaledCancer = hIntExpInd$interScaled,
                                interScaledHealthy = interScaledH, gAlign = gAlignWweights)
    #calculate mean error & interpolated error along the trajectory
    interError = interWeights(expDataBatch = rbind(c(abs(tuMapPt - hIntExpInd$trajAvH)),
                                                   c(abs(tuMapPt - hIntExpInd$trajAvH))),
                              trajCond = hIntExpInd$trajAvH,
                              winSz = 0.05, numPts = 200)
    return(interError$interpolatedVals[1,])
  }))
  return(apply(errorAlongTraj,2,median))
}
