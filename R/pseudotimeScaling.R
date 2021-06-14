#' Title Generates the tuMap pseudotime axis based on tuMap alignment
#'
#' @param ptCancer original pseudotime values for each cell in the cancer sample
#' @param interScaledCancer interpolated scaled data of the cancer sample
#' @param interScaledHealthy interpolated scaled data of the averaged healthy trajectory
#' @param gAlign tuMap's global alignment using the selected set of markers
#'
#' @return scaled pseudotime values
#' @export
#' @import FNN
#'
pseudotimeScaling <- function(ptCancer, interScaledCancer, interScaledHealthy, gAlign){
  realTraj = ptCancer
  intTraj = interScaledCancer$traj
  intAssign = FNN::get.knnx(intTraj, realTraj, k = 1)$nn.index
  gAlign4Map = gAlign$align[[1]]

  #convert each intAssign to the aligned interpolated point
  hTraj = sapply(intAssign, function(iAs){return(interScaledHealthy$traj[median(gAlign4Map$index2[gAlign4Map$index1 == iAs])])})
  return(hTraj)
}
