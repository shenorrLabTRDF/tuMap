require(ggplot2)
require(pheatmap)
require(stringr)
require(Biobase)
require(reshape2)
require(cellAlign)
require(glmnet)
require(survival)
require(survAUC)
library(survminer)
require(tumap)
require(psupertime)
require(caret)
require(emdist)
require(flowCore)
require(SingleCellExperiment)

#' Title Generates an expression set from a single or multiple FCS files
#'
#' @param fileNames a vector with the files' names as strings
#' @param fileDir a string, full path of the directory where files are located
#' @param transform logical, whether or not to perform arcsinh transformation
#'
#' @return an expression set with the expression data
#' @export
#'
#' @examples readFCSeset(fileName = 'Healthy2_Basal', fileDir = file.path(dataDir,'fcsFilesHealthy'), transform = T)

readFCSeset <- function(fileNames, fileDir, transform = T){

  #make a large expression set with expression values across all files:
  exprsAllFilesList = lapply(fileNames, function(fileName){
    print(fileName)
    relFiles = list.files(fileDir)[grep(fileName, list.files(fileDir))]

    fcsFiles = read.flowSet(files = file.path(fileDir, relFiles), transformation = FALSE)

    #read expression data of all the FCS files in the folder:
    exprsFiles = do.call('cbind', lapply(1:length(fcsFiles), function(i){
      popName = str_replace(str_replace(str_extract(relFiles[i], '\\_[a-zA-Z0-9]+\\.fcs'),'\\.fcs',''),'\\_','')
      expData = exprs(fcsFiles[[i]])

      #convert rownames into epitopes:
      colnames(expData) = sapply(colnames(expData), function(channel){
        return(pData(parameters(fcsFiles[[1]]))[pData(parameters(fcsFiles[[1]]))[,'name'] == channel,'desc'])
      })

      expData = cbind(expData, popName = popName, fileName = fileName)
      print(dim(expData))
      return(t(expData))
    }))

    return(exprsFiles)
  })

  #merge the expression matrices using the common set of markers:
  commonMarkers = Reduce(intersect, lapply(exprsAllFilesList, function(x){return(rownames(x))}))

  #merge the expression data matrices using the common set of markers:
  exprsAllFiles = do.call('cbind', lapply(exprsAllFilesList, function(x){return(x[commonMarkers,])}))

  #generate phenotypic data:
  pData = data.frame(popName = exprsAllFiles['popName',], fileName = exprsAllFiles['fileName',])
  phenoData <- new("AnnotatedDataFrame", data=pData)

  #convert expression matrix to numeric:
  exprsFiles = exprsAllFiles[!(rownames(exprsAllFiles) %in% c('popName','fileName')),]
  class(exprsFiles) <- "numeric"

  #remove duplicated cells:
  dupCells = which(duplicated(t(exprsFiles)))
  if(length(dupCells) > 0){
    exprsFiles = exprsFiles[,-dupCells]
    phenoData = phenoData[-dupCells,]
  }

  #organize data as an expression set
  fcsEset = new("ExpressionSet", exprs = exprsFiles,
                phenoData = phenoData)

  if(!transform){return(fcsEset)}

  #arcsinh transformation:
  markersTransform = rownames(fcsEset)[!(rownames(fcsEset) %in% c('Time','Cell_length'))]
  transMat = do.call('rbind', lapply(markersTransform, function(marker){
    print(marker)
    return(asinh(exprs(fcsEset)[marker,]/5))
  }))
  exprs(fcsEset)[markersTransform,] = transMat

  return(fcsEset)
}

#' Title Developmental classification for AML CyTOF dataset
#'
#' @param transEsetBQuery an expressioset of CyTOF data of the sample (markers along rows and cells along columns)
#' @param refMarkers a vector of markers used for classification
#' @param popMeans Mean expression values of the markers used for classification per cell population
#' @param popSx Covariance matrix of the markers used for classification per cell population
#'
#' @return a vector of the predicted cell type per cell
#' @export
#'
#' @examples
devMapperAML <- function(transEsetBQuery, refMarkers, popMeans, popSx){
  #calculate Mahanabolis distance between each cell and the reference populations:
  mahDistPerCell = do.call('cbind',lapply(names(popMeans), function(pop){
    mahDist <- sqrt(mahalanobis(t(exprs(transEsetBQuery))[,refMarkers],
                                popMeans[[pop]][refMarkers],
                                popSx[[pop]][refMarkers,refMarkers], inverted = TRUE))
    return(mahDist)
  }))
  colnames(mahDistPerCell) = names(popMeans)
  #per cell, classify according to the minimal distance:
  predictedCellType = names(popMeans)[apply(mahDistPerCell,1,which.min)]
  predictedCellType[apply(mahDistPerCell,1,min) > length(popMeans[[1]])] = 'unclassified'
  return(predictedCellType)
}

#' Title plots the pseudotime of sequential developmental stages in the sample
#'
#' @param eset an expressionset with pseudotime information in its phenotypic data
#' @param method pseudotime method used for plot (name of the pseudotime value in the phenotypic data)
#' @param sampleName name of the sample asshould appear in the plot
#'
#' @return
#' @export
#'
#' @examples plotDev(eset, method = 'SCORP', sampleName = 'healthy)
plotDev <- function(eset, method, sampleName = '', orderedPops){
  df = pData(eset)
  method = method
  df$popName = factor(df$popName, levels = orderedPops)
  ggplot(df, aes_string(x = 'popName', y = method), environment = environment()) + geom_boxplot() + coord_flip() +
    ggtitle(paste(sampleName, 'cell populations by pseudotime')) + theme_classic()
}

#' Title Interpolates markers' expression data along developmental trajectories
#'
#' @param transEsetB an expressioset with markers' expression levels by single cells
#' @param markers a vector of markers whose expression will be interpolated
#' @param method the trajectory inference methodology along which the markers expression levels will be interpolated
#'
#' @return
#' @export
#'
#' @examples
cellAlignInter <- function(transEsetB, markers, method = 'dpt', scale = T){
  #interpolate expression:
  interData = interWeights(expDataBatch = exprs(transEsetB)[markers,],
                           trajCond = pData(transEsetB)[,method],
                           winSz = 0.05, numPts = 200)
  if(!scale){return(interData)}
  interDataScaled = scaleInterpolate(interData)
  return(interDataScaled)
}


#' Title Generates pseudotime trajectory based on cellular annotation into cell populations
#'
#' @param esetALLDev an expressionset of the cancer expression data
#' @param markers Developmental markers based on which the pseudotime trajectory is generated
#' @param orderedPops a vector of the ordered cell populations
#'
#' @return an expressionset with an additional pseudotime value per cell
#' @export
#' @import psupertime
#'
#' @examples
BuildDevTrajHealthyPops <- function(esetALLDev, markers,
                                    orderedPops = c('HSC','Progenitor1','Progenitor2','Progenitor3','PreProB','ProB1','ProBII','PreBI','PreBII','ImmatureBI','immatureBII','B')){

  #if the number of cells in the SCE exceeds some threshold, sample only threshold (in order to save running time):
  message('subsampling of cells')
  szThresh = 10000
  if(ncol(esetALLDev) > szThresh){esetALLDevSub = esetALLDev[,sample(x = 1:ncol(esetALLDev), size = szThresh)]
  }else{esetALLDevSub = esetALLDev}

  pretend.cell.labels <- esetALLDevSub$popName
  counts <- exprs(esetALLDevSub)

  sce <- SingleCellExperiment(list(logcounts=counts),
                              colData=DataFrame(label=pretend.cell.labels))

  #determine the order of the populations based on the information given by the user:
  ordering = sapply(esetALLDevSub$popName, function(pop){
    return(which(orderedPops == pop))
  })

  ordering = as.numeric(ordering)

  #apply psupertime algorithm:
  message('calculate psupertime values')
  psuper_obj  = psupertime(x = sce[markers,], y = ordered(ordering, levels = 1:max(ordering)),
                           sel_genes = 'all', smooth = F, penalization = 'best', score = 'xentropy')

  #if sampling was performed, project the other cells onto the same space using the set of model's coefficients:
  if(ncol(esetALLDev) > szThresh){
    scaledExp = apply(t(exprs(esetALLDev)[as.character(psuper_obj$beta_dt$symbol),]), 2, scale)
    esetALLDev$psuper = scaledExp %*% matrix(as.numeric(psuper_obj$beta_dt$beta), ncol=1)
  }else{esetALLDev$psuper = psuper_obj$proj_dt$psuper}

  #Trajectory trimming:
  message('trajectory trimming & scaling')
  threshDens = 0.005
  densPsupertime = density(esetALLDev$psuper)
  minTrim = densPsupertime$x[min(which(densPsupertime$y > threshDens))]
  maxTrim = densPsupertime$x[max(which(densPsupertime$y > threshDens))]

  esetALLDev = esetALLDev[,(esetALLDev$psuper > minTrim) & (esetALLDev$psuper < maxTrim)]
  esetALLDev$psuperScaled = (esetALLDev$psuper - min(esetALLDev$psuper))/(max(esetALLDev$psuper) - min(esetALLDev$psuper))

  return(esetALLDev)
}
