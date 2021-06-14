## The healthy trajectory:
fcsFilesHDir = '~/MDPhD/AML/revisionAnalysis2/AML/Data/fcsFilesH6Pops'
fileNames = c('036','049','050','067','092','132','138','247','277')
esetH = readFCSeset(fileNames, fcsFilesHDir, transform = T)

devPops = c('monoblasts','promonocytes','monocytes')
devMarkers = c('CD34','CD117','CD33','CD64','CD13','CD14','CD11b')
esetHDev = esetH[devMarkers,esetH$popName %in% devPops]

#sample cells with inverse relation to the population's frequency:
sampProb = table(esetHDev$popName)/ncol(esetHDev)
esetHDevSamp = esetHDev[,sample(1:ncol(esetHDev), size = 10000, prob = sapply(esetHDev$popName, function(pop){return(1/sampProb[pop])}))]

#apply PCA:
pcaRes = prcomp(t(exprs(esetHDevSamp)), center = T)$x[,c('PC1','PC2')]
pcaResDf = data.frame(pcaRes, popName = esetHDevSamp$popName)
pcaResDf$popName = factor(pcaResDf$popName, levels = c('monoblasts','promonocytes','monocytes'))
ggplot(pcaResDf, aes(x = PC1, y = PC2, color = popName)) + geom_point() + theme_classic() +
  scale_color_manual(values = wes_palette(n = 3, name = "GrandBudapest1")) +
  ggtitle('trajectory of sampled data with inverse relation to population frequency')

esetHDev = BuildDevTrajHealthyPops(esetHDev, markers = devMarkers, orderedPops = c('monoblasts','promonocytes', 'monocytes'))
#Plot the developmental populations ordering along the trajectory:
plotDev(esetHDev, method = 'psuperScaled', orderedPops = c('monoblasts','promonocytes', 'monocytes'))

#markers expression dynamics:
interScaledH = cellAlignInter(esetHDev, markers = devMarkers, method = 'psuperScaled')
pheatmap(interScaledH$scaledData, cluster_rows = F, cluster_cols = F,
         main = 'markers expression dynamics along the healthy trajectory',
         color = colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))

#density along the trajectory:
ptXrossInd = data.frame(ind = esetHDev$fileName, pt = esetHDev$psuperScaled)
ggplot() + geom_density(ptXrossInd, mapping = aes(x = pt, y = ..density.., group = ind), color = 'grey') +
  geom_density(ptXrossInd, mapping = aes(x = pt, y = ..density..)) +
  theme_classic() +
  ggtitle('cellular density along the healthy trajectory across patients')

classMarkers = c('CD34','CD123','CD33','CD64','CD11c','CD14','CD13','CD117','CD71','CD11b')
#define popMeans-  list of elements, in each the mean expression value of specific markers:
popMeans = lapply(unique(esetH$popName), function(pop){
  return(apply(exprs(esetH)[classMarkers,esetH$popName == pop],1,mean))
})
names(popMeans) = unique(esetH$popName)
#calculate the markers' covariance matrix per population:
popSx = lapply(unique(esetH$popName), function(pop){
  covRes = cov(t(exprs(esetH)[classMarkers,esetH$popName == pop]))
  return(covRes)
})
names(popSx) = unique(esetH$popName)

esetH = esetH[classMarkers,]

#classify the cells based on the data:
esetH$predCellTypes = devMapperAML(transEsetBQuery = esetH, refMarkers = classMarkers, popMeans, popSx)
confMatRes = caret::confusionMatrix(data = factor(esetH$predCellTypes, levels = unique(esetH$popName)),
                                    reference = factor(esetH$popName, levels = unique(esetH$popName)))
confMatRes$overall['Accuracy']
confMatRes$table
NormedConfMat = apply(confMatRes$table,2,function(x){return(x/sum(x))})
pheatmap(NormedConfMat, cluster_cols = F, cluster_rows = F,
         main = 'confusion matrix\nreference labels on columns and prediction on rows')

## Developmental trajectory assembly for AML samples
#define markers required for classification and trajectory building, and the developmental cell populations:
classMarkers = c('CD34','CD123','CD33','CD64','CD11c','CD14','CD13','CD117','CD71','CD11b')
devMarkers = c('CD34','CD117','CD33','CD64','CD13','CD14','CD11b')
devPops = c('monoblasts','promonocytes','monocytes')

#Classify each cell to the closest population:
esetAML = readFCSeset(fileNames = sample, fileDir = AMLDir, transform = T)
esetAML$popName = devMapperAML(esetAML, classMarkers, popMeans, popSx)

#Generate a trajectory using the developmental cell populations:
esetAMLDev = esetAML[,esetAML$popName %in% devPops]
esetAMLDev = BuildDevTrajHealthyPops(esetAMLDev, devMarkers, orderedPops = devPops)

##Application of the tuMap algorithm on the AML samples
AMLwTrajDir = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/AMLSamplesWTraj'
AMLwTrajFiles = list.files(AMLwTrajDir)
tuMapTrajDir = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/AMLsamplesDevTrajwTumap'

markersWeightsXrossPatients = do.call('rbind', lapply(AMLwTrajFiles, function(file){
  print(file)
  esetAMLDev = readRDS(file = file.path(AMLwTrajDir,file))
  sampleID = str_extract(file, '[0-9]+\\_[0-9]+')

  # calculate the weighted dissimilarity matrix between our trajectory and the healthy:
  interScaledAML = cellAlignInter(esetAMLDev, devMarkers, method = 'psuperScaled')

  wDisMat = tumap::WeightedDissimilarityMatrix(x = interScaledAML$scaledData[devMarkers,],
                                               y = interScaledH$scaledData[devMarkers,], const = 1)

  #global alignment with the weighted distance matrix:
  gAlignWweights = cellAlign::globalAlign(wDisMat$dissimilarity_matrix, step.pattern = symmetric1)
  esetAMLDev$tuMapPt = tumap::pseudotimeScaling(ptCancer = esetAMLDev$psuperScaled, interScaledCancer = interScaledAML,
                                                interScaledHealthy = interScaledH, gAlign = gAlignWweights)

  saveRDS(esetAMLDev, file = file.path(tuMapTrajDir, paste0(sampleID, 'wtuMap.rds')))

  return(data.frame(sampleID = sampleID, t(wDisMat$weights)))
}))
rownames(markersWeightsXrossPatients) = markersWeightsXrossPatients$sampleID

#show the markers' weights as violin plots:
markersWeightsXrossPatients_melt = melt(markersWeightsXrossPatients, id.vars = 'sampleID')
meanWeights = sapply(devMarkers, function(marker){return(mean(markersWeightsXrossPatients_melt$value[markersWeightsXrossPatients_melt$variable == marker]))})
ordMarkers = names(meanWeights)[order(meanWeights)]
markersWeightsXrossPatients_melt$variable = factor(markersWeightsXrossPatients_melt$variable, levels = ordMarkers)
ggplot(markersWeightsXrossPatients_melt, aes(x = variable, y = value, fill = variable)) + geom_jitter(width = 0.1, color = 'black', size = 0.5) + geom_violin(alpha = 0.2) +
  theme_classic() + ggtitle('markers weights across samples')

#For tuMap alignment error estimation:
hInds = unique(esetHDev$fileName)
intDataHInd = lapply(hInds, function(hInd){
  esetHInd = esetHDev[devMarkers,esetHDev$fileName == hInd]
  esetHInd$psuperOrg = esetHInd$psuperScaled
  esetHInd = BuildDevTrajHealthyPops(esetHInd, devMarkers, orderedPops = c('monoblasts','promonocytes', 'monocytes'))
  return(list(trajInd = esetHInd$psuperScaled,
              trajAvH = esetHInd$psuperOrg,
              interScaled = cellAlignInter(esetHInd, devMarkers, method = 'psuperScaled', scale = T)))
})
saveRDS(intDataHInd, file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/savedDataHealthy/intDataHInd.rds')

#read the expression matrix for one AML sample:
esetAMLDev = readRDS(file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/AMLSamplesWTraj/1_1421wTraj.rds')

# calculate the weighted dissimilarity matrix between our trajectory and the healthy:
interScaledAML = cellAlignInter(esetAMLDev, devMarkers, method = 'psuperScaled')

wDisMat = tumap::WeightedDissimilarityMatrix(x = interScaledAML$scaledData[devMarkers,],
                                             y = interScaledH$scaledData[devMarkers,], const = 1)

tuMapError = tumap::estimatetuMapErr(wDisMat, intDataHInd, interScaledH)
plot(tuMapError)

##Cellular density:
samplesData = read.csv(file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/clinicalData/samplesData.csv', stringsAsFactors = F)
tuMapTrajDir = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/AMLsamplesDevTrajwTumap'
tuMapTrajFiles = list.files(tuMapTrajDir)

#calculate the cellular density along the tuMap trajectory across the AML samples:
tumapDensXrossSamples = do.call('rbind', lapply(tuMapTrajFiles, function(tuMapFile){
  esetAML = readRDS(file = file.path(tuMapTrajDir, tuMapFile))
  tumapDens = density(esetAML$tuMapPt, from = 0, to = 1, n = 128, adjust = 2)
  return(tumapDens$y)
}))
rownames(tumapDensXrossSamples) = str_extract(tuMapTrajFiles,'[0-9]+\\_[0-9]+')

#calculate the healthy cellular density along the pseudotime axis:
hDens = density(esetHDev$psuperScaled, from = 0, to = 1, n = 128, adjust = 2)

#calculate per sample the EMD distance of cellular density from the healthy:
EMDtuMapH = sapply(tuMapTrajFiles, function(tuMapFile){
  esetAML = readRDS(file = file.path(tuMapTrajDir, tuMapFile))
  tumapDens = density(esetAML$tuMapPt, from = 0, to = 1, n = 128, adjust = 2)
  return(emdist::emd(cbind(tumapDens$y, tumapDens$x),
                     cbind(hDens$y, hDens$x)))
})
names(EMDtuMapH) = str_extract(names(EMDtuMapH),'[0-9]+\\_[0-9]+')

#show the dynamics of EMDH of patient2 over time:
EMDHDynamicsPt2 = do.call('rbind', lapply(samplesData$SampleID[samplesData$Patient == 2], function(sampleID){
  return(data.frame(patient = 2, sampleID = sampleID,
                    timePoint = samplesData$Timepoint[samplesData$SampleID == sampleID],
                    EMDH = EMDtuMapH[sampleID]))
}))
EMDHDynamicsPt2$timePoint = factor(EMDHDynamicsPt2$timePoint, levels = c('d0','d14','rem','rel'))
ggplot(EMDHDynamicsPt2, aes(x = timePoint, y = EMDH, group = patient)) + geom_line(linetype = 'dashed') + geom_point() +
  theme_classic() + ggtitle('EMDH of longitudinal samples of patient 2')

#order the individuals by the EMD distance from the healthy:
ordInd = names(EMDtuMapH)[order(-EMDtuMapH)]
tumapDensXrossSamplesOrd = tumapDensXrossSamples[ordInd,]

#add annotations per sample describing the time point in which the sample was taken:
annotRows = do.call('rbind', lapply(rownames(tumapDensXrossSamplesOrd), function(sampleID){
  return(data.frame(patient = samplesData$Patient[samplesData$SampleID == sampleID],
                    timePoint = samplesData$Timepoint[samplesData$SampleID == sampleID],
                    EMDH = EMDtuMapH[sampleID]))
}))
rownames(annotRows) = rownames(tumapDensXrossSamplesOrd)

#add the healthy density at the bottom of the table:
hDens = readRDS(file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/savedDataHealthy/densityTrajH.rds')
tumapDensXrossSamplesOrd = rbind(tumapDensXrossSamplesOrd, hDens$y)
rownames(tumapDensXrossSamplesOrd)[nrow(tumapDensXrossSamplesOrd)] = 'H'
annotRows = rbind(annotRows, data.frame(patient = 'H', timePoint = 'H', EMDH = 0))
rownames(annotRows)[nrow(annotRows)] = 'H'
pheatmap(tumapDensXrossSamplesOrd, cluster_cols = F, cluster_rows = F, border_color = NA,
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), show_rownames = F, gaps_row = nrow(tumapDensXrossSamplesOrd) - 1,
         main = 'cellular density along the tuMap pseudotime axis', annotation_row = annotRows[,c('timePoint', 'EMDH'), drop = F])

#paired comparison between d0 and d14 healthy EMD:
EMDHD0D14 = dcast(annotRows, formula = patient ~ timePoint, value.var = 'EMDH')
EMDHD0D14 = EMDHD0D14[!is.na(EMDHD0D14$d14),c('d0','d14')]
t.test(EMDHD0D14$d0, EMDHD0D14$d14, paired = T)

#paired comparison between relapse and d14 healthy EMD:
EMDHD0Rel = dcast(annotRows, formula = patient ~ timePoint, value.var = 'EMDH')
EMDHD0Rel = EMDHD0Rel[!is.na(EMDHD0Rel$d14) & !is.na(EMDHD0Rel$rel), c('d14','rel')]
t.test(EMDHD0Rel$rel, EMDHD0Rel$d14, paired = T)

#boxplot of all the data:
ggplot(subset(annotRows, timePoint %in% c('d0','d14','rel')), aes(x = timePoint, y = EMDH)) + geom_boxplot() + theme_classic() +
  ggtitle('EMD from the healthy distribution along tuMap trajectory')
t.test(annotRows$EMDH[annotRows$timePoint == 'd0'], annotRows$EMDH[annotRows$timePoint == 'd14'])

##Correlation with outcome:
samplesData = read.csv(file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/clinicalData/samplesData.csv', stringsAsFactors = F)
patientsData = read.csv(file = '~/MDPhD/AML/revisionAnalysis2/AML/savedData/clinicalData/patientsData.csv', stringsAsFactors = F)

#the Surv() function requires the data to be organized as: 0-alive, 1 - dead (https://www.rdocumentation.org/packages/survival/versions/2.11-4/topics/Surv)
patientsData$status.alive = 1 - patientsData$is.alive

#get only those patients with d14 samples:
samplesData = samplesData[samplesData$SampleID %in% rownames(tumapDensXrossSamples),]
d14Patients = samplesData$Patient[samplesData$Timepoint == 'd14']
patientsDataD14 = patientsData[d14Patients,]

#organize predictors data:
densityD0D14 = do.call('rbind', lapply(d14Patients, function(pt){
  #diagnosis age:
  diagAge = patientsData$age.diagnosis[patientsData$PatientID == pt]
  #sample ID of the sample taken at day 0 and its cellular density along the tuMap axis:
  d0Sample = samplesData$SampleID[samplesData$Patient == pt & samplesData$Timepoint == 'd0']
  d0Density = tumapDensXrossSamples[d0Sample,]
  #sample ID of the sample taken at day 14 and its cellular density along the tuMap axis:
  d14Sample = samplesData$SampleID[samplesData$Patient == pt & samplesData$Timepoint == 'd14']
  d14Density = tumapDensXrossSamples[d14Sample,]
  #concatenate the data as a vector:
  predPt = c(diagAge, d0Density, d14Density)
  names(predPt) = c('diagAge', paste0('d0_',1:length(d0Density)), paste0('d14_',1:length(d14Density)))
  return(predPt)
}))
rownames(densityD0D14) = paste0('pt',d14Patients)

#generate a glmnet regularized COX model:
set.seed(1)
cv.fitTraj <- cv.glmnet(x = densityD0D14,
                        Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']), family="cox",
                        alpha = 1, nfolds = 5)

#choose the best lambda:
best.lambdaTraj = cv.fitTraj$lambda.min

#get the coefficients of the model with the best lambda:
CoefficientsTraj <- coef(cv.fitTraj, s = best.lambdaTraj)
Active.IndexTraj <- which(CoefficientsTraj != 0)
Active.CoefficientsTraj <- CoefficientsTraj[Active.IndexTraj]
coefBestLambdaTraj = data.frame(fet = colnames(densityD0D14)[Active.IndexTraj],
                                coef = CoefficientsTraj[Active.IndexTraj])

#fit a full model with all the data
lasso.modTraj <- glmnet(x = densityD0D14,
                        Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']), family="cox",
                        alpha = 1, lambda = best.lambdaTraj)

#generate a Kaplan Meier survival curve for the best lambda:
lpTraj = predict(lasso.modTraj, newx = densityD0D14, s = best.lambdaTraj)

#put a threshold on the median lpAll:
threshTraj = median(lpTraj[,1])

#calculate log-rank test for 2 groups stratified by the median risk:
patientsDataD14$groupsModelTraj = (lpTraj[,1] <= threshTraj)

# plot the Kaplam Meier plot (Figure 4d)
survDiffObj = survdiff(Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']) ~ patientsDataD14$groupsModelTraj)
pv = 1 - pchisq(survDiffObj$chisq, length(survDiffObj$n) - 1)
fit = survfit(Surv(patientsDataD14[,'followup.alive'], patientsDataD14[,'status.alive']) ~ patientsDataD14$groupsModelTraj)
ggsurvplot(fit, data = patientsDataD14, pval = TRUE)

#show the distribution of d14 cellular density along the trajectory for the low vs. high risk - Supp. Fig. 9d:
densityD14 = densityD0D14[,str_detect(colnames(densityD0D14),'d14')]
densityD14_m = melt(densityD14)
densityD14_m$ind = as.numeric(str_replace(densityD14_m$Var2,'d14_',''))
densityD14_m$traj = (densityD14_m$ind - 1)*(1/127)
densityD14_m$riskGrp = sapply(densityD14_m$Var1, function(ptid){
  return(patientsDataD14$groupsModelTraj[patientsDataD14$PatientID == str_extract(ptid,'[0-9]+')])})
ggplot(densityD14_m, aes(x = traj, y = value, group = Var1, color = riskGrp)) + geom_line() + theme_classic() +
  ggtitle('day 14 density distributions for the high vs. low risk patients')
