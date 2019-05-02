# tuMap
Alignment of tumors' trajectories

## Introduction
Abnormal differentiation is a key feature of cancer, yet currently there is no framework that enables a comparative analysis of differentiation processes across patients while preserving their individual-level resolution. For this purpose, we developed devMap, an algorithm that uses high-dimensional trajectory alignment to anchor cancer-related developmental processes to a common backbone process, thus allowing for their systematic comparison. 
Here we provide a tutorial exemplifying how to use devMap R package to align an AML trajectory to the healthy myeloid differentiation process.

## 1. devMap installation
You can install the devMap R package using the package archieve file.

Now you can load the package functionalities using the following code:
```
library(devMap)
```

devMap package depends on dtw and cellALign package. Make sure you have these packages installed!
```
library(dtw)
library(cellAlign)
```

## 2. Upload data
The data used in this tutorial is attached to the package and is available upon installation.
There are some kinds of data used in this tutorial:

#### Raw CyTOF expression data of healthy myeloid cells:
```
data(expSetH)
```
To visualize the healthy trajectory, one can use the following code:

```
require(ggplot2)
ggplot(expSetH, aes(x = DC1, y = DC2, color = dpt)) + geom_point()
```

To visualize markers expression dynamics along the healthy trajectory, you can use the following code:
```
require(pheatmap)
ptBins = seq(0,1,by = 0.05)
markers4plot = c('CD34','CD117','CD38','HLADR','CD123','CD45RA','CD64','CD14','CD11b','CD33','CD45','CD13')

#per bin, calculate the median expression of each marker across cells assigned to this bin:
expPerBin = do.call('rbind', lapply(2:length(ptBins), function(i){
  expDynMarkers = expSetH[(expSetH$dpt < ptBins[i]) & (expSetH$dpt >= ptBins[i-1]), markers4plot]
  return(apply(expDynMarkers,2,median))
}))
pheatmap(t(expPerBin)[markers4plot,], 
         cluster_cols = F, cluster_rows = F, colorRampPalette(rev(brewer.pal(n = 7, name = "Spectral")))(100), border_color = NA,
         main = 'markers expression dynamics along pseudotime')

```

#### Interpolated scaled expression of the averaged healthy trajectory:
```
data(interDataScaledHAveraged)
```

#### Interpolated scaled expression of each healthy individual:
```
data(interDataScaledHInd)
```

#### Raw CyTOF expression data of an AML patient:
```
data(rawDataAML)
```

#### Interpolated scaled expression data of an AML patient:
```
data(interDataScaledAML)
```

## 3. Markers selection step
In this step, devMap identifies the set of markers exhibiting conserved expression dynamics along the trajectory as compared to the healthy backbone by a step-wise exclusion of non-conserved markers to optimize the alignment of these two trajectories.

```
markers = c('CD38','CD34','CD123','CD33','CD45RA','CD64','CD14','CD13','CD117','CD11b','CD45')
fetSelRes = markersOrdering(markers, 
                            interScaledCancer = interDataScaledAML, 
                            interScaledHealthy = interDataScaledHAveraged, 
                            interScaledHPerInd = interDataScaledHInd)
```

Plot the resulting exclusion ordering. The x-axis denotes the markers exclusion ordering. "dist" parameter denotes the average alignment cost resulting from the alignment of AML to the averaged healthy trajectory whereas "refCost" refers to the average alignment cost resulting from alignment of 2 healthy individuals. devMap chooses the largest set of markers that yields a better alignment cost as compared to alignment of healthy trajectories.
```
require(reshape2)
fetSelResm = melt(fetSelRes, id.vars = 'markers')
ggplot(fetSelResm, aes(x = markers, y = value, group = variable, color = variable)) + 
  geom_point() + geom_line() + theme_classic() + ggtitle('markers selection ordering')
```

Now we can extract the set of alignment markers chosen for this AML sample:
```
alignMarkers = markersSelection(fetSelRes)
```

## 4. Alignment of the AML and the averaged healthy trajectories using the selected set of conserved markers
```
gAlign = cellAlign::globalAlign(interDataScaledAML$scaledData[alignMarkers,], 
                                interDataScaledHAveraged$scaledData[alignMarkers,], step.pattern = symmetric1)
plotAlign(gAlign)
```

## 5. Scaling of the pseudotime axis of the AML sample based on the alignment
Follwing this step, the pseudotime axis of the cancer samlpe is alignment to the healthy trajectory. 
```
scaledPt = pseudotimeScaling(rawDataAML$traj, interDataScaledAML, interDataScaledHAveraged, gAlign)
plot(scaledPt, rawDataAML$traj)
```

## 6. Application of the resulting scaling to show an improvement in markers expression dynamics coherency between the AML and the healthy data
```
#show the improvement in coherency between datasets following scaling:
dfCD11b = data.frame(exp = c(rawDataAML$CD11b, expSetH$CD11b),
                     traj = c(rawDataAML$traj, expSetH$dpt),
                     scaledTraj = c(scaledPt, expSetH$dpt),
                     class = c(rep('AML', nrow(rawDataAML)), rep('Healthy', nrow(expSetH))))
ggplot(dfCD11b, aes(x = traj, y = exp, color = class)) + geom_point(alpha = 0.1) + theme_classic() + ggtitle('before devMap application')
ggplot(dfCD11b, aes(x = scaledTraj, y = exp, color = class)) + geom_point(alpha = 0.1) + theme_classic() + ggtitle('following devMap application')
```
