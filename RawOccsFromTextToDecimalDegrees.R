library(sp)

#melo2010
melo2010 <- read.csv("data/OccurrencesFromLitAndOtherDatabases/Melo2010Raw.csv", header = T)
newLong <- char2dms(melo2010$Long, chd = "d", chm = "m", chs = "s")
melo2010$Long <- as.numeric(newLong)
newLat <- char2dms(melo2010$Lat, chd = "d", chm = "m", chs = "s")
melo2010$Lat <- as.numeric(newLat)

depth <- melo2010$Depth
depth <- strsplit(depth, split = "m")
depthUncertainty <- rep(times = length(depth), NA)
singleDepth <- rep(times = length(depth), NA)
for(x in 1:length(depth)){
  if(length(depth[[x]]) == 0){
    singleDepth[[x]] <- NA
    depthUncertainty[[x]] <- NA
  } else if(grepl("–", depth[[x]])){
    minDepth <- as.numeric(unlist(strsplit(split = "–", depth[[x]])))[1]
    maxDepth <- as.numeric(unlist(strsplit(split = "–", depth[[x]])))[2]
    depthUncertainty[[x]] <- (maxDepth - minDepth)/2
    singleDepth[[x]] <- as.numeric(depthUncertainty[[x]]) + minDepth
  } else {
    singleDepth[[x]] <- depth[[x]]
    depthUncertainty[[x]] <- NA
  }
}
melo2010$Depth <- singleDepth
melo2010$DepthUncertainty <- depthUncertainty
write.csv(melo2010, file = "data/OccurrencesFromLitAndOtherDatabases/melo2010Processed.csv", row.names = F)


#pires2019
pires2019 <- read.csv("data/OccurrencesFromLitAndOtherDatabases/Pires2019Raw.csv", header = T)
newLong <- char2dms(pires2019$Long, chd = "d", chm = "m", chs = "s")
pires2019$Long <- as.numeric(newLong)
newLat <- char2dms(pires2019$Lat, chd = "d", chm = "m", chs = "s")
pires2019$Lat <- as.numeric(newLat)

depth <- pires2019$Depth
depth <- strsplit(depth, split = "m")
depthUncertainty <- rep(times = length(depth), NA)
singleDepth <- rep(times = length(depth), NA)
for(x in 1:length(depth)){
  if(length(depth[[x]]) == 0){
    singleDepth[[x]] <- NA
    depthUncertainty[[x]] <- NA
  } else if(grepl("–", depth[[x]])){
    minDepth <- as.numeric(unlist(strsplit(split = "–", depth[[x]])))[1]
    maxDepth <- as.numeric(unlist(strsplit(split = "–", depth[[x]])))[2]
    depthUncertainty[[x]] <- (maxDepth - minDepth)/2
    singleDepth[[x]] <- as.numeric(depthUncertainty[[x]]) + minDepth
  } else {
    singleDepth[[x]] <- depth[[x]]
    depthUncertainty[[x]] <- NA
  }
}
pires2019$Depth <- singleDepth
pires2019$DepthUncertainty <- depthUncertainty
write.csv(pires2019, file = "data/OccurrencesFromLitAndOtherDatabases/pires2019Processed.csv", row.names = F)
