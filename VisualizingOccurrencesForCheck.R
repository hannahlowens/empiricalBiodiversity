library(ggplot2)
library(ggmap)
library(rnaturalearth)
library(dplyr)
library(grid)

# Map occurrences
occMap <- function(occ_dat){
  
  bb_occ <- make_bbox(lon = occ_dat$decimalLongitude, lat = occ_dat$decimalLatitude, f = 0.05)
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  obis_map <- ggplot(world) +
    geom_sf(color = "gray", fill = "gray") +
    coord_sf(xlim = bb_occ[c(1,3)], ylim = bb_occ[c(2,4)]) +
    theme(panel.background = element_rect(fill = "steelblue")) +
    coord_equal() +
    theme(panel.grid = element_blank())
  
  # Now add the occurrence points
  cols <- gsub(x = gsub("OBIS", replacement = "darkorange", x = occ_dat$source), 
               "GBIF", replacement = "darkred")
  
  obis_map <- obis_map + 
    geom_point(data = occ_dat, aes(x = decimalLongitude, y = decimalLatitude),
               colour = cols, shape = 20, alpha = 2/3) + 
    coord_sf(xlim = bb_occ[c(1,3)], ylim = bb_occ[c(2,4)]) +
    xlab("Longitude") +
    ylab("Latitude") +
    ggtitle(paste0(occ_dat$scientificName[1], ", ", nrow(occ_dat), " points"))
  return(obis_map)
}

# Map depth v latitude
occTrend <- function(occ_dat){
  occ_dat <- occ_dat[!(occ_dat$depth==0),]
  cols <- gsub(x = gsub("OBIS", replacement = "darkorange", x = occ_dat$source), 
               "GBIF", replacement = "darkred")
  occtr <- ggplot(occ_dat, aes(x = abs(decimalLatitude), y = depth)) + 
    geom_smooth(method = "lm", colour = "steelblue") +
    geom_point(colour = cols, alpha = 2/3, shape = 20) + 
    ggtitle(paste0(occ_dat$scientificName[1], ", ", nrow(occ_dat), " points")) +
    xlab("Absolute Latitude") +
    ylab("Depth") +
    theme_classic()
  return(occtr)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

pdf(file = "occurrencePreview.pdf")
occs <- list.files(path = "data/MergedOccurrences/", pattern = ".csv", full.names = T)
rawCounts <- vector(mode = "list", length = length(occs))
occCounts <- vector(mode = "list", length = length(occs))
depthCounts <- vector(mode = "list", length = length(occs))
for (sp in occs){
  oc <- read.csv(sp)
  rawCounts[[match(sp, occs)]] <- nrow(oc)
  if (nrow(oc) > 1) oc <- oc %>% distinct(scientificName, decimalLatitude, decimalLongitude, depth, source)
  occCounts[[match(sp, occs)]] <- nrow(oc)
  if (nrow(oc) > 1){
    p1 <- occMap(oc)
    oc <- oc[complete.cases(oc),]
    depthCounts[[match(sp, occs)]] <- nrow(oc)
    p2 <- occTrend(oc)
    multiplot(p1, p2, cols=1)
  }
  print(paste0(oc[1,1], " has ", nrow(oc), " points."))
}
dev.off()

postVisOccs <- cbind(stringr::str_extract(occs, pattern = "\\w*\\s\\w*"),rawCounts,occCounts, depthCounts)
write.csv(postVisOccs, "data/postVisualizationOccurrenceCounts.csv", row.names = F)
