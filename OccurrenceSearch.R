library(occCite)
library(rfishbase)
library(robis)
library(dplyr)
library(rgbif)
library(rgdal)
library(raster)
library(ggplot2)
library(viridis)
library(rphylopic)

# Get lists of taxa from FishBase and OBIS -----
# A convenience function
'%!in%' <- function(x,y)!('%in%'(x,y))

# Getting species lists
gadiformListFB <- species_list(Order = "Gadiformes")
gadiformListOBIS <- checklist("Gadiformes")
gadiformListOBIS <- gadiformListOBIS %>% filter(taxonRank == "Species", taxonomicStatus == "accepted")
inFBnotOBIS <- gadiformListFB[gadiformListFB %!in% gadiformListOBIS$acceptedNameUsage]
inOBISnotFB <- gadiformListOBIS[gadiformListOBIS$acceptedNameUsage %!in% gadiformListFB,]
inBoth <- gadiformListFB[gadiformListFB %in% gadiformListOBIS$acceptedNameUsage]

scombriformListFB <- species_list(Family = c("Amarsipidae", "Nomeidae", "Ariommatidae", 
                                           "Stromateidae", "Pomatomidae", "Centrolophidae", 
                                           "Icosteidae", "Arripidae", "Tetragonuridae", 
                                           "Chiasmodontidae", "Scombridae", "Caristiidae", 
                                           "Bramidae", "Scombrolabracidae", "Scombropidae", 
                                           "Gempylidae", "Trichiuridae")) # From Betancur-R et al 2017; Eschmeyer was missing Scombropidae
scombriformListOBIS <- checklist(c("Amarsipidae", "Nomeidae", "Ariommatidae", 
                                   "Stromateidae", "Pomatomidae", "Centrolophidae", 
                                   "Icosteidae", "Arripidae", "Tetragonuridae", 
                                   "Chiasmodontidae", "Scombridae", "Caristiidae", 
                                   "Bramidae", "Scombrolabracidae", "Scombropidae", 
                                   "Gempylidae", "Trichiuridae"))
scombriformListOBIS <- scombriformListOBIS %>% filter(taxonRank == "Species", taxonomicStatus == "accepted")
scombInFBnotOBIS <- scombriformListFB[scombriformListFB %!in% scombriformListOBIS$acceptedNameUsage]
scombInOBISnotFB <- scombriformListOBIS[scombriformListOBIS$acceptedNameUsage %!in% scombriformListFB,]
scombInBoth <- scombriformListFB[scombriformListFB %in% scombriformListOBIS$acceptedNameUsage]

belonListFB <- species_list(Family = c("Scomberesocidae",	"Belonidae", "Hemiramphidae", 
                                       "Zenarchopteridae", "Exocoetidae", "Adrianichthyidae"))
belonListOBIS <- checklist(c("Scomberesocidae",	"Belonidae", "Hemiramphidae", 
                              "Zenarchopteridae", "Exocoetidae", "Adrianichthyidae")) # Eschmeyer and Betancur-R agree (low boostrap support for the group in BR, though)
belonListOBIS <- belonListOBIS %>% filter(taxonRank == "Species", taxonomicStatus == "accepted")
belonInFBnotOBIS <- belonListFB[belonListFB %!in% belonListOBIS$acceptedNameUsage]
belonInOBISnotFB <- belonListOBIS[belonListOBIS$acceptedNameUsage %!in% belonListFB,]
belonInBoth <- belonListFB[belonListFB %in% belonListOBIS$acceptedNameUsage]

# Next, hand-curated lists, compared to Eschmeyer's Catalog, accessed September 2020

# Filter list for species found in the Atlantic -----
masterList <- read.csv("data/TaxonomicResolution.csv")

# Atlantic distribution according to Fishbase
fishbaseDistributions <- distribution(masterList$FBName)
fishbaseDistributions <- table(fishbaseDistributions[,c("Species", "FAO")])
fishbaseDistributions <- fishbaseDistributions[,c("Atlantic, Antarctic", "Atlantic, Eastern Central", 
                                                  "Atlantic, Northeast", "Atlantic, Northwest", 
                                                  "Atlantic, Southeast", "Atlantic, Southwest", 
                                                  "Atlantic, Western Central", 
                                                  "Mediterranean and Black Sea")] # Area of interest
fishbaseDistributions <- apply(fishbaseDistributions, 1, sum)
fbATLPresent <- fishbaseDistributions[masterList$FBName] > 0 # Presence/absence in area of interest
masterList <- cbind(masterList,fbATLPresent)

# Atlantic distribution according to OBIS
faoShapefile <- readOGR("data/FAO Fishing Areas 2005/FAO_AREAS.shp")
atlanticShapefile <- aggregate(faoShapefile[faoShapefile$OCEAN=="Atlantic",], dissolve = T)
OCEAN <- "Atlantic"
atlanticShapefile <- SpatialPolygonsDataFrame(atlanticShapefile, as.data.frame(OCEAN))
writeOGR(obj = atlanticShapefile, dsn = "data/", layer = "atlantic", driver = "ESRI Shapefile")
atlanticSimp <- rgeos::gSimplify(atlanticShapefile, tol = 3)
atlanticWKT <- wicket::sp_convert(atlanticSimp, group = T)
atlSpp <- NULL
index <- 1
# Because OBIS gets angry when you ask for more than 100 names at a time
while(length(unique(masterList$OBISChecklist)) > index){
  atlSpp <- append(atlSpp, checklist(scientificname = unique(masterList$OBISChecklist)[index:(index + 99)], 
                                     geometry = atlanticWKT)$scientificName)
  print(index)
  index <- index + 100
}
atlSpp <- gsub(pattern = "(\\w+\\s\\w+)(\\s\\w+)", atlSpp, perl = T, replacement = "\\1") # Toss subspecific epithets
OBIS_atlPresent <- masterList$OBISChecklist %in% atlSpp
masterList <- cbind(masterList, OBIS_atlPresent)

fbOBagree <- masterList$fbATLPresent == masterList$OBIS_atlPresent
masterList <- cbind(masterList, fbOBagree)
write.csv(masterList, "data/taxaWithinAreaOfInterest.csv", row.names = F)

# At this point, I did a manual rectification of the list in cases where there was no data or sources disagreed
# This list was finalized as of Sept 1, 2020

# Getting GBIF taxonomy to verify
atlChecklist <- read.csv("data/taxaWithinAreaOfInterest.csv")
GBIFnames <- vector(mode = "character", length = nrow(atlChecklist))
for(i in 1:nrow(atlChecklist)){
  nameCheck <- as.character(atlChecklist$EschmeyerName[i])
  GBIFname <- name_lookup(query = nameCheck, limit = 1)
  if(length(GBIFname$data)==0){
    GBIFnames[[i]] <- NA
  } else{
    GBIFnames[[i]] <- GBIFname$data$scientificName
  }
}
atlChecklist <- cbind(atlChecklist, GBIFnames)
write.csv(atlChecklist, "data/taxaWithinAreaOfInterest.csv", row.names = F)

# Visualizing the checklist ----
atlChecklist <- read.csv("data/taxaWithinAreaOfInterest.csv")

# All taxa
Source <- c(rep("Fishbase", 2), rep("OBIS", 2), rep("Curated", 2))
Location <- rep(c("Atlantic", "Not Atlantic"), 3)
Count <- c(length(unique(atlChecklist$FBName[atlChecklist$fbATLPresent])),
           length(unique(atlChecklist$FBName[!atlChecklist$fbATLPresent])),
           length(unique(atlChecklist$OBISChecklist[atlChecklist$OBIS_atlPresent])),
           length(unique(atlChecklist$OBISChecklist[!atlChecklist$OBIS_atlPresent])), 
           length(unique(atlChecklist$EschmeyerName[atlChecklist$isAtlantic])),
           length(unique(atlChecklist$EschmeyerName[!atlChecklist$isAtlantic])))
data <- data.frame(Source, Location, Count)

ggplot(data, aes(fill=Location, y=Count, x=Source)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Checklist of Cods, Tunas, Flyingfishes, and Their Allies") +
  scale_fill_brewer(palette="Dark2") +
  theme_minimal(base_size = 18)

# Gadiformes
gadChecklist <- atlChecklist[atlChecklist$Group == "G",]
Source <- c(rep("Fishbase", 2), rep("OBIS", 2), rep("Curated", 2))
Location <- rep(c("Atlantic", "Not Atlantic"), 3)
Count <- c(length(unique(gadChecklist$FBName[gadChecklist$fbATLPresent])),
           length(unique(gadChecklist$FBName[!gadChecklist$fbATLPresent])),
           length(unique(gadChecklist$OBISChecklist[gadChecklist$OBIS_atlPresent])),
           length(unique(gadChecklist$OBISChecklist[!gadChecklist$OBIS_atlPresent])), 
           length(unique(gadChecklist$EschmeyerName[gadChecklist$isAtlantic])),
           length(unique(gadChecklist$EschmeyerName[!gadChecklist$isAtlantic])))
data <- data.frame(Source, Location, Count)
cod <- image_data("bba1800a-dd86-451d-a79b-c5944cfe5231", size = 256)[[1]]

ggplot(data, aes(fill=Location, y=Count, x=Source)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Checklist of Cods and Their Allies") +
  scale_fill_brewer(palette="Dark2") +
  theme_minimal(base_size = 18) +
  add_phylopic(cod, alpha = 1, x = 3, y = sum(data[1:2,3])-27.5+5, ysize = 55)

# Scombriformes
scombChecklist <- atlChecklist[atlChecklist$Group == "S",]
Source <- c(rep("Fishbase", 2), rep("OBIS", 2), rep("Curated", 2))
Location <- rep(c("Atlantic", "Not Atlantic"), 3)
Count <- c(length(unique(scombChecklist$FBName[scombChecklist$fbATLPresent])),
           length(unique(scombChecklist$FBName[!scombChecklist$fbATLPresent])),
           length(unique(scombChecklist$OBISChecklist[scombChecklist$OBIS_atlPresent])),
           length(unique(scombChecklist$OBISChecklist[!scombChecklist$OBIS_atlPresent])), 
           length(unique(scombChecklist$EschmeyerName[scombChecklist$isAtlantic])),
           length(unique(scombChecklist$EschmeyerName[!scombChecklist$isAtlantic])))
data <- data.frame(Source, Location, Count)
tuna <- image_data("16989401-0080-4502-828d-e85a45a262be", size = 256)[[1]]

ggplot(data, aes(fill=Location, y=Count, x=Source)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Checklist of Tunas and Their Allies") +
  scale_fill_brewer(palette="Dark2") +
  theme_minimal(base_size = 18) +
  add_phylopic(tuna, alpha = 1, x = 3, y = sum(data[5:6,3])-(35/2) + 5, ysize = 35)

# Beloniformes
belonChecklist <- atlChecklist[atlChecklist$Group == "B",]
Source <- c(rep("Fishbase", 2), rep("OBIS", 2), rep("Curated", 2))
Location <- rep(c("Atlantic", "Not Atlantic"), 3)
Count <- c(length(unique(belonChecklist$FBName[belonChecklist$fbATLPresent])),
           length(unique(belonChecklist$FBName[!belonChecklist$fbATLPresent])),
           length(unique(belonChecklist$OBISChecklist[belonChecklist$OBIS_atlPresent])),
           length(unique(belonChecklist$OBISChecklist[!belonChecklist$OBIS_atlPresent])), 
           length(unique(belonChecklist$EschmeyerName[belonChecklist$isAtlantic])),
           length(unique(belonChecklist$EschmeyerName[!belonChecklist$isAtlantic])))
data <- data.frame(Source, Location, Count)
flyingFish <- image_data("b6626eff-b00e-428a-b3fa-51679a0cfaa2", size = 256)[[1]]

ggplot(data, aes(fill=Location, y=Count, x=Source)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Checklist of Flyingfishes and Their Allies") +
  scale_fill_brewer(palette="Dark2") +
  theme_minimal(base_size = 18) +
  add_phylopic(flyingFish, alpha = 1, x = 3, y = sum(data[5:6,3])-7.5, ysize = 15)

# Query GBIF for occurrence data ----
atlChecklist <- read.csv("data/taxaWithinAreaOfInterest.csv", stringsAsFactors = F)
atlChecklist <- atlChecklist[atlChecklist$isAtlantic,]
login <- GBIFLoginManager(user = "hannah0wens",
                          email = "hannah.owens@gmail.com",
                          pwd = "Llab7a3m!");

GBIFsearchList <- c(atlChecklist$GBIFnames, atlChecklist$GBIF.Synonym)
GBIFsearchList <- unique(GBIFsearchList[!is.na(GBIFsearchList)])
GBIFsearchList <- studyTaxonList(GBIFsearchList, datasources = "GBIF Backbone Taxonomy")

myBridgeTreeObject <- occQuery(x = GBIFsearchList, GBIFLogin = login, datasources = "gbif", 
                               GBIFDownloadDirectory = "data/GBIFDownloads/")
saveRDS(myBridgeTreeObject, file = "data/GBIFDownloads/myBridgeTreeObject")
myBridgeTreeObject <- readRDS(file = "data/GBIFDownloads/myBridgeTreeObject")
#Get and write citations
myOccCitations <- occCitation(myBridgeTreeObject)
rawCitations <- file("data/rawCitations.txt")
writeLines(print(myOccCitations), rawCitations)
close(rawCitations)
