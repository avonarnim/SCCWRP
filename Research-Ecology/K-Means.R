#Project #3 - Make a k-means clustering for the data
#Adam von Arnim 11/19/2019
library("plyr")
library(dplyr)
library(tidyr)
library("geosphere")


#reads all data
metadata <- read.csv("~/Downloads/MetagenomicSampleSiteMetadata (1).csv")
stCodes <- read.delim("~/Downloads/SampleStationCodesID.txt")

#remove station codes w/ 202/230, merge station codes with metadata
stCodes <- stCodes[!(stCodes$SampleNum=="230" | stCodes$SampleNum=="202"), ]
metadata <- merge(stCodes, metadata)

#does some formatting and creates land usage column
metadata <- suppressWarnings(metadata[which(metadata$SampleNum %in% as.numeric(colnames(masterComInput))),])
metadata <- arrange(metadata, SampleNum)
metadata$elev_range <- as.numeric(metadata$elev_range)
metadata$max_elev <- as.numeric(metadata$max_elev)
metadata$LU <- metadata$Ag_2011_5K+metadata$URBAN_2011_5K+metadata$CODE_21_2011_5K

#takes in n&p data
npLabData <- read.csv("~/Downloads/NandP_labdata.csv")
npLabData <- npLabData %>% mutate_all(na_if,"")
npLabData$MaxOrthoP <- as.numeric(as.character(npLabData$MaxOrthoP))
npLabData$MaxP <- as.numeric(as.character(npLabData$MaxP))

#adds n/p data to metadata
metadata <- merge(metadata, npLabData)

#takes in watershed data
watersheds <- read.csv("~/Downloads/MetagenomicSitesWithWatersheds.csv")
llhuc <- watersheds[, (names(watersheds) %in% c("Latitude", "Longitude", "HUC8"))]
#finds huc6/4
llhuc$HUC6 <- trunc((llhuc$HUC8)/100)
llhuc$HUC4 <- trunc((llhuc$HUC6)/100)
metadata <- merge(metadata, llhuc)

#finds Haversine distances between coordinates
as.data.frame(distances)
distances = distm(metadata[,c("Longitude","Latitude")]) #haversine comes automatically when inputting long/lat
#used to hold all variances
variances <- list()

#finds minimum variance threshold
for (i in 1:1000)
{ clusters <- kmeans(distances, centers=4)
  variances <- c(variances, sd(table(clusters$cluster))/mean(table(clusters$cluster)))
}
threshold = min(variances)
#finds a clustering that is similar to the threshold
repeat{
  clusters <- kmeans(distances, centers=4)
  this_variance = sd(table(clusters$cluster))/mean(table(clusters$cluster))
  if (this_variance <= threshold)
    break
}
