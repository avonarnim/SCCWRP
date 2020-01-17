#Project 4 - subsetting clusters and running zetadiv
#Adam von Arnim
library(zetadiv)
library("plyr")
library(dplyr)
library(tidyr)
library("geosphere")
library(raster)

#reads all data
metadata <- read.csv("~/Downloads/MetagenomicSampleSiteMetadata (1).csv")
stCodes <- read.delim("~/Downloads/SampleStationCodesID.txt")
communityInput <- read.delim("~/Downloads/18SV9P2TableWithTaxonomy (1).txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
commInput2 <- read.delim("~/Downloads/18SV9P1TableWithTaxonomy.txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")

#creates filter for only data from desired OTUIDs
communityInput <- dplyr::rename(communityInput, OTUID = `OTU ID`)
commInput2 <- dplyr::rename(commInput2, OTUID = `OTU ID`)
allComInput <- as.data.frame(unique(c(communityInput$OTUID, commInput2$OTUID)))
names(allComInput) <- "OTUID"
communityInput$OTUID <- as.character(communityInput$OTUID)
#gets rid of unnecessary information
communityInput$DNAStandard <- NULL
communityInput$`Ext-Blank1` <- NULL
communityInput$`Ext-Blank2` <- NULL
communityInput$FB <- NULL
communityInput$NTC <- NULL
communityInput$SNAStandardII <- NULL
commInput2$DNAstandardI <- NULL
commInput2$DNAstandardII <- NULL
#combines the two plates
allComInput <- left_join(allComInput, communityInput, by=c("OTUID"))
allComInput <- left_join(allComInput, commInput2, by=c("OTUID"))

#Get a list of unique OTU names from count tables and convert to a data frame.
#uniqueOTUs <- as.data.frame(unique(c(communityInputRawPlate1$OTUID,communityInputRawPlate2$OTUID)))
uniqueOTUs <- rbind(communityInput[,c("OTUID","ConsensusLineage")],commInput2[,c("OTUID","ConsensusLineage")])
#colnames(uniqueOTUs) <- c("OTUID")
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)
#Split OTU names into Domain through Genus+Species.
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Domain", "Kingdom","Phylum","Class","Order","Family","GenusSpecies"),sep=";", extra="drop"))
uniqueOTUs$GenusSpecies <- trimws(uniqueOTUs$GenusSpecies,which="left") #Remove starting blank space from genus names
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'GenusSpecies',c("Genus","Species"),sep=" ", extra="warn"))
uniqueOTUs <- uniqueOTUs[!duplicated(uniqueOTUs$OTUID),]
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Domain!="Unassigned" & uniqueOTUs$Genus!="Ambiguous_taxa",]

uniqueOTUs$Genus <- gsub("^.*__", "", uniqueOTUs$Genus) #gsub("D_([0-9)]__", ""\", ConsensusLineage)

masterComInput <- left_join(uniqueOTUs, communityInput, by="OTUID")
masterComInput <- left_join(masterComInput, commInput2, by="OTUID")
masterComInput$ConsensusLineage.x <- NULL
masterComInput$ConsensusLineage.y <- NULL

taxLevels <- c("Domain", "Kingdom","Phylum","Class","Order","Family","Genus", "Species", "OTUID")
taxLevel <- c("Genus")
taxIgnore <- taxLevels[taxLevels != taxLevel]

masterComInput <- masterComInput[, -which(names(masterComInput) %in% taxIgnore)]

#sums up the presence frequency data from all rows with a shared Genus
masterComInput <- aggregate(.~Genus, masterComInput, FUN=sum, na.action=na.omit)

#coverts to presence/absence
#must first assign GenusSpecies info to row names
row.names(masterComInput) <- masterComInput$Genus
masterComInput$Genus <- NULL
masterComInput[masterComInput <= 2] <- 0
masterComInput[masterComInput > 2] <- 1

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

metadata <- metadata[!duplicated(metadata$SampleNum),]

#finds Haversine distances between coordinates
as.data.frame(distances)
distances <- distm(metadata[,c("Longitude","Latitude")]) #haversine comes automatically when inputting long/lat
#used to hold all variances
variances <- c()

#finds minimum variance threshold
for (i in 1:50)
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

#each row is a separate test run, each column is one of the coefficients + one of the variables
aggregateData <- c()

for (j in 1:10){
  for (i in unique(clusters$cluster))
  {
      #take a random sample from the clusters, subset metadata/genus presence appropriately
      clusterIndices <- which(clusters[["cluster"]]==i)
      clusterGroup <- sample(clusterIndices, size=10)
      sampleMetadata <- metadata[clusterGroup,]
      sampleComInput <- masterComInput[,which(names(masterComInput) %in% sampleMetadata$SampleNum)]
      
      data.spec <- as.data.frame(t(sampleComInput))
      data.xy = sampleMetadata[,c("Latitude", "Longitude")]
      zetaDistance <- Zeta.ddecay(xy=data.xy, data.spec=data.spec, distance.type="ortho", order=3)
  
      zetaDecay <-  Zeta.decline.ex(sampleComInput, orders=1:6)
      dat <- data.frame()
      dat[1,1] <- zetaDecay$zeta.exp$coefficients[1]  #zeta exponential decay intercept
      dat[1,2] <- zetaDecay$zeta.exp$coefficients[2]  #zeta exponential decay exponent
      dat[1,3] <- zetaDecay$zeta.pl$coefficients[1]  #zeta power law exponent
      dat[1,4] <- zetaDecay$zeta.pl$coefficients[2]  #zeta power law exponent
      dat[1,5] <- mean(sampleMetadata$LU, na.rm=TRUE)
      dat[1,6] <- mean(sampleMetadata$MaxN, na.rm=TRUE)
      dat[1,7] <- mean(sampleMetadata$MaxOrthoP, na.rm=TRUE)
      dat[1,8] <- mean(sampleMetadata$MaxP, na.rm=TRUE)
      dat[1,9] <- mean(sampleMetadata$site_elev, na.rm=TRUE)
      dat[1,10] <- length(unique(sampleMetadata$HUC8)) #number watersheds
      dat[1,11] <- zetaDecay$aic$AIC[1]  #aic exponential decay coefficient
      dat[1,12] <- zetaDecay$aic$AIC[2]  #aic power law coefficient
      #dat[1,13] <- mean(distance[])    #mean distance
      
      aggregateData <- rbind(aggregateData,dat)
  }
}
#make graphs showing the zeta diversity compared to other variables

