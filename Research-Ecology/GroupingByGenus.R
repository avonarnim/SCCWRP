#Project 2 - Make a presence/absence matrix for genuses
#Adam von Arnim 11/5/2019
library("plyr")
library(dplyr)
library(tidyr)

#reads all data
communityInput <- read.delim("~/Downloads/18SV9P2TableWithTaxonomy (1).txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
commInput2 <- read.delim("~/Downloads/18SV9P1TableWithTaxonomy.txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
metadata <- read.csv("~/Downloads/MetagenomicSampleSiteMetadata (1).csv")
stCodes <- read.delim("~/Downloads/SampleStationCodesID.txt")

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




