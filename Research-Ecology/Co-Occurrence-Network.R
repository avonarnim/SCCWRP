#Project 5 - creating co-occurance networks
#Adam von Arnim
require("plyr")
require(dplyr)
require(geosphere)
require(tidyr)
library(igraph)
library(network)
library(data.table)
library(netassoc)

#reads all data
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

masterComInput <- masterComInput[, colSums(masterComInput != 0) > 0] #gets rid of empty columns

#creates network
net <- make_netassoc_network(masterComInput, vegan::permatfull(masterComInput,fixedmar="both",mtype="prab",times=100)$perm[[1]],method="partial_correlation",args=list(method="shrinkage"),p.method='fdr', numnulls=1000, plot=FALSE,alpha=1e-4,verbose=FALSE)
dev.off()
networkgraph <- as.undirected(net$network_all,mode="collapse") #Generate full graph from co-occurrence patterns.
networkgraph_pos <- as.undirected(net$network_pos,mode="collapse") #Generate positive edge graph from co-occurrence patterns.
networkgraph_neg <- as.undirected(net$network_neg,mode="collapse") #Generate negative edge graph from co-occurrence patterns.
plot_netassoc_network(networkgraph)

