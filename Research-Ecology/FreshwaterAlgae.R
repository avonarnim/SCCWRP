#Project 6 - filtering the 18SV9 OTU table to only contain freshwater algae reads
#Adam von Arnimn
require("plyr")
require(dplyr)
require(tidyr)
require(naniar)
require(taxize)

options(ENTREZ_KEY="72b83e2d468969a5c9f1870be0806b8b4308")

#reads in 18SV9 tables
communityInput <- read.delim("~/Downloads/18SV9P2TableWithTaxonomy (1).txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
commInput2 <- read.delim("~/Downloads/18SV9P1TableWithTaxonomy.txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")

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

#Get a list of unique OTU names
uniqueOTUs <- rbind(communityInput[,c("OTUID","ConsensusLineage")],commInput2[,c("OTUID","ConsensusLineage")])
uniqueOTUs$OTUID <- as.character(uniqueOTUs$OTUID)
#Remove superfluous strings from taxonomic labels.
uniqueOTUs$ConsensusLineage <- gsub("D_[0-9]+__","",uniqueOTUs$ConsensusLineage)
uniqueOTUs$ConsensusLineage <- gsub("g__","",uniqueOTUs$ConsensusLineage)
#Split OTU names into Domain through Genus+Species.
uniqueOTUs$FullTaxonomy <- uniqueOTUs$ConsensusLineage
uniqueOTUs <- suppressWarnings(separate(uniqueOTUs,'ConsensusLineage',c("Rank1", "Rank2","Rank3","Rank4","Rank5","Rank6","Rank7"),sep=";", extra="drop"))
uniqueOTUs$Rank7 <- trimws(uniqueOTUs$Rank7,which="left") #Remove starting blank space from genus names
uniqueOTUs <- uniqueOTUs[!duplicated(uniqueOTUs$OTUID),]
#Filter out ambiguous taxonomies
uniqueOTUs <- uniqueOTUs[uniqueOTUs$Rank1!="Unassigned" & uniqueOTUs$Rank1!="Ambiguous_taxa",]
ambiguousList <- c("Incertae Sedis","metagenome","sp.","environmental","eukaryote","uncultured","soil","Ambiguous_taxa","group","cf.","aff.","gen.","marine","cf","unidentified","Uncultured","invertebrate")
ambiguousList <- as.list(ambiguousList)
uniqueOTUs <- data.frame(lapply(uniqueOTUs, trimws), stringsAsFactors = FALSE)
uniqueOTUs <- replace_with_na_all(data=uniqueOTUs,condition=~.x %in% as.list(ambiguousList))
#Get the furthest resolved taxonomic level.
uniqueOTUs$LeafTaxa <- apply(uniqueOTUs[,!colnames(uniqueOTUs) %in% c("FullTaxonomy")], 1, function(x) tail(na.omit(x), 1))

#combines the two plates
allComInput <- left_join(uniqueOTUs, communityInput, by="OTUID")
allComInput <- left_join(allComInput, commInput2, by="OTUID")
allComInput$ConsensusLineage.x <- NULL
allComInput$ConsensusLineage.y <- NULL

#reads in algae table
algaeList <- read.table("~/Downloads/algae_STE.csv", header=TRUE, sep=",",as.is=T,skip=0,fill=TRUE,quote="",check.names=FALSE, encoding = "UTF-8",na.strings=c("","NA"))
algaeList <- algaeList[algaeList$Kingdom =="Eukaryota",]
algaeList <- algaeList[, !duplicated(colnames(algaeList))]
algaeList <- dplyr::rename(algaeList,LeafTaxa = FinalIDassigned)
algaeList$Kingdom <- NULL
algaeList <- algaeList[!duplicated(algaeList$LeafTaxa),]
algaeList$FinalID <- NULL

algaeComInput <- allComInput[which(allComInput$LeafTaxa %in% c(algaeList$LeafTaxa, algaeList$Genus, algaeList$Phylum, algaeList$Class, algaeList$Order, algaeList$Family)), ]
algaePhyla <- na.omit(unique(algaeList$Phylum))

#i need to find all the taxonomies in allComInput where there is a gap between the leaftaxa & the other taxa
#then throw them through the following for loop
#and within the for loop, extract a filler for the gap or just redo the ranking table in allComInput

#subsetComInput gives everything in allComInput which has a matching phylum
subsetComInput <- allComInput %>% filter_at(vars(Rank1,Rank2,Rank3,Rank4,Rank5,Rank6,Rank7), any_vars(. %in% algaePhyla))


#finding the full taxonomies for all organisms which at least meet the phylum criteria
drawingInput <- unique(subsetComInput$LeafTaxa)
#cbind.fill allows for creating a data frame of all taxonomies, held in uniqueClassifications
cbind.fill <- function(...){
  transpoted <- lapply(list(...),t)
  transpoted_df <- lapply(transpoted, as.data.frame)
  return (data.frame(t(rbind.fill(transpoted_df))))
}
uniqueClassifications <- data.frame()
for (name in drawingInput){
  thisClassification <- classification(name, db="ncbi",rows=1)[[name]]
  if (!(is.na(thisClassification)))
  {
    thisClassification <- as.data.frame(thisClassification[["name"]])
    uniqueClassifications <- cbind.fill(uniqueClassifications, thisClassification)
  }
}
uniqueClassifications <- t(as.data.frame(uniqueClassifications))
algaeComInput <- allComInput[which(allComInput$LeafTaxa %in% c(algaeList$LeafTaxa, algaeList$Genus, algaeList$Phylum, algaeList$Class, algaeList$Order, algaeList$Family)), ]

