#Project 1 - Make a presence/absence matrix for data
#Adam von Arnim 10/19/2019
library(zetadiv)
library("plyr")
library(dplyr)
communityInput <- read.delim("~/Downloads/18SV9P2TableWithTaxonomy (1).txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
commInput2 <- read.delim("~/Downloads/18SV9P1TableWithTaxonomy.txt", header=T, as.is=T,skip=0,fill=TRUE,check.names=FALSE,encoding="UTF-8")
sampleIDs <- read.delim("~/Downloads/SampleStationCodesID.txt")
metadata <- read.csv("~/Downloads/MetagenomicSampleSiteMetadata (1).csv")
metadata <- merge(sampleIDs, metadata)

communityInput <- dplyr::rename(communityInput, OTUID = `OTU ID`)
commInput2 <- dplyr::rename(commInput2, OTUID = `OTU ID`)
allComInput <- as.data.frame(unique(c(communityInput$OTUID, commInput2$OTUID)))
names(allComInput) <- "OTUID"
communityInput$OTUID <- as.character(communityInput$OTUID)

communityInput$DNAStandard <- NULL
communityInput$`Ext-Blank1` <- NULL
communityInput$`Ext-Blank2` <- NULL
communityInput$FB <- NULL
communityInput$NTC <- NULL
communityInput$SNAStandardII <- NULL
communityInput$ConsensusLineage <- NULL
commInput2$DNAstandardI <- NULL
commInput2$DNAstandardII <- NULL
commInput2$ConsensusLineage <- NULL

allComInput <- left_join(allComInput, communityInput, by=c("OTUID"))
allComInput <- left_join(allComInput, commInput2, by=c("OTUID"))

row.names(communityInput) <- communityInput[, 1]
communityInput$OTUID <- NULL
communityInput[communityInput > 0] <- 1
row.names(commInput2) <- commInput2[,1]
commInput2$OTUID <- NULL

row.names(allComInput) <- allComInput[,1]
allComInput$`Row.names` <- NULL
allComInput$`230.x` <- NULL
allComInput$`230.y` <- NULL

metadata <- metadata[which(metadata$SampleNum %in% as.numeric(colnames(allComInput))),]
metadata <- arrange(metadata, SampleNum)
allComInput <- allComInput[, as.character(metadata$SampleNum)]



data.spec <- as.data.frame(t(communityInput))
data.xy = metadata[,c("Latitude", "Longitude")]
zetaDistance <- Zeta.ddecay(xy=data.xy, data.spec=data.spec, distance.type="ortho", order=3)

zetaDecay <-  Zeta.decline.ex(communityInput, orders=1:10)
dat <- data.frame()
dat[1,1] <- zetaDecay$zeta.exp$coefficients[1]  #zeta exponential decay intercept
dat[1,2] <- zetaDecay$zeta.exp$coefficients[2]  #zeta exponential decay exponent
dat[1,3] <- zetaDecay$zeta.pl$coefficients[1]  #zeta power law exponent
dat[1,4] <- zetaDecay$zeta.pl$coefficients[2]  #zeta power law exponent
dat