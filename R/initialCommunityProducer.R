initialCommunityProducer <- function(speciesLayers, speciesPresence, studyArea, rstStudyArea) {
  specieslayerInStudyArea <- crop(speciesLayers, studyArea)  ## TODO: Ceres: this is probably no longer necessary
  #if(isTRUE(tryCatch(getCluster(), error=function(x) TRUE, silent=TRUE))) beginCluster()
  specieslayerInStudyArea <- specieslayerInStudyArea * rstStudyArea  ## TODO: Ceres: this is probably no longer necessary
  names(specieslayerInStudyArea) <- names(speciesLayers)  ## TODO: Ceres: this is probably no longer necessary

  #specieslayerInStudyArea <- specieslayerInStudyArea*(!is.na(rstStudyArea))
  # specieslayerInStudyArea <- suppressWarnings(fastMask(specieslayerInStudyArea, #                                                  studyArea))
  speciesNames <- names(specieslayerInStudyArea)[which(maxValue(specieslayerInStudyArea) >= speciesPresence)]   ## TODO: Ceres: this is probably no longer necessary
  specieslayerBySpecies <- raster::subset(specieslayerInStudyArea, speciesNames[1])    ## TODO: Ceres: this is probably no longer necessary
  specieslayerBySpecies[which(is.na(specieslayerBySpecies[]) & specieslayerBySpecies[] <= 5)] <- 0    ##  Ceres: this is weird, shouldn't this be OR?
  # specieslayerBySpecies[Which(is.na(specieslayerBySpecies) & specieslayerBySpecies<=5,
  #                             cells = TRUE)] <- 0 # 5% or less presence removed
  speciesComMap <- as.logical(specieslayerBySpecies)
  rm(specieslayerBySpecies)
  k <- 1
  for (species in speciesNames[2:length(speciesNames)]) {
    specieslayerBySpecies <- raster::subset(specieslayerInStudyArea, species)
    specieslayerBySpecies[which(is.na(specieslayerBySpecies[]) & specieslayerBySpecies[] <= 5)] <- 0
    # specieslayerBySpecies[Which(is.na(specieslayerBySpecies) & specieslayerBySpecies <= 5,
    #                             cells = TRUE)] <- 0
    speciesMap <- as.logical(specieslayerBySpecies)
    speciesComMap <- speciesMap*(10^k) + speciesComMap
    k <- k + 1
    rm(specieslayerBySpecies, speciesMap)
  }
  # set the non-forested area as NA
  #speciesComMap1 <- speciesComMap
  #speciesComMap[Which(speciesComMap == 0, cells = TRUE, na.rm = FALSE)] <- NA
  speciesComMap[which(speciesComMap[] == 0)] <- NA

  initialCommunities <- data.table(mapcode = sort(unique(getValues(speciesComMap))))
  initialCommunities[, mapCodeStr := as.character(mapcode)]
  initialCommunities[, NofStr := nchar(mapCodeStr)]   ## shouldn't 10 and 11 be distinguished as one refers to dominance and the other co-dominance?
  for (i in 1:(length(speciesNames) - 1)) {
    initialCommunities[NofStr == i, mapCodeFull := paste(paste(rep("0",(length(speciesNames) - i)),
                                                               collapse = ""),
                                                         mapCodeStr,
                                                         sep = "")]
  }
  initialCommunities[NofStr == length(speciesNames), mapCodeFull := mapCodeStr]
  output <- data.table(mapcode = numeric(), speciesPresence = character(),
                       species = character())
  for (i in 1:nrow(initialCommunities)) {
    outputAdd <- data.table(mapcode = initialCommunities$mapcode[i],
                            speciesPresence = substring(initialCommunities$mapCodeFull[i],
                                                        seq(1, length(speciesNames), 1),
                                                        seq(1, length(speciesNames), 1)),   ## why not simply use strsplit?
                            species = speciesNames[length(speciesNames):1])

    output <- rbind(output, outputAdd)
  }
  initialCommunities <- output[speciesPresence != "0",]
  initialCommunities[,newMapCode := as.numeric(as.factor(mapcode))]
  mapcodeconnection <- unique(initialCommunities[, .(mapcode, newMapCode)], by = "mapcode")
  indexTable <- data.table(pixelIndex = 1:ncell(speciesComMap),
                           mapcode = getValues(speciesComMap))
  indexTable <- indexTable[!is.na(mapcode), ]
  indexTable <- setkey(indexTable, mapcode)[setkey(mapcodeconnection, mapcode),
                                            nomatch = 0]
  speciesComMap[indexTable$pixelIndex] <- indexTable$newMapCode
  initialCommunities[, ':='(mapcode = newMapCode, newMapCode = NULL, speciesPresence = NULL)]
  
  return(list(initialCommunityMap = speciesComMap,
              initialCommunity = initialCommunities))
}
