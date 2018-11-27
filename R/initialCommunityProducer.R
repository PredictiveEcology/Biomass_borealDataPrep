initialCommunityProducer <- function(speciesLayers, speciesPresence, rstStudyArea, standAgeMap) {
  sppNames <- names(speciesLayers)
  speciesLayers <- speciesLayers * rstStudyArea  ## TODO: Ceres: this is probably no longer necessary
  names(speciesLayers) <- sppNames

  speciesNames <- names(speciesLayers)[which(maxValue(speciesLayers) >= speciesPresence)]

  k <- 0
  n <- length(speciesNames)
  digits <- 2
  mapCodes <- rep("", ncell(speciesLayers[[1]]))

  for (species in speciesNames[1:n]) { ## TODO: use raster::calc or similar
    specieslayerBySpecies <- raster::subset(speciesLayers, species)
    specieslayerBySpecies[which(is.na(specieslayerBySpecies[]) |
                                  specieslayerBySpecies[] <= 5)] <- 0L ## why 5% cover?
    speciesCodes <- paddedFloatToChar(round(specieslayerBySpecies[], -1) / 10, padL = digits, padR = 0)

    mapCodes <- paste0(mapCodes, speciesCodes)
    k <- k + 1
    rm(speciesCodes)
  }

  ## append stand ages
  ageMap <- standAgeMap[]
  ageMap[which(is.na(ageMap))] <- 0L
  ageCodes <- paddedFloatToChar(round(ageMap, -1) / 10, padL = digits, padR = 0)
  mapCodes <- paste0(mapCodes, ageCodes)
  rm(ageCodes)

  ## convert mapCodes for use with data.table and as integer raster
  mapCodesInt <- strtoi(mapCodes, base = 10) ## much faster than as.integer
  ids <- which(mapCodesInt == 0L)
  mapCodesInt[ids] <- NA_integer_
  mapCodes[ids] <- NA_character_
  speciesComMap <- speciesLayers[[1]]
  speciesComMap[] <- mapCodesInt

  initialCommunities <- data.table(mapcode = sort(unique(mapCodesInt)))
  initialCommunities[, mapCodeStr := sort(unique(mapCodes))]
  output <- data.table(mapcode = numeric(), speciesPresence = character(),
                       species = character(), age = numeric())
  for (i in 1:nrow(initialCommunities)) {
    outputAdd <- data.table(mapcode = initialCommunities$mapcode[i],
                            speciesPresence = substring(initialCommunities$mapCodeStr[i],
                                                        seq(1, digits * n, 2),
                                                        seq(2, digits * n, 2)),
                            species = speciesNames[1:length(speciesNames)],
                            age1 = strtoi(rep(substring(initialCommunities$mapCodeStr[i],
                                                        2 + digits * n - 1,
                                                        2 + digits * n), 3),
                                          base = 10) * 10)

    output <- rbind(output, outputAdd)
  }
  initialCommunities <- output[speciesPresence != "00",]
  initialCommunities[, newMapCode := as.numeric(as.factor(mapcode))]
  mapcodeconnection <- unique(initialCommunities[, .(mapcode, newMapCode)], by = "mapcode")
  indexTable <- data.table(pixelIndex = 1:ncell(speciesComMap),
                           mapcode = getValues(speciesComMap))
  indexTable <- indexTable[!is.na(mapcode), ]
  indexTable <- setkey(indexTable, mapcode)[setkey(mapcodeconnection, mapcode), nomatch = 0]
  speciesComMap[indexTable$pixelIndex] <- indexTable$newMapCode
  initialCommunities[, ':='(mapcode = newMapCode, newMapCode = NULL, speciesPresence = NULL)]

  return(list(initialCommunityMap = speciesComMap,
              initialCommunity = initialCommunities))
}
