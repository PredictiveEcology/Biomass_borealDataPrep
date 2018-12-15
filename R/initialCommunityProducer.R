initialCommunityProducer <- function(speciesLayers, #speciesPresence,
                                     #rstStudyArea,
                                     standAgeMap,
                                     percentileGrps = 10) {
  sppNames <- names(speciesLayers)
  # speciesLayers <- speciesLayers * rstStudyArea  ## Eliot rm'd this. Don't want NAs here. prev comment: Ceres: this is probably no longer necessary
  # names(speciesLayers) <- sppNames

  speciesNames <- names(speciesLayers)#[which(maxValue(speciesLayers) >= speciesPresence)]

  n <- length(speciesNames)
  pctRound <- 100 / percentileGrps
  digits <- nchar(percentileGrps) # use quintiles, i.e., 100 pct / 20 ==> 1, 2, 3, 4, 5
  mapCodes <- rep("", ncell(speciesLayers[[1]]))
  mapCodesNAs <- which(is.na(speciesLayers[[1]][]))
  mapCodes <- mapCodes[-mapCodesNAs]

  # matrix of all pct covers
  xtileAbundNum <- round(round(speciesLayers[][-mapCodesNAs,], -1) / pctRound, 0) # use
  sppMatrix <- matrix(rep(speciesNames, each = NROW(xtileAbundNum)), ncol = ncol(xtileAbundNum))
  sppMatrix[xtileAbundNum==0] <- ""
  xtileAbund <- if (digits > 1)
    matrix(paddedFloatToChar(xtileAbundNum, digits, 0), ncol = ncol(xtileAbundNum))
  else
    xtileAbundNum
  mapCodesSpp <- apply(xtileAbund, 1, paste, collapse = "")
  sppVector <- apply(sppMatrix, 1, paste, collapse = ".")
  sppVector2 <- gsub(pattern = "\\.{2,}", replacement = "\\.", sppVector)
  sppVector2 <- gsub(pattern = "^\\.", replacement = "", sppVector2)
  sppVector2 <- gsub(pattern = "\\.$", replacement = "", sppVector2)
  #b <- apply(xtileAbund, 1, function(x) paddedFloatToChar(paste(x, collapse = ""), padL = digits, padR = 0))

  # for (species in speciesNames[1:n]) { ## TODO: use raster::calc or similar
  #   specieslayerBySpecies <- speciesLayers[[species]]
  #   #specieslayerBySpecies[#which(is.na(specieslayerBySpecies[]) |
  #   #                             specieslayerBySpecies[] <= 5)] <- 0L ## why 5% cover?
  #   speciesCodes <- paddedFloatToChar(round(round(specieslayerBySpecies[][-mapCodesNAs], -1) / pctRound, 0), padL = digits, padR = 0)
  #
  #   mapCodes <- paste0(mapCodes, speciesCodes)
  #   rm(speciesCodes)
  # }
  ## append stand ages
  ageMap <- standAgeMap[]
  #ageMapNAs <- which(is.na(ageMap))

  maxAgeOnMap <- max(ageMap[-mapCodesNAs], na.rm = TRUE)
  # ageMap[which(is.na(ageMap))] <- 0L
  if (maxAgeOnMap > 1e3)
    stop("This module currently assumes that the maximum age is 999. It is not. Please adjust ages or this module")
  digits <- nchar(ceiling(maxAgeOnMap / pctRound))
  roundedAge <- round(round(ageMap[-mapCodesNAs], -1) / pctRound, 0)
  xTileAge <- paddedFloatToChar(roundedAge,
                                padL = digits, padR = 0)
  #ids <- which(as.numeric(mapCodes) == 0)
  mapCodeGrps <- paste0(mapCodesSpp, xTileAge)
  rm(xTileAge)

  ## convert mapCodes for use with data.table and as integer raster
  mapCodesFac <- factor(mapCodeGrps)
  mapCodesInt <- as.integer(mapCodesFac)
  allPixels <- rep(NA, ncell(speciesLayers))
  allPixels[-mapCodesNAs] <- TRUE
  speciesComMap <- raster(speciesLayers[[1]])
  #speciesComMap[] <- NA
  speciesComMap[allPixels] <- mapCodesInt # integer is OK now that it is factor

  message("There are ", NROW(unique(mapCodesFac)), " initial, unique communities, \n  based on ",
          "species abundance (rounded to ", pctRound, "-percentile groups -- i.e., ",
          percentileGrps, " groups) ", " and \n  initial stand age (rounded to ", pctRound, "-year groups)\n",
          "  Modify percentileGrps to change the number of initial communties")


  initialCommunities <- data.table(mapcode = mapCodesSpp,
                                   mapCodeInt = mapCodesInt,
                                   mapCodeFac = mapCodesFac,
                                   speciesPresence = as.vector(t(xtileAbundNum)) * percentileGrps,
                                   species = as.vector(t(sppMatrix)),
                                   age1 = roundedAge * percentileGrps,
                                   pixelIndex = which(allPixels))
  initialCommunities <- initialCommunities[speciesPresence > 0]
  # LIKELY DELETE BELOW THIS
  # if (FALSE) {
  #   mapCodesInt <- as.numeric(mapCodeGrps)
  #   #mapCodesInt[ids] <- NA_integer_
  #   #mapCodes[ids] <- NA_character_
  #   #speciesComMap <- speciesLayers[[1]]
  #   #speciesComMap[] <- mapCodesInt
  #
  #   initialCommunities <- data.table(mapcode = sort(unique(na.omit(mapCodesInt))))
  #   initialCommunities[, mapCodeStr := sort(unique(na.omit(mapCodesSpp)))]
  #   output <- data.table(mapcode = numeric(), speciesPresence = character(),
  #                        species = character(), age1 = numeric())
  #   for (i in 1:nrow(initialCommunities)) {
  #     outputAdd <- data.table(mapcode = initialCommunities$mapcode[i],
  #                             speciesPresence = substring(initialCommunities$mapCodeStr[i],
  #                                                         seq(1, digits * n, 2),
  #                                                         seq(2, digits * n, 2)),
  #                             species = speciesNames[1:length(speciesNames)],
  #                             age1 = as.numeric(rep(substring(initialCommunities$mapCodeStr[i],
  #                                                             2 + digits * n - 1,
  #                                                             2 + digits * n), n)) * 10)
  #
  #     output <- rbind(output, outputAdd)
  #   }
  #
  #   initialCommunities <- output[speciesPresence != "00",]
  #   initialCommunities[, newMapCode := as.numeric(as.factor(mapcode))]
  #   mapcodeconnection <- unique(initialCommunities[, .(mapcode, newMapCode)], by = "mapcode")
  #   indexTable <- data.table(pixelIndex = 1:ncell(speciesComMap),
  #                            mapcode = getValues(speciesComMap))
  #   indexTable <- indexTable[!is.na(mapcode), ]
  #   indexTable <- setkey(indexTable, mapcode)[setkey(mapcodeconnection, mapcode), nomatch = 0]
  #   speciesComMap[indexTable$pixelIndex] <- indexTable$newMapCode
  #   initialCommunities[, ':='(mapcode = newMapCode, newMapCode = NULL, speciesPresence = NULL)]
  # }
  return(list(initialCommunityMap = speciesComMap,
              initialCommunity = initialCommunities))
}
