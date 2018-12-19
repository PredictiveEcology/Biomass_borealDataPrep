initialCommunityProducer <- function(speciesLayers, #speciesPresence,
                                     ecoregionMap,
                                     #rstStudyArea,
                                     standAgeMap,
                                     percentileGrps = 10) {

  speciesNames <- names(speciesLayers)

  n <- length(speciesNames)
  pctRound <- 100 / percentileGrps
  digits <- nchar(percentileGrps) # use quintiles, i.e., 100 pct / 20 ==> 1, 2, 3, 4, 5
  mapCodes <- rep("", ncell(speciesLayers[[1]]))
  mapCodesNAs <- which(is.na(speciesLayers[[1]][]))
  mapCodes <- mapCodes[-mapCodesNAs]

  ############################################################
  # Get a full length logical vector of the non-NA cells
  ############################################################
  allPixels <- rep(NA, ncell(speciesLayers))
  allPixels[-mapCodesNAs] <- TRUE

  ############################################################
  # matrix of all pct covers, collapse into single string
  ############################################################
  xtileAbundNum <- round(round(speciesLayers[][-mapCodesNAs,], -1) / pctRound, 0) # use
  xtileAbund <- if (digits > 1)
    matrix(paddedFloatToChar(xtileAbundNum, digits, 0), ncol = ncol(xtileAbundNum))
  else
    xtileAbundNum
  mapCodesSpp <- apply(xtileAbund, 1, paste, collapse = "")

  # Note that a dataset may have a pixel with zero cover, but
  #   also have ages that are non-zero. This is likely
  #   a bug in the data. Check here and provide message
  mapCodesZeros <- mapCodesSpp == "0000000000"
  pixelsWithZeroCover <- which(allPixels)[mapCodesZeros]

  ###################
  # Species matrix, simply names. This is used later
  ##########
  sppMatrix <- matrix(rep(speciesNames, each = NROW(xtileAbundNum)), ncol = ncol(xtileAbundNum))
  sppMatrix[xtileAbundNum==0] <- ""

  ############################################################
  ## stand age
  ############################################################
  ageMap <- standAgeMap[]
  ageMapOnlyNAs <- which(is.na(ageMap[-mapCodesNAs]) & mapCodesZeros)
  if (length(ageMapOnlyNAs) > 0)
    message("##############################\n",
            "There are ", length(ageMapOnlyNAs), " pixels on the standAgeMap ",
            "that are NA, yet have non-zero percent cover ",
            "on the speciesLayers stack. Proceeding anyway.")

  maxAgeOnMap <- max(ageMap[-mapCodesNAs], na.rm = TRUE)
  # ageMap[which(is.na(ageMap))] <- 0L
  if (maxAgeOnMap > 1e3)
    stop("This module currently assumes that the maximum age is 999. It is not. Please adjust ages or this module")
  digits <- nchar(ceiling(maxAgeOnMap / pctRound))
  roundedAge <- round(round(ageMap[-mapCodesNAs], -1) / pctRound, 0)

  pixelsWithNonZeroAge <- which(allPixels)[roundedAge > 0]
  pixelsWithZeroCoverAndAgeGT0 <- pixelsWithNonZeroAge[pixelsWithNonZeroAge %in% pixelsWithZeroCover]

  if (length(pixelsWithZeroCoverAndAgeGT0) > 0) {
    message("There are ", length(pixelsWithZeroCoverAndAgeGT0),
            " pixels with zero tree cover, but non-zero age. This is ",
            "possibly a bug in the speciesLayers raster stack. ",
            "Proceeding anyway, assuming these pixels have no trees")
    # set ages to 0 for these pixels
    whNonZeroAgeShortVec <- which(allPixels) %in% pixelsWithZeroCoverAndAgeGT0
    roundedAge[whNonZeroAgeShortVec] <- 0
  }

  xTileAge <- paddedFloatToChar(roundedAge,
                                padL = digits, padR = 0)
  xTileAge[ageMapOnlyNAs] <- "NA"

  ############################################################
  ## ecoregion
  ############################################################
  ecoregion <- ecoregionMap[]
  ecoregionOnlyNAs <- which(is.na(ecoregion[-mapCodesNAs]) & mapCodesZeros)
  if (length(ecoregionOnlyNAs) > 0)
    message("##############################\n",
            "There are ", length(ecoregionOnlyNAs), " pixels on the ecoregionMap ",
            "that are NA, yet have non-zero percent cover ",
            "on the speciesLayers stack. Proceeding anyway.")


  if (raster::is.factor(ecoregionMap)) {
    ecoregionValues <- factorValues2(ecoregionMap, ecoregionMap[][-mapCodesNAs], att = "ecoregion")
    uniqueEcoregionCodes <- unique(ecoregion, na.rm = TRUE)
    ecoregionValuesChar <- as.character(ecoregionValues)
    xTileEcoregion <- ecoregionValuesChar
  } else {
    ecoregionValues <- ecoregion[][-mapCodesNAs]
    digits <- max(nchar(ecoregionValues))
    xTileEcoregion <- paddedFloatToChar(ecoregionValues, padL = digits, padR = 0)
    uniqueEcoregionCodes <- unique(ecoregionValues, na.rm = TRUE)
  }

  xTileEcoregion[ecoregionOnlyNAs] <- "NA"

  ############################################################
  ##### Make unique map code
  ############################################################
  mapCodeGrps <- paste0(mapCodesSpp, "_", xTileAge, "_", xTileEcoregion)
  rm(xTileAge)

  ## convert mapCodes for use with data.table and as integer raster
  mapCodesFac <- factor(mapCodeGrps)
  mapCodesInt <- as.integer(mapCodesFac)

  message("There are ", NROW(unique(mapCodesFac)), " initial, unique communities, \n  based on ",
          "species abundance (rounded to ", pctRound, "-percentile groups -- i.e., ",
          percentileGrps, " groups) ", "\n  initial stand age (rounded to ", pctRound, "-year groups),\n",
          "  and ecoregionMap.",
          "  Modify percentileGrps to change the number of initial communties")

  ############################################################
  # Create speciesComMap
  ############################################################
  speciesComMap <- raster(speciesLayers[[1]])
  speciesComMap[allPixels] <- mapCodesInt # integer is OK now that it is factor


  ############################################################
  # Create initialCommunities object
  ############################################################
  initialCommunitiesWide <- data.table(mapcode = mapCodesFac, #mapCodesSpp,
                                   mapCodeInt = mapCodesInt,
                                   speciesPresence = xtileAbundNum * percentileGrps,
                                   species = sppMatrix,
                                   age1 = roundedAge * percentileGrps,
                                   pixelIndex = which(allPixels))

  speciesNames1 <- grep("species\\.", colnames(initialCommunitiesWide), value = TRUE)
  speciesPresence1 <- grep("speciesPresence", colnames(initialCommunitiesWide), value = TRUE)
  initialCommunities <- data.table::melt(
    initialCommunitiesWide,
    measure.vars = list("species" = speciesNames1,
                        "speciesPresence" = speciesPresence1))

  # Remove any lines where there is no cover
  initialCommunities <- initialCommunities[speciesPresence > 0]

  return(list(initialCommunityMap = speciesComMap,
              initialCommunity = initialCommunities))
}
