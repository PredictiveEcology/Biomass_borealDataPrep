obtainSEP <- function(speciesLayers, possibleEcoregionSrcs, speciesEcoregionTable) {

  ##################################################################
  # Initiate a data.table with all values
  ##################################################################
  attempt <- 0

  stillHave0sForSEP <- TRUE
  ind <- 0
  while (stillHave0sForSEP) {
    ind <- ind + 1
    ecoregionName <- names(possibleEcoregionSrcs)[ind]
    ecoregionMap <- possibleEcoregionSrcs[[ind]]
    ecoregionCodes <- factorValues2(ecoregionMap,
                                    getValues(ecoregionMap),
                                    att = "ecoregion")
    speciesEcoregionTemp <- data.table(ecoregionCode = ecoregionCodes,
                                       getValues(speciesLayers))
    speciesEcoregionTemp[, ecoregionInt := as.integer(ecoregionCode)]
    speciesEcoregionTemp <- speciesEcoregionTemp[complete.cases(speciesEcoregionTemp)]

    speciesEcoregionTemp <- melt.data.table(speciesEcoregionTemp,
                                            value.name = "SEP",
                                            measure.vars = names(speciesLayers),
                                            variable.name = "species")
    speciesEcoregionTemp[SEP == 0, presence := 0]
    speciesEcoregionTemp[SEP > 0, presence := 1]
    speciesEcoregionTemp <- speciesEcoregionTemp[,.(SEP = sum(presence)/length(SEP),
                                                    ecoregionInt = unique(ecoregionInt)),
                                                 by = c("ecoregionCode", "species")]

    # remove the cases where all species in a ecoregion are zero
    speciesEcoregionTemp <- speciesEcoregionTemp[, list(species, SEP, ecoregionInt, hasAllZero = all(SEP == 0)),
                                                 by = "ecoregionCode"]
    speciesEcoregionTemp <- speciesEcoregionTemp[hasAllZero == FALSE]
    speciesEcoregionTemp[, hasAllZero := NULL]
    speciesEcoregionTemp[, SEP := round(SEP, 4)]

    # merge species names as factors
    ecoregionLabels <- raster::levels(ecoregionMap)[[1]]$ecoregion
    speciesEcoregionTemp[, `:=`(ecoregionCode = factor(speciesEcoregionTemp$ecoregionCode,
                                                       labels = ecoregionLabels,
                                                       levels = ecoregionLabels))]

    if (ind == 1) {
      speciesEcoregionNoZeros <- speciesEcoregionTemp[SEP > 0]
    } else {
      stillHave0sForSEP <- !isTRUE(ind == length(possibleEcoregionSrcs))
      noSEP <- speciesEcoregionTable[!speciesEcoregionNoZeros, on = c("ecoregionCode", "species")]
      newSEP <- noSEP[speciesEcoregionTemp, on = c(paste0(ecoregionName, "==ecoregionCode"), "species"), nomatch = 0]
      speciesEcoregionTemp <- rbindlist(list(speciesEcoregionNoZeros, newSEP[, names(speciesEcoregionNoZeros), with = FALSE]))
    }

  }
  # Fill in 0 for maxBiomass and maxANPP when SEP was estimated to be 0
  #speciesEcoregionTemp[SEP == 0, ':='(maxBiomass = 0, maxANPP = 0)]
  speciesEcoregionTemp[, ecoregionInt := NULL]
  return(speciesEcoregionTemp)
}


createSEPStack <- function(speciesEcoregionTable, ecoregionMap, destinationPath) {
  ##################################################################
  # This next section creates the stack
  ##################################################################
  speciesLevels <- unique(speciesEcoregionTable$species)
  speciesEcoregionTable[, ecoregionInt := as.integer(ecoregionCode)]
  abundanceMapStack <- stack()
  names(ecoregionMap) <- "ecoregion"
  names(speciesLevels) <- as.character(speciesLevels)
  abundanceList <- lapply(speciesLevels, function(sp) {
  #for (i in seq(speciesLevels)) {
    speciesEcoregionBySpecies <- speciesEcoregionTable[species == sp, ][
      , species := NULL]
    abundanceMap <- rasterizeReduced(speciesEcoregionBySpecies, ecoregionMap, "SEP",
                                     mapcode = "ecoregionInt")
    names(abundanceMap) <- as.character(sp)

    abundanceMap <- writeRaster(abundanceMap,
                                filename =
                                  file.path(destinationPath,
                                            paste0("SpeciesEstabProb_", sp, ".tif")),
                                overwrite = TRUE)
  })

  return(stack(abundanceList))

}
