obtainSEP <- function(speciesLayers, ecoregionMap, speciesEcoregionTable,
                      destinationPath) {

  ##################################################################
  # Initiate a data.table with all values
  ##################################################################
  speciesTableEcoregion <- data.table(ecoregion = factorValues2(ecoregionMap,
                                                                getValues(ecoregionMap),
                                                                att = "ecoregion"),
                                      getValues(speciesLayers))
  speciesTableEcoregion[, ecoregionInt := as.integer(ecoregion)]
  speciesTableEcoregion <- speciesTableEcoregion[complete.cases(speciesTableEcoregion)]

  speciesTableEcoregion <- melt.data.table(speciesTableEcoregion,
                                           value.name = "SEP",
                                           measure.vars = names(speciesLayers),
                                           variable.name = "species")
  speciesTableEcoregion[SEP == 0, presence := 0]
  speciesTableEcoregion[SEP > 0, presence := 1]
  speciesTableEcoregion <- speciesTableEcoregion[,.(SEP = sum(presence)/length(SEP),
                                                    ecoregionInt = unique(ecoregionInt)),
                                                 by = c("ecoregion", "species")]
  speciesTableEcoregion[, SEP := round(SEP, 4)]

  # merge species names as factors
  ecoregionLabels <- raster::levels(ecoregionMap)[[1]]$ecoregion
  speciesTableEcoregion[, `:=`(ecoregion = factor(speciesTableEcoregion$ecoregion,
                       labels = ecoregionLabels,
                       levels = ecoregionLabels))]

  speciesTableEcoregion <- speciesEcoregionTable[speciesTableEcoregion, on = c(ecoregionCode = "ecoregion", "species")]
  # Fill in 0 for maxBiomass and maxANPP when SEP was estimated to be 0
  speciesTableEcoregion[SEP == 0, ':='(maxBiomass = 0, maxANPP = 0)]

  ##################################################################
  # This next section creates the stack
  ##################################################################
  speciesLevels <- unique(speciesTableEcoregion$species)
  abundanceMapStack <- stack()
  names(ecoregionMap) <- "ecoregion"
  for (i in 1:length(speciesLevels)) {
    speciesTableEcoregionBySpecies <- speciesTableEcoregion[species == speciesLevels[i], ][
      , species := NULL]
    abundanceMap <- rasterizeReduced(speciesTableEcoregionBySpecies, ecoregionMap, "SEP",
                                     mapcode = "ecoregionInt")
    names(abundanceMap) <- as.character(speciesLevels[i])

    abundanceMap <- writeRaster(abundanceMap,
                                filename =
                                  file.path(destinationPath,
                                            paste0("SpeciesEstabProb_", speciesLevels[i], ".tif")),
                                overwrite = TRUE)
    abundanceMapStack <- stack(abundanceMapStack, abundanceMap)

  }

  return(list(SEPTable = speciesTableEcoregion,
         SEPMap = abundanceMapStack))
}
