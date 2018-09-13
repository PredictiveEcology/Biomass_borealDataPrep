obtainSEP <- function(speciesLayers, ecoregionMap, SEPMinThresh = 0) {
  #speciesLayerEcoregion <- crop(speciesLayers, ecoregionMap)
  #speciesLayerEcoregion <- suppressWarnings(mask(speciesLayerEcoregion, ecoregionMap))
  speciesTableEcoregion <- cbind(data.table(ecoregion = getValues(ecoregionMap)),
                                 data.table(getValues(speciesLayers)))
  speciesTableEcoregion <- speciesTableEcoregion[complete.cases(speciesTableEcoregion)]
  speciesTableEcoregion <- melt.data.table(speciesTableEcoregion,
                                           measure.vars = names(speciesTableEcoregion)[-1],
                                           variable.name = "species")
  speciesTableEcoregion[value <= SEPMinThresh, presence := 0]
  speciesTableEcoregion[value > SEPMinThresh, presence := 1]
  speciesTableEcoregion <- speciesTableEcoregion[,.(SEP = sum(presence)/length(value)),
                                                 by = c("ecoregion", "species")]
  
  # This next section is just for visualizing
  if (FALSE) {
    speciesLevels <- unique(speciesTableEcoregion$species)
    abundanceMapStack <- stack()
    names(ecoregionMap) <- "ecoregion"
    for (i in 1:length(speciesLevels)) {
      speciesTableEcoregionBySpecies <- speciesTableEcoregion[species == speciesLevels[i], ][
        , species := NULL]
      abundanceMap <- rasterizeReduced(speciesTableEcoregionBySpecies, ecoregionMap, "SEP")
      names(abundanceMap) <- as.character(speciesLevels[i])
      abundanceMapStack <- stack(abundanceMapStack, abundanceMap)
    }
  }
  return(speciesAbundanceTable = speciesTableEcoregion)
}
