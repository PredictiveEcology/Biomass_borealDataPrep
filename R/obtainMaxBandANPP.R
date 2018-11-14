obtainMaxBandANPP <- function(speciesLayers, biomassLayer, SALayer, ecoregionMap,
                              pctCoverMinThresh = 0) {
  # These commented lines are superceded by prepInputs
  # speciesinstudyarea <- crop(speciesLayers, ecoregionMap) # slowest
  # speciesinstudyarea <- suppressWarnings(mask(speciesinstudyarea, ecoregionMap))
  # biomassinstudyarea <- crop(biomassLayer, ecoregionMap)
  # biomassinstudyarea <- suppressWarnings(mask(biomassinstudyarea, ecoregionMap))
  # SAinstudyarea <- crop(SALayer, ecoregionMap)
  # SAinstudyarea <- suppressWarnings(mask(SAinstudyarea, ecoregionMap))

  speciesTable <- data.table(biomass = getValues(biomassLayer))

  speciesTable[, ':='(SA = getValues(SALayer), ecoregion = getValues(ecoregionMap))]
  outputPartial <- data.table(ecoregion = numeric(), species = character(),
                              maxBiomass = numeric(), maxANPP = numeric())
  speciess <- names(speciesLayers)

  for (species in speciess) {
    indispeciesraster <- raster::subset(speciesLayers, species)
    speciesTable[, percentage := getValues(indispeciesraster)]
    speciesTable_narmed <- speciesTable[!is.na(ecoregion), ]
    speciesTable_narmed[, speciesBiomass := biomass*percentage*0.01] ## TODO: see #9
    speciesTable_narmed <- speciesTable_narmed[percentage > pctCoverMinThresh, ] ## TODO: see #10
    speciesTable_narmed[,species := species]

    ## TODO -- why 100* and why 0.8 ... Ceres and Eliot think that we should use max(speciesBiomass),
    ## as long as there are enough, or perhaps sort(max(speciesBiomass))[50] or something.
    ## quantiles are VERY sensitive to the sample size when distribution is highly skewed, so probably a bad idea.
    speciesTable_narmed <- speciesTable_narmed[, .(maxBiomass = 100 * quantile(speciesBiomass, 0.8, na.rm = TRUE)),
                                               by = c("ecoregion", "species")]
    speciesTable_narmed[, maxANPP := maxBiomass / 30]
    outputPartial <- rbind(outputPartial, speciesTable_narmed)
  }
  output <- data.table(expand.grid(ecoregion = as.numeric(unique(getValues(ecoregionMap))),
                                   species = speciess))[!is.na(ecoregion),][,species := as.character(species)]
  output <- outputPartial[output, on = c("ecoregion", "species"), nomatch = NA]
  #output <- dplyr::left_join(output, outputPartial, by = c("ecoregion", "species")) %>%
  #  data.table()
  return(speciesBiomass = output)
}

