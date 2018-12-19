obtainMaxBandANPP <- function(speciesLayers, biomassLayer, #SALayer,
                              ecoregionMap,
                              minNumPixelsToEstMaxBiomass = 100,
                              quantileForMaxBiomass = 0.99) {
                              #pctCoverMinThresh = 0) {
    # These commented lines are superceded by prepInputs
  # speciesinstudyarea <- crop(speciesLayers, ecoregionMap) # slowest
  # speciesinstudyarea <- suppressWarnings(mask(speciesinstudyarea, ecoregionMap))
  # biomassinstudyarea <- crop(biomassLayer, ecoregionMap)
  # biomassinstudyarea <- suppressWarnings(mask(biomassinstudyarea, ecoregionMap))
  # SAinstudyarea <- crop(SALayer, ecoregionMap)
  # SAinstudyarea <- suppressWarnings(mask(SAinstudyarea, ecoregionMap))

  speciesTable <- data.table(biomass = getValues(biomassLayer),
                             #SA = getValues(SALayer),
                             ecoregion = getValues(ecoregionMap))
  speciesTableWPct <- data.table(speciesTable,
                                 percentage = speciesLayers[],
                                 totalpct = apply(speciesLayers[], 1, sum)) # they don't all add up to 100%, so standardize
  #speciesTableWPct[totalpct == 0, totalpct]
  speciesTableWPct <- speciesTableWPct[!is.na(totalpct) & totalpct != 0]

  speciess <- names(speciesLayers)

  for (species in speciess) {
    pctSp <- paste0("percentage.", species)
    biomassSp <- paste0("biomass", species)

    # rescaling pct cover so all sum to 100
    set(speciesTableWPct, NULL, pctSp, as.integer(100 * speciesTableWPct[[pctSp]]/speciesTableWPct$totalpct)) # integer precision is enough

    # assigne species-specific biomass, based on biomass * pctCover
    set(speciesTableWPct, NULL, biomassSp, speciesTableWPct$biomass * speciesTableWPct[[pctSp]] / 100)
    #set(speciesTableWPct, NULL, pctSp, NULL) # get rid of pct column
  }

  # Calculate max Biomass
  biomassSpp <- paste0("biomass", speciess)
  outputPartial <- speciesTableWPct[, lapply(.SD, function(b) {
    #100 * quantile(b, 0.8, na.rm = TRUE)
    logBgtZero <- log(b[b > 0])
    100 * exp(quantile(logBgtZero, quantileForMaxBiomass)) # take log to change from a skewed distribution
    }),
    by = c("ecoregion"), .SDcols = biomassSpp]

  # This calculates sample sizes by species-ecoregion
  NbyEcoregionSpecies <- speciesTableWPct[, lapply(.SD, function(b) {
    N <- sum(b > 0)
    if (N > minNumPixelsToEstMaxBiomass)
      return(N)
    else
      return (NA_integer_)
    }),
    by = c("ecoregion"), .SDcols = biomassSpp]

  # Remove biomass entries that have NA in the sample size matrix -- i.e., fewer than minNumPixelsToEstMaxBiomass
  biomassWithEnoughN <- as.matrix(outputPartial) * as.matrix(NbyEcoregionSpecies > 0)
  outputPartial <- as.data.table(biomassWithEnoughN)

  # convert biomass to integer as more precision is unnecessary
  outputPartial <- outputPartial[, append(list(ecoregion = ecoregion),
                                          lapply(.SD, function(col) as.integer(round(col)))),
                                 .SDcols = biomassSpp]

  setnames(outputPartial, old = biomassSpp, new = gsub("biomass", "", biomassSpp))

  # Convert to long format
  output <- data.table::melt(outputPartial,
                   value.name = "maxBiomass",
                   measure.vars = names(outputPartial)[-1],
                   variable.name = "species")

  # Convert cases where maxBiomass is 0 to NA
  output <- output[maxBiomass == 0, maxBiomass := NA]
  output[, ecoregionCode := factorValues2(ecoregionMap, output$ecoregion, att = "ecoregion" )]

  # Create maxANPP as maxBiomass / 30
  set(output, NULL, "maxANPP", as.integer(round(output$maxBiomass / 30, 0)))

  return(speciesBiomass = output)
}

