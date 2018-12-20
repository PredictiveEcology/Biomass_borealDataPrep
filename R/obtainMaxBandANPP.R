#' Obtain maximum biomass and ANPP from bigger ecological area
#'
#' @param speciesLayers TODO: description needed
#' @param biomassLayer TODO: description needed
#' @param ecoregionMap TODO: description needed
#' @param minNumPixelsToEstMaxBiomass TODO: description needed
#' @param quantileForMaxBiomass TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table as.data.table data.table melt set setnames
#' @importFrom pemisc factorValues2
#' @importFrom raster getValues
obtainMaxBandANPP <- function(speciesLayers, biomassLayer, #SALayer,
                              ecoregionMap,
                              minNumPixelsToEstMaxBiomass = 100,
                              quantileForMaxBiomass = 0.99) {
                              #pctCoverMinThresh = 0) {

  speciesTable <- data.table(biomass = getValues(biomassLayer),
                             #SA = getValues(SALayer),
                             ecoregionCode = factorValues2(ecoregionMap,
                                                           getValues(ecoregionMap), att = 5))
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
    by = c("ecoregionCode"), .SDcols = biomassSpp]

  # This calculates sample sizes by species-ecoregion
  NbyEcoregionSpecies <- speciesTableWPct[, lapply(.SD, function(b) {
    N <- sum(b > 0)
    if (N > minNumPixelsToEstMaxBiomass)
      return(N)
    else
      return (NA_integer_)
    }),
    by = c("ecoregionCode"), .SDcols = biomassSpp]

  ## Remove biomass entries that have NA in the sample size matrix
  ## i.e., fewer than minNumPixelsToEstMaxBiomass
  biomassWithEnoughN <- matrix(nrow = nrow(outputPartial), ncol = ncol(outputPartial)-1)
  biomassWithEnoughN <- as.matrix(outputPartial[,-1]) * as.matrix(NbyEcoregionSpecies[,-1] > 0)
  outputPartial <- data.table(ecoregionCode = outputPartial$ecoregionCode,
                              as.data.table(biomassWithEnoughN))

  # convert biomass to integer as more precision is unnecessary
  outputPartial[, (biomassSpp) := lapply(.SD, function(col) as.integer(round(col))),
                .SDcols = biomassSpp]

  setnames(outputPartial, old = biomassSpp, new = gsub("biomass", "", biomassSpp))

  # Convert to long format
  output <- data.table::melt(outputPartial,
                             value.name = "maxBiomass",
                             measure.vars = names(outputPartial)[-1],
                             variable.name = "species")

  # Convert cases where maxBiomass is 0 to NA
  output <- output[maxBiomass == 0, maxBiomass := NA]

  # Create maxANPP as maxBiomass / 30
  set(output, NULL, "maxANPP", as.integer(round(output$maxBiomass / 30, 0)))

  return(output)
}
