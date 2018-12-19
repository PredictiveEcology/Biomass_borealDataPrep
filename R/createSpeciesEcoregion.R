createSpeciesEcoregion <- function(possibleEcoregionSrcs, rasterToMatch, speciesLayers, biomassMap,
                                   minNumPixelsToEstMaxBiomass, quantileForMaxBiomass) {

  stillHaveNAsForMaxBiomass <- TRUE
  attempt <- 0
  speciesEcoregionTable <- data.table()
  while (stillHaveNAsForMaxBiomass) {
    attempt <- attempt + 1
    if (!is(possibleEcoregionSrcs[[attempt]], "Raster")) {
      # if (is(possibleEcoregionSrcs[[attempt]], "Spatial")) {
      #   possibleEcoregionSrcs[[attempt]] <- fasterizeFromSp(possibleEcoregionSrcs[[attempt]], rasterToMatch,
      #                                                       grep("ECO", names(possibleEcoregionSrcs[[attempt]]), value = TRUE)[1])
      #   setColors(possibleEcoregionSrcs[[attempt]]) <- "Blues"
      #   possibleEcoregionSrcs[[attempt]] <- postProcess(possibleEcoregionSrcs[[attempt]],
      #                                                   rasterToMatch = speciesLayers[[1]],
      #                                                   maskWithRTM = TRUE, filename2 = NULL)
      #   #fasterize::fasterize(possibleEcoregionSrcs[[attempt]],
      #   #                     raster = rasterToMatch)
      # } else {
        stop("All ecoRegion, ecoDistrict, ecoZone maps must be Raster objects")
      #}
    }

    speciesEcoregionTableNew <- Cache(obtainMaxBandANPP,
                                      speciesLayers = speciesLayers,
                                      biomassLayer = biomassMap,
                                      #SALayer = standAgeMap,
                                      #ecoregionMap = ecoregionFiles$ecoregionMap,
                                      ecoregionMap = possibleEcoregionSrcs[[attempt]],
                                      minNumPixelsToEstMaxBiomass = minNumPixelsToEstMaxBiomass,
                                      quantileForMaxBiomass = quantileForMaxBiomass,
                                      #pctCoverMinThresh = 50,
                                      userTags = "stable")
    speciesEcoregionTableNewNAs <- speciesEcoregionTableNew[is.na(maxBiomass),]
    stillHaveNAsForMaxBiomass <- NROW(speciesEcoregionTableNewNAs) > 0

    message(sum(!is.na(speciesEcoregionTableNew$maxBiomass)), " out of ",
            NROW(speciesEcoregionTableNew), " ", names(possibleEcoregionSrcs)[attempt],
            "-species combinations ",
            "had more than ", minNumPixelsToEstMaxBiomass, " data points in: ",
            "the ", names(possibleEcoregionSrcs)[attempt], " raster ",
            #print(possibleEcoregionSrcs[[attempt]])
            "and could therefore estimate maxBiomass and maxANPP")


    speciesEcoregionTable <- if (NROW(speciesEcoregionTable) > 0) { # already exists
      a <- data.table(first = factorValues2(possibleEcoregionSrcs[[1]], possibleEcoregionSrcs[[1]][], att = 5),
                      new = factorValues2(possibleEcoregionSrcs[[attempt]], possibleEcoregionSrcs[[attempt]][], att = 5))
      a <- unique(a)[!is.na(first)]
      colNames <- names(possibleEcoregionSrcs[attempt])
      if ("ecoregionCode" %in% colnames(speciesEcoregionTableNew))
        setnames(speciesEcoregionTableNew, "ecoregionCode", colNames)

      speciesEcoregionTableMerge <- speciesEcoregionTable[a, on = "ecoregionCode==first", nomatch = 0]
      speciesEcoregionTableMerge <- speciesEcoregionTableMerge[speciesEcoregionTableNew, on = c("new==ecoRegion", "species"),
                                                               nomatch = 0]
      speciesEcoregionTableMerge[is.na(maxBiomass), `:=`(maxBiomass = i.maxBiomass, maxANPP = i.maxANPP)]
      speciesEcoregionTable <- speciesEcoregionTableMerge[, `:=`(i.maxBiomass = NULL, i.maxANPP = NULL)]
      # speciesEcoregionTable <- rbindlist(list(speciesEcoregionTable[!is.na(maxBiomass)],
      #                                         speciesEcoregionTableMerge), fill = TRUE)
      setnames(speciesEcoregionTable, "new", "ecoRegion")
      if (isTRUE(stillHaveNAsForMaxBiomass)) {
        message("--------------------------------")
        message("Out of ", length(levels(speciesEcoregionTable$ecoregionCode)), " possible combinations per species, ",
                "there are still species-ecoregion combinations not estimable:")
        print(speciesEcoregionTable[is.na(maxBiomass), list(numEcoregionsWNoMaxBiomass = .N), by = species])
        print(speciesEcoregionTable[is.na(maxBiomass), list(totalMissingMaxBiomass = .N)])
      }
      stillHaveNAsForMaxBiomass <- FALSE
      speciesEcoregionTable
    } else {
      speciesEcoregionTableNew
    }
  }
  #speciesEcoregionTable[ , ecoregion := as.factor(ecoregion)]
  return(speciesEcoregionTable)
}
