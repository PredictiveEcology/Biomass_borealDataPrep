#' Produce initial communities
#'
#' TODO: description needed
#'
#' @param speciesLayers TODO: description needed
#' @param ecoregionMap  TODO: description needed
#' @param standAgeMap  TODO: description needed
#' @param percentileGrps  TODO: description needed
#'
#' @return TODO: description needed
#'
#' @export
#' @importFrom data.table melt set
#' @importFrom pemisc factorValues2
#' @importFrom raster is.factor
#' @importFrom SpaDES.core paddedFloatToChar
initialCommunityProducer <- function(speciesLayers, #speciesPresence,
                                     ecoregionMap,
                                     #rstStudyArea,
                                     standAgeMap,
                                     percentileGrps = 10) {

  speciesNames <- names(speciesLayers)
  names(speciesNames) <- speciesNames

  n <- length(speciesNames)
  pctRound <- 100 / percentileGrps

  pixels <- seq(ncell(speciesLayers))
  initialCommunities <- lapply(speciesNames, function(sp) {
    dt <- na.omit(data.table(speciesCode = sp,
                             cover = round(speciesLayers[[sp]][]/ pctRound, 0) * pctRound,
                             pixelIndex = pixels))
  })
  initialCommunities <- rbindlist(initialCommunities)
  initialCommunities[, `:=`(age = asInteger(round(standAgeMap[][pixelIndex] / pctRound, 0) * pctRound),
                            ecoregionGroup = factorValues2(ecoregionMap,
                                                           ecoregionMap[][pixelIndex], att = 5))]

  # Case of pixels with species cover, but they NA for age
  ageMapOnlyNAs <- initialCommunities[is.na(age), unique(pixelIndex)]
  if (length(ageMapOnlyNAs) > 0) {
    message("##############################\n",
            "There are ", length(ageMapOnlyNAs), " pixels (",
            round(length(ageMapOnlyNAs)/length(unique(initialCommunities$pixelIndex)), 3),
            " proportion of pixels with data) on the standAgeMap ",
            "that are NA, yet have non-zero percent cover ",
            "on the speciesLayers stack. Setting their age to the mean of their 8 neighbours")
    b <- adj(standAgeMap, cells = ageMapOnlyNAs, pairs = TRUE, include = FALSE, returnDT = TRUE)
    setnames(b, old = "from", new = "pixelIndex")
    dd <- b[, list(newage = asInteger(round(mean(standAgeMap[][to], na.rm = TRUE) / pctRound, 0) * pctRound)),
            by = "pixelIndex"]
    initialCommunities <- dd[initialCommunities, on = "pixelIndex"]
    initialCommunities[is.na(age), age := newage]
    initialCommunities[, newage := NULL]
  }

  # Case of zero species cover, but they have age
  pixelsWithZeroCoverAndAgeGT0 <- initialCommunities[age > 0, all(cover == 0), by = "pixelIndex"][V1 == TRUE]
  pixelsWithZeroCoverAndAgeGT0 <- unique(pixelsWithZeroCoverAndAgeGT0, by = "pixelIndex")
  pixelsWithZeroCoverAndAgeGT0[, `:=`(V1 = NULL, newage = 0)]
  if (NROW(pixelsWithZeroCoverAndAgeGT0) > 0) {
    message("There are ", NROW(pixelsWithZeroCoverAndAgeGT0),
            " pixels with zero tree cover, but non-zero age.",
            " This is possibly a bug in the speciesLayers raster stack.",
            " Proceeding anyway, assuming these pixels have no trees.")
  }
  # set ages to 0 for these pixels
  initialCommunities <- pixelsWithZeroCoverAndAgeGT0[initialCommunities, on = "pixelIndex"]

  initialCommunities <- initialCommunities[cover > 0 | newage == 0]
  initialCommunities[newage == 0, `:=`(age = 0, speciesCode = "empty")]
  initialCommunities <- unique(initialCommunities,
                               by = c("pixelIndex", "speciesCode", "cover", "age", "ecoregionGroup"))
  initialCommunities[, newage := NULL]

  # Check standAge is OK
  mapCodesNAs <- which(is.na(speciesLayers[[1]][]))
  maxAgeOnMap <- max(standAgeMap[][-mapCodesNAs], na.rm = TRUE)
  # ageMap[which(is.na(ageMap))] <- 0L
  if (maxAgeOnMap > 1e3)
    stop("This module currently assumes that the maximum age is 999.",
         " It is not. Please adjust ages or this module.")

  message("running landR::addPixelGroup to create initial communities")
  initialCommunities <- addPixelGroup(initialCommunities, 0,
                                      columns = c("ecoregionGroup", "speciesCode", "age", "cover"),
                                      successionTimestep = 10)

  setnames(initialCommunities, old = "pixelGroup", new= "mapcode")
  message("There are ", NROW(unique(initialCommunities$mapcode)), " initial, unique communities",
          " (in ", length(unique(initialCommunities$pixelIndex))," pixels),",
          "\n  based on ",
          "species abundance (rounded to ", pctRound, "-percentile groups -- i.e., ",
          percentileGrps, " groups) ", "\n  initial stand age (rounded to ", pctRound, "-year groups),\n",
          "  and ecoregionMap.",
          "  Modify percentileGrps to change the number of initial communties")

  # remove unneeded columns of initialCommunities based on sim$species
  initialCommunities[, `:=`(speciesPresence = as.integer(cover),
                            description = NA)]
  initialCommunities[, cover2 := NULL]


  return(initialCommunities)
}
