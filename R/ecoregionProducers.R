ecoregionProducer <- function(ecoregionMap, ecoregionName,
                              ecoregionActiveStatus, rasterToMatch) {
  # change the coordinate reference for all spatialpolygons
  message("ecoregionProducer 1: ", Sys.time())
  #ecoregionMapInStudy <- raster::intersect(ecoregionMapFull, fixErrors(aggregate(studyArea)))

  # Alternative
  message("ecoregionProducer fastRasterize: ", Sys.time())
  ecoregionMap <- fasterize::fasterize(sf::st_as_sf(ecoregionMap), raster(rasterToMatch), field = "ECODISTRIC")
  ecoregionMap[is.na(rasterToMatch[])] <- NA

  ecoregionFactorValues <- na.omit(unique(ecoregionMap[]))

  ecoregionTable <- data.table(
    mapcode = seq_along(ecoregionFactorValues[!is.na(ecoregionFactorValues)]),
    ecoregion = as.numeric(ecoregionFactorValues[!is.na(ecoregionFactorValues)])
  )
  message("ecoregionProducer mapvalues: ", Sys.time())
  ecoregionMap[] <- plyr::mapvalues(ecoregionMap[], from = ecoregionTable$ecoregion, to = ecoregionTable$mapcode)
  ecoregionActiveStatus[, ecoregion := as.character(ecoregion)]
  ecoregionTable <- ecoregionTable[!is.na(mapcode),][, ecoregion := as.character(ecoregion)]
  message("ecoregionProducer dplyr_leftjoin: ", Sys.time())
  ecoregionTable <- dplyr::left_join(ecoregionTable,
                                     ecoregionActiveStatus,
                                     by = "ecoregion") %>%
    data.table()
  ecoregionTable[is.na(active), active := "no"]
  ecoregionTable <- ecoregionTable[,.(active, mapcode, ecoregion)]
  return(list(ecoregionMap = ecoregionMap,
              ecoregion = ecoregionTable))
}

# nonActiveEcoregionProducer <- function(nonactiveRaster, activeStatus, ecoregionMap,
#                                        ecoregion, initialCommunityMap, initialCommunity) {
#   nonactiveRasterSmall <- crop(nonactiveRaster, ecoregionMap)
#   nonecomapcode <- activeStatus[active == "no", ]$mapcode
#   whNANonActiveRasterSmall <- which(nonactiveRaster[] %in% nonecomapcode)
#   initialCommunityMap[whNANonActiveRasterSmall] <- NA
#   ecoregionMap[whNANonActiveRasterSmall] <- NA
#
#   initialCommunity <- initialCommunity[mapCodeFac %in% sort(unique(getValues(initialCommunityMap))), ]
#   ecoregion <- ecoregion[mapcode %in% sort(unique(getValues(ecoregionMap))), ]
#   return(list(ecoregionMap = ecoregionMap,
#               ecoregion = ecoregion,
#               initialCommunityMap = initialCommunityMap,
#               initialCommunity = initialCommunity))
# }
