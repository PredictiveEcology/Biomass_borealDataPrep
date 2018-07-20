createInitCommMap <- function(initCommMap, values, filename) {
  map <- raster::setValues(initCommMap, values = values)
  raster::writeRaster(map, overwrite = TRUE, filename = filename, datatype = "INT2U")
}

toSentenceCase <- function(x) {
  newNames <- tolower(x)
  substr(newNames, 1, 1) <- toupper(substr(newNames, 1, 1))
  newNames
}
