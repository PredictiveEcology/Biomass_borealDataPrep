createInitCommMap <- function(initCommMap, values, filename) {
  map <- raster::setValues(initCommMap, values = values)
  raster::writeRaster(map, overwrite = TRUE, filename = filename, datatype = "INT2U")
}

toSentenceCase <- function(x) {
  newNames <- tolower(x)
  substr(newNames, 1, 1) <- toupper(substr(newNames, 1, 1))
  newNames
}

## ------------------------------------------------------------------------
## FASTERIZE POLYGONS USING FASTERIZE

## sp: a shapefile to rasterize
## raster: the template raster to use
## fieldname: the field to use for the rasterizing.
##  Will be ignored if the shapefile has no fields

fasterizeFromSp <- function(sp, raster, fieldName) {
  ## check if projections are the same
  if(!identical(crs(sp), crs(raster))) 
    stop("fasterize will probably be wrong, as shp and raster projections do not match")
  
  tempSf <- sf::st_as_sf(sp)
  
  if(all(names(tempSf) == "geometry")) {
    ## if there are no fields, ignore fieldname
    fasterize::fasterize(tempSf, raster)
  } else 
    fasterize::fasterize(tempSf, raster, field = fieldName)
}