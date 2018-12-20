#' Garbage collection
#'
#' Although R and your operating system should be managing memory for you,
#' the \pkg{sp} and \pkg{raster} packages sometimes suck at memory management.
#'
#' This provides a simple wrapper arounnd \code{\link{gc}} to try to force memory cleanup.
#'
#' @export
.gc <- function() for (i in 1:10) gc()

#' Fasterize polygons using \code{fasterize}
#'
#' @param sp a shapefile to rasterize
#' @param raster the template raster to use
#' @param fieldname the field to use for the rasterizing (will be ignored if the
#'                  shapefile has no fields)
#'
#' @return TODO: is it a \code{RasterLayer}?
#'
#' @export
#' @importFrom fasterize fasterize
#' @importFrom sf st_as_sf
fasterizeFromSp <- function(sp, raster, fieldName) {
  ## check if projections are the same
  if (!identical(crs(sp), crs(raster)))
    stop("fasterize will probably be wrong, as shp and raster projections do not match")

  tempSf <- sf::st_as_sf(sp)

  if (all(names(tempSf) == "geometry")) {
    ## if there are no fields, ignore fieldname
    fasterize::fasterize(tempSf, raster)
  } else
    fasterize::fasterize(tempSf, raster, field = fieldName)
}
