#' @importFrom sf st_cast st_transform
#' @importFrom fasterize fasterize
#' @importFrom raster crs
#' @importFrom LandR asInteger
#' @importFrom reproducible Cache prepInputs

prepInputsStandAgeMap <- function(..., ageURL, ageFun, maskWithRTM,
                                  method, datatype, filename2,
                                  fireURL, fireFun,
                                  rasterToMatch, fireField, startTime) {
  standAgeMap <- Cache(prepInputs, ...,
                       maskWithRTM = maskWithRTM, method = method,
                       datatype = datatype, filename2 = filename2,
                       url = ageURL, fun = ageFun, rasterToMatch = rasterToMatch)
  standAgeMap[] <- asInteger(standAgeMap[])
  if (!(missing(fireURL) || is.null(fireURL) || is.na(fireURL))) {
    fireYear <- Cache(prepInputsFireYear, ...,
                      url = fireURL,
                      fun = fireFun,
                      rasterToMatch = rasterToMatch,
                      field = fireField
    )
    toChange <- !is.na(fireYear[]) & fireYear[] <= asInteger(startTime)
    standAgeMap[] <- asInteger(standAgeMap[])
    standAgeMap[toChange] <- asInteger(startTime) - asInteger(fireYear[][toChange])
  }
  standAgeMap
}

prepInputsFireYear <- function(..., rasterToMatch, field) {
  a <- Cache(prepInputs, ...)
  gg <- st_cast(a, "MULTIPOLYGON") # collapse them into a single multipolygon
  d <- st_transform(gg, crs(rasterToMatch))
  fasterize(d, raster = rasterToMatch, field = field)
}
