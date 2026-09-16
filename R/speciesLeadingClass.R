## Land-cover class of a pixel from its species cover, for the class-240 fill.
##
## 210, 220 and 230 are statements about composition -- a pixel cannot be pine-dominated and
## carry a broadleaf class -- so species cover can decide them. The reverse is not true: a map
## that says "broadleaf" cannot say which species. That is why the fill takes composition from
## here rather than from another year's land cover or from a neighbouring pixel's class.
##
## Site is NOT decided here. Whether a pixel is wet (80/81) is a property of the ground, and
## no amount of species cover expresses it; the caller keeps site from the land-cover record
## and uses this only to choose the composition class within an upland pixel.
##
## Thresholds match `LandR::vegTypeMapGenerator(vegLeadingProportion = 0.8, mixedType = 2)`,
## measured directly: conifer share of tree cover >= 0.8 is conifer-leading, <= 0.2 is
## deciduous-leading, and anything between is mixed.
##
## `coverDT` has one row per pixel and one numeric column per species, in percent. Column names
## are species codes in the `sppEquivCol` convention, optionally prefixed "cover." as
## `makePixelTable()` leaves them. Species absent from `sppEquiv` are ignored.
##
## Returns an integer vector, one per row: 210, 220, 230, or NA where the pixel has no tree
## cover at all and therefore nothing to infer from -- those pixels are left for the caller's
## fallback rather than being forced into a class.
speciesLeadingClass <- function(coverDT, sppEquiv, sppEquivCol = "LandR",
                                vegLeadingProportion = 0.8) {
  stopifnot(vegLeadingProportion > 0.5, vegLeadingProportion <= 1)

  sppCols <- grep("^cover\\.", colnames(coverDT), value = TRUE)
  bare <- if (length(sppCols)) sub("^cover\\.", "", sppCols) else colnames(coverDT)
  if (!length(sppCols)) {
    sppCols <- bare[bare %in% sppEquiv[[sppEquivCol]]]
    bare <- sppCols
  }
  keep <- bare %in% sppEquiv[[sppEquivCol]]
  sppCols <- sppCols[keep]
  bare <- bare[keep]
  if (!length(sppCols)) {
    return(rep(NA_integer_, NROW(coverDT)))
  }

  ## Conifer or deciduous, from the species table rather than a hard-coded list.
  type <- sppEquiv[["Type"]][match(bare, sppEquiv[[sppEquivCol]])]
  isConifer <- !is.na(type) & type == "Conifer"
  isDecid <- !is.na(type) & type == "Deciduous"

  m <- as.matrix(coverDT[, sppCols, with = FALSE])
  m[is.na(m)] <- 0
  conCover <- if (any(isConifer)) rowSums(m[, isConifer, drop = FALSE]) else rep(0, nrow(m))
  decCover <- if (any(isDecid)) rowSums(m[, isDecid, drop = FALSE]) else rep(0, nrow(m))
  total <- conCover + decCover

  out <- rep(NA_integer_, nrow(m))
  has <- total > 0
  frac <- conCover[has] / total[has]

  ## Both thresholds are inclusive, matching `vegTypeMapGenerator()`: measured against it, a
  ## conifer share of exactly 0.80 is conifer-leading and exactly 0.20 is deciduous-leading,
  ## with mixed strictly between. Compare with a tolerance rather than exactly, because
  ## `1 - 0.8` is 0.19999999999999996 while `20/100` is 0.20000000000000001, so an exact
  ## `<=` sent a true 20/80 split to mixed instead of deciduous.
  eps <- sqrt(.Machine$double.eps)
  out[has] <- data.table::fifelse(
    frac >= vegLeadingProportion - eps, 210L,
    data.table::fifelse(frac <= (1 - vegLeadingProportion) + eps, 220L, 230L)
  )
  out
}
