#' Is `studyArea` inside `studyArea_biomassParam`?
#'
#' `sf::st_within()` returns one list element per feature of `x`, so `length()` of it is the
#' number of features (always 1 here) and the old check could never fire; `lengths()` gives the
#' number of `y` features that contain each `x` feature. `y` is merged first, so a `studyArea`
#' that spans two polygons of `studyArea_biomassParam` still counts as inside.
#'
#' @param studyArea `sf` with one feature.
#' @param studyArea_biomassParam `sf` with one or more features; projected to `studyArea`'s CRS
#'   if it differs.
#'
#' @return `TRUE` if `studyArea` lies within `studyArea_biomassParam`.
#' @noRd
studyAreaIsWithin <- function(studyArea, studyArea_biomassParam) {
  if (sf::st_crs(studyArea_biomassParam) != sf::st_crs(studyArea)) {
    studyArea_biomassParam <- sf::st_transform(studyArea_biomassParam, sf::st_crs(studyArea))
  }
  isTRUE(lengths(sf::st_within(studyArea, sf::st_union(studyArea_biomassParam)))[[1]] > 0)
}
