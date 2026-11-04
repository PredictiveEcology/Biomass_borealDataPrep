## The check that `studyArea` lies inside `studyArea_biomassParam` used
## `length(st_within(...))[[1]] == 0`. `length()` of the result is the number of `studyArea`
## features (1), so the check never fired.

squarePoly <- function(xmin, xmax, ymin = 0, ymax = 1000) {
  sf::st_sf(geometry = sf::st_sfc(sf::st_polygon(list(rbind(
    c(xmin, ymin), c(xmax, ymin), c(xmax, ymax), c(xmin, ymax), c(xmin, ymin)
  ))), crs = 3978))
}

test_that("studyAreaIsWithin is TRUE when studyArea is inside, FALSE when it sticks out", {
  sa_bp <- squarePoly(0, 1000)
  expect_true(studyAreaIsWithin(squarePoly(100, 900, 100, 900), sa_bp))
  expect_true(studyAreaIsWithin(sa_bp, sa_bp)) ## the default: studyArea_biomassParam <- studyArea
  expect_false(studyAreaIsWithin(squarePoly(500, 1500), sa_bp))
  expect_false(studyAreaIsWithin(squarePoly(2000, 3000), sa_bp))
})

test_that("a studyArea spanning two studyArea_biomassParam polygons counts as inside", {
  sa_bp <- rbind(squarePoly(0, 1000), squarePoly(1000, 2000))
  expect_true(studyAreaIsWithin(squarePoly(500, 1500, 100, 900), sa_bp))
})

test_that("studyAreaIsWithin projects studyArea_biomassParam to studyArea's CRS", {
  sa_bp <- squarePoly(0, 1000)
  sa_bp_ll <- sf::st_transform(sa_bp, 4326)
  expect_true(studyAreaIsWithin(squarePoly(100, 900, 100, 900), sa_bp_ll))
  expect_false(studyAreaIsWithin(squarePoly(500, 1500), sa_bp_ll))
})
