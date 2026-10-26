## Cohort B is rounded half up twice on its way from the raw biomass map (to whole g/m2, then to
## pixelGroupBiomassClass), so the check that every cohort is within the map rounds maxRawB the
## same way. It used round(maxRawB, -2), which stopped NorthernAlbertaUplands: 329.4970703 t/ha
## gave a cohort of 33000 against a limit of 32900.

moduleBody <- function(name) {
  exprs <- parse(testthat::test_path("..", "..", "Biomass_borealDataPrep.R"), keep.source = FALSE)
  def <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
                  identical(x[[2]], as.name(name)), exprs)
  stopifnot(length(def) == 1L)
  body(eval(def[[1]][[3]]))
}

## every call to `fn` anywhere in `expr` (a `|>` pipe is already a nested call once parsed)
callsTo <- function(expr, fn) {
  out <- list()
  if (!is.call(expr)) return(out)
  if (identical(expr[[1]], as.name(fn))) out <- list(expr)
  for (i in seq_along(expr)) {
    if (is.call(expr[[i]])) out <- c(out, callsTo(expr[[i]], fn))
  }
  out
}

## a cohort's B from its pixel's t/ha, as createBiomass_coreInputs() makes it: LandR::makePixelTable()
## sets totalBiomass, LandR::partitionBiomass() gives the cohort its `share`, then the class rounding
cohortB <- function(tPerHa, pixelGroupBiomassClass, share = 1) {
  totalBiomass <- LandR::asInteger(tPerHa * 100)
  LandR::asInteger(totalBiomass * share / pixelGroupBiomassClass) * pixelGroupBiomassClass
}

test_that("NorthernAlbertaUplands' largest cohort is within the limit; round(maxRawB, -2) was not", {
  maxRawB <- 329.4970703 * 100
  expect_identical(cohortB(329.4970703, 100), 33000)
  expect_identical(cohortBiomassLimit(maxRawB, 100), 33000)
  expect_lte(cohortB(329.4970703, 100), cohortBiomassLimit(maxRawB, 100))
  expect_identical(round(maxRawB, -2), 32900) ## the old limit
})

test_that("WesternAlbertaUpland's maximum rounds down, so the limit is the maximum itself", {
  maxRawB <- 304.4286499 * 100
  expect_identical(cohortB(304.4286499, 100), 30400)
  expect_identical(round(maxRawB, -2), 30400) ## the old limit, which passed here
  expect_identical(cohortBiomassLimit(maxRawB, 100), maxRawB)
  ## a young cohort re-estimated inside a fire perimeter: capped at maxRawB by
  ## youngBiomassLimit(), and not rounded to the class
  expect_lte(30442L, cohortBiomassLimit(maxRawB, 100))
})

test_that("the limit follows pixelGroupBiomassClass", {
  ## 32949.70703 -> 32950 -> 33000 in 1000 g/m2 classes, 32950 in 10 g/m2 classes
  expect_identical(cohortBiomassLimit(32949.70703, 1000), 33000)
  expect_identical(cohortB(329.4970703, 10), 32950)
  expect_identical(cohortBiomassLimit(32949.70703, 10), 32950)
  ## 32600 -> 33000 in 1000 g/m2 classes, 400 above round(maxRawB, -2)
  expect_identical(cohortB(326, 1000), 33000)
  expect_identical(cohortBiomassLimit(32600, 1000), 33000)
  ## 32400 -> 32000, below the young-cohort cap, maxRawB
  expect_identical(cohortBiomassLimit(32400, 1000), 32400)
})

test_that("no cohort made from the map is above the limit of the map's maximum", {
  withr::local_seed(20261026)
  ## the vector maxima step by 0.1 g/m2 over 10 t/ha, so they cross every rounding tie
  maxima <- seq(300, 310, by = 0.001)
  ## the last is the old LandWeb rule at 120 m pixels, 1000 / (250 / 120)^2 = 230.4
  for (cl in c(10, 100, 250, 1000, 1000 / (250 / 120)^2)) {
    excess <- vapply(maxima, function(m) {
      v <- c(m, stats::runif(49, 0, m))
      share <- c(1, stats::runif(49)) ## the largest pixel has a single species
      max(cohortB(v, cl, share)) - cohortBiomassLimit(max(v) * 100, cl)
    }, numeric(1))
    expect_lte(max(excess), 0, label = paste("largest excess with class", cl))
  }
})

test_that("createBiomass_coreInputs() checks cohort B against cohortBiomassLimit()", {
  bod <- moduleBody("createBiomass_coreInputs")
  checks <- Filter(function(x) identical(x[[2]], quote(na.omit(pixelCohortData$B))),
                   callsTo(bod, "inRange"))
  expect_length(checks, 1L)
  expect_identical(checks[[1]][[4]],
                   quote(cohortBiomassLimit(maxRawB, P(sim)$pixelGroupBiomassClass)))
  expect_length(Filter(function(x) identical(x, quote(round(maxRawB, -2))), callsTo(bod, "round")), 0L)
})
