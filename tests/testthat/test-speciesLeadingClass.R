## Composition of a class-240 pixel comes from its species cover, not from another year's land
## cover or a neighbour's class:
##   conifer share of tree cover >= threshold  -> conifer leading
##                              <= 1 - threshold -> deciduous leading
##                                in between   -> mixed
##
## The DEFAULT threshold is 0.75, which is NTEMS' own: EOSD (Wulder & Nelson 2003) defines
## coniferous/broadleaf as 75% or more of total basal area and mixed wood as neither reaching
## 75%. The rule has the same shape as LandR::vegTypeMapGenerator(mixedType = 2), and the
## cross-check below measures that agreement directly -- but at 0.8, which is what
## vegTypeMapGenerator means by leading vegetation. Tests that exercise the boundary therefore
## pass `vegLeadingProportion` explicitly; a test that omits it is testing the 0.75 default.

sppEq <- data.table::data.table(
  LandR = c("Pice_mar", "Pice_gla", "Pinu_ban", "Abie_bal", "Lari_lar", "Betu_pap", "Popu_tre"),
  Type = c(rep("Conifer", 5), rep("Deciduous", 2))
)

test_that("conifer, deciduous and mixed are separated at the 0.8 / 0.2 thresholds", {
  withr::local_package("data.table")
  ## conifer share: 1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.2, 0.1, 0.0
  cov <- data.table(
    Pice_mar = c(100, 90, 80, 70, 50, 30, 20, 10, 0),
    Popu_tre = c(0, 10, 20, 30, 50, 70, 80, 90, 100)
  )
  expect_identical(
    speciesLeadingClass(cov, sppEq, vegLeadingProportion = 0.8),
    c(210L, 210L, 210L, 230L, 230L, 230L, 220L, 220L, 220L)
  )
})

test_that("both thresholds are inclusive, and survive floating point", {
  withr::local_package("data.table")
  ## `1 - 0.8` is 0.19999999999999996 but `20/100` is 0.20000000000000001, so an exact `<=`
  ## sends a true 20/80 split to mixed. vegTypeMapGenerator() calls it deciduous-leading.
  cov <- data.table(
    Pice_mar = c(81, 80, 79, 21, 20, 19),
    Popu_tre = c(19, 20, 21, 79, 80, 81)
  )
  expect_identical(
    speciesLeadingClass(cov, sppEq, vegLeadingProportion = 0.8),
    c(210L, 210L, 230L, 230L, 220L, 220L)
  )
})

test_that("a pixel with no tree cover yields NA, for the caller's fallback", {
  withr::local_package("data.table")
  cov <- data.table(Pice_mar = c(0, 40), Popu_tre = c(0, 10))
  expect_identical(speciesLeadingClass(cov, sppEq), c(NA_integer_, 210L))
})

test_that("several conifers and deciduous species are pooled by Type", {
  withr::local_package("data.table")
  ## 30 + 30 conifer vs 30 + 10 deciduous -> 0.6 conifer -> mixed
  cov <- data.table(Pice_mar = 30, Pinu_ban = 30, Betu_pap = 30, Popu_tre = 10)
  expect_identical(speciesLeadingClass(cov, sppEq), 230L)
})

test_that("the 'cover.' prefix that makePixelTable leaves is handled", {
  withr::local_package("data.table")
  cov <- data.table(pixelIndex = 1:2, cover.Pice_mar = c(90, 5), cover.Popu_tre = c(10, 95))
  expect_identical(speciesLeadingClass(cov, sppEq), c(210L, 220L))
})

test_that("columns that are not species are ignored", {
  withr::local_package("data.table")
  cov <- data.table(pixelIndex = 1:2, lcc = c(240, 240),
                    Pice_mar = c(90, 5), Popu_tre = c(10, 95))
  expect_identical(speciesLeadingClass(cov, sppEq), c(210L, 220L))
})

test_that("no recognised species at all gives NA rather than a wrong class", {
  withr::local_package("data.table")
  cov <- data.table(Quer_rub = c(50, 50))
  expect_identical(speciesLeadingClass(cov, sppEq), c(NA_integer_, NA_integer_))
})

test_that("it agrees with vegTypeMapGenerator on the same splits", {
  skip_if_not_installed("LandR")
  withr::local_package("data.table")
  withr::local_package("terra")
  fracs <- c(1.0, 0.9, 0.81, 0.8, 0.79, 0.7, 0.5, 0.3, 0.21, 0.2, 0.19, 0.1, 0.0)
  ## 0.8 explicitly: this measures that the RULE matches vegTypeMapGenerator, which means
  ## leading vegetation at 0.8. The 0.75 default is a separate claim, tested on its own below.
  mine <- speciesLeadingClass(
    data.table(Pice_mar = fracs * 100, Popu_tre = (1 - fracs) * 100), sppEq,
    vegLeadingProportion = 0.8
  )
  pgm <- terra::rast(nrows = 1, ncols = 1, vals = 1L)
  theirs <- vapply(fracs, function(f) {
    x <- data.table(pixelGroup = 1L, speciesCode = c("Pice_mar", "Popu_tre"),
                    B = c(f * 100, (1 - f) * 100))
    o <- suppressMessages(LandR::vegTypeMapGenerator(
      x, pixelGroupMap = pgm, vegLeadingProportion = 0.8, mixedType = 2, doAssertion = FALSE
    ))
    lv <- terra::levels(o)[[1]]
    lab <- as.character(lv[[2]][match(terra::values(o, mat = FALSE)[1], lv[[1]])])
    if (lab == "Mixed") 230L else if (lab == "Pice_mar") 210L else 220L
  }, integer(1))
  expect_identical(mine, theirs)
})

## `deciduousCoverWeight` makes deciduous cover conifer-equivalent before the ratio is taken.
## partitionBiomass() applies it as `cover * c(1, x)[decid + 1]` -- deciduous is MULTIPLIED, i.e.
## shrunk -- which was checked directly: equal 50/50 cover gives conifer/deciduous biomass of
## 1.1878 = 1/0.8419. Dividing would inflate deciduous instead.
test_that("the deciduous discount shrinks deciduous cover, pushing the conifer share up", {
  withr::local_package("data.table")
  disc <- 0.8418911
  ## 45 conifer / 55 deciduous: raw share 0.45 (mixed); discounted 45/(45+46.3) = 0.493, still mixed
  cov <- data.table(Pice_mar = 45, Popu_tre = 55)
  expect_identical(speciesLeadingClass(cov, sppEq), 230L)
  expect_identical(speciesLeadingClass(cov, sppEq, deciduousCoverWeight = disc), 230L)

  ## A pixel just below the conifer threshold on raw cover crosses it once deciduous is shrunk.
  ## Pinned at 0.8 so the demonstration is about the discount, not about the default threshold:
  ## raw 0.780 is already conifer-leading under the 0.75 default.
  cov2 <- data.table(Pice_mar = 78, Popu_tre = 22)          ## raw 0.780 -> mixed at 0.8
  expect_identical(speciesLeadingClass(cov2, sppEq, vegLeadingProportion = 0.8), 230L)
  expect_identical(
    speciesLeadingClass(cov2, sppEq, vegLeadingProportion = 0.8, deciduousCoverWeight = disc),
    210L
  )

  ## the default discount is a no-op
  expect_identical(
    speciesLeadingClass(cov2, sppEq, vegLeadingProportion = 0.8, deciduousCoverWeight = 1),
    speciesLeadingClass(cov2, sppEq, vegLeadingProportion = 0.8)
  )
})

## The default is the product's own threshold, not LandR's. Guards against it drifting back to
## 0.8 silently: 0.78 conifer is mixed wood to vegTypeMapGenerator but coniferous to NTEMS.
test_that("the default threshold is NTEMS' 0.75, not vegTypeMapGenerator's 0.8", {
  withr::local_package("data.table")
  ## Pinned, so this tests the fallthrough rather than whatever the session happens to have
  ## set. The default is a call now -- getOption nested over getOption -- so it needs eval().
  withr::local_options(list(NTEMS.mixedwoodProp = NULL, LandR.lccLeadingProportion = NULL))
  expect_identical(eval(formals(speciesLeadingClass)$vegLeadingProportion), 0.75)

  ## conifer share 0.78 and 0.74, and their deciduous mirrors
  cov <- data.table(Pice_mar = c(78, 74, 26, 22), Popu_tre = c(22, 26, 74, 78))
  expect_identical(speciesLeadingClass(cov, sppEq), c(210L, 230L, 230L, 220L))
  expect_identical(speciesLeadingClass(cov, sppEq, vegLeadingProportion = 0.8),
                   c(230L, 230L, 230L, 230L))
})

## The module helper honours the same one-knob option as LandR, so a user who re-points the
## ecosystem does not have to remember this function separately.
test_that("NTEMS.mixedwoodProp re-points the fill without touching the call site", {
  withr::local_package("data.table")
  cov <- data.table(Pice_mar = 74, Popu_tre = 26)          ## 0.74: mixed under the 0.75 default
  withr::with_options(list(NTEMS.mixedwoodProp = NULL, LandR.lccLeadingProportion = NULL),
                      expect_identical(speciesLeadingClass(cov, sppEq), 230L))
  withr::with_options(list(NTEMS.mixedwoodProp = 0.7),
                      expect_identical(speciesLeadingClass(cov, sppEq), 210L))
})
