## A maxB the fit clamped to 0 is replaced by the 95th percentile of the biomass observed for the
## species in its group; genuine fits are never touched (an upper quantile of many cohorts is
## still an extreme statistic). A percentile, not the maximum, so one freak cohort cannot set the
## ceiling -- the same shape as the module's `quantile(age, 0.99)` longevity rule.
library(data.table)

se <- data.table(ecoregionGroup = factor(c("1_081", "1_081", "1_210")),
                 speciesCode = factor(c("Betu_pap", "Pice_mar", "Pice_mar")),
                 maxB = c(0L, 3200L, 5400L), maxANPP = c(0L, 107L, 180L), establishprob = 0.5)
cd <- data.table(ecoregionGroup = c("1_081", "1_081", "1_081", "1_210"),
                 speciesCode = c("Betu_pap", "Betu_pap", "Pice_mar", "Pice_mar"),
                 B = c(300L, 700L, 2900L, 6000L))

test_that("only a clamped (0) maxB is replaced, with maxANPP to match", {
  out <- floorMaxBAtObserved(se, cd)
  q <- LandR::asInteger(quantile(c(300, 700), 0.95, names = FALSE))   # 680
  expect_identical(out$maxB, c(q, 3200L, 5400L))       # 5400 < observed 6000 but is a real fit
  expect_identical(out$maxANPP, c(LandR::asInteger(q / 30), 107L, 180L))
  expect_identical(attr(out, "nRaised"), 1L)
})

test_that("a single wild outlier does not set the ceiling", {
  ## A realistic group: many cohorts around 300-800, plus one freak 90,000.
  many <- data.table(ecoregionGroup = "1_081", speciesCode = "Betu_pap",
                     B = c(seq(300L, 800L, by = 10L), 90000L))
  out <- floorMaxBAtObserved(se, rbind(cd[speciesCode != "Betu_pap"], many))
  expect_lt(out[speciesCode == "Betu_pap", maxB], 1000L)      # the outlier is ignored
  expect_gt(out[speciesCode == "Betu_pap", maxB], 700L)       # but the upper tail is respected

  ## The maximum, which this replaced, would have taken the outlier whole.
  outMax <- floorMaxBAtObserved(se, rbind(cd[speciesCode != "Betu_pap"], many), probs = 1)
  expect_identical(outMax[speciesCode == "Betu_pap", maxB], 90000L)
})

test_that("with few observations the percentile necessarily sits near the maximum", {
  ## Not a defect, a property of sample quantiles: with 3 values the 95th percentile is close to
  ## the largest of them. The protection this gives is real only for well-populated groups.
  cdFew <- data.table(ecoregionGroup = "1_081", speciesCode = "Betu_pap", B = c(300L, 700L, 9000L))
  out <- floorMaxBAtObserved(se, rbind(cd[speciesCode != "Betu_pap"], cdFew))
  expect_gt(out[speciesCode == "Betu_pap", maxB], 8000L)
})

test_that("probs is configurable, and probs = 1 is the maximum", {
  out <- floorMaxBAtObserved(se, cd, probs = 1)
  expect_identical(out[speciesCode == "Betu_pap", maxB], 700L)
  expect_error(floorMaxBAtObserved(se, cd, probs = c(0.5, 0.95)))
})

test_that("a genuine fit below the observed maximum is never raised", {
  out <- floorMaxBAtObserved(se, cd)
  expect_identical(out[ecoregionGroup == "1_210", maxB], 5400L)
})

test_that("rows already above the observations, and the other columns, are untouched", {
  out <- floorMaxBAtObserved(se, cd)
  expect_identical(out[speciesCode == "Pice_mar" & ecoregionGroup == "1_081", maxB], 3200L)
  expect_identical(out$establishprob, se$establishprob)
  expect_identical(names(out), names(se))
  expect_identical(levels(out$ecoregionGroup), levels(se$ecoregionGroup))
})

test_that("a group-species with one observation floors at that observation", {
  cdOne <- cd[!(speciesCode == "Betu_pap" & B == 300L)]
  out <- floorMaxBAtObserved(se, cdOne)
  expect_identical(out[speciesCode == "Betu_pap", maxB], 700L)
})

test_that("a group-species with no observations keeps its fitted value", {
  out <- floorMaxBAtObserved(se, cd[speciesCode != "Betu_pap"])
  expect_identical(out[speciesCode == "Betu_pap", maxB], 0L)
})

test_that("the input table is not modified", {
  before <- copy(se)
  invisible(floorMaxBAtObserved(se, cd))
  expect_identical(se, before)
})
