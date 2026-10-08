## A landscape built with a known deciduous cover weight, so the estimator can be checked against
## the answer rather than against itself.
##
##   totalBiomass = k(height, closure, ecoregion) * (C + w * D)
##
## `k` is deliberately made to depend on structure AND to be correlated with deciduous cover --
## that correlation is the site-quality confound the structural controls exist to remove, and an
## estimator without them returns something well above `w` on exactly this data (verified below).
makeDecidLandscape <- function(w, n = 4000L, sd = 0.05, siteConfound = 1.0, seed = 42L,
                               dRange = c(0, 1)) {
  set.seed(seed)
  d <- stats::runif(n, dRange[1L], dRange[2L])  ## deciduous share of cover
  ## richer sites carry more deciduous: structure rises with d
  height  <- pmin(26, pmax(2.5, 6 + 14 * siteConfound * d + stats::rnorm(n, 0, 2)))
  closure <- pmin(90, pmax(6, 35 + 35 * siteConfound * d + stats::rnorm(n, 0, 6)))
  eco     <- sample(c("2", "3", "4"), n, replace = TRUE)
  age     <- as.integer(pmax(5, stats::rnorm(n, 90, 25)))
  k <- exp(2 + 1.1 * log(height) + 0.9 * closure / 100 + c("2" = 0, "3" = 0.2, "4" = -0.15)[eco])
  totalBiomass <- as.integer(round(k * (100 * (1 - d) + w * 100 * d) *
                                     exp(stats::rnorm(n, 0, sd))))
  pcd <- data.table::rbindlist(list(
    data.table::data.table(pixelIndex = seq_len(n), speciesCode = "Pice_mar",
                           initialEcoregionCode = paste0(eco, "_210"), age = age,
                           totalBiomass = totalBiomass, cover = 100 * (1 - d)),
    data.table::data.table(pixelIndex = seq_len(n), speciesCode = "Popu_tre",
                           initialEcoregionCode = paste0(eco, "_210"), age = age,
                           totalBiomass = totalBiomass, cover = 100 * d)
  ))[cover > 0]
  ras <- function(v) {
    r <- terra::rast(nrows = 1, ncols = n, xmin = 0, xmax = n, ymin = 0, ymax = 1)
    terra::values(r) <- v
    r
  }
  list(pcd = pcd[], height = ras(height), closure = ras(closure))
}

test_that("the fit recovers a known deciduous cover weight", {
  skip_if_not_installed("terra")
  for (w in c(0.6, 0.85, 1.3)) {
    L <- makeDecidLandscape(w)
    est <- suppressMessages(
      deciduousCoverWeightFn(L$pcd, L$height, L$closure, minDeciduousPixels = 100L)
    )
    expect_equal(est, w, tolerance = 0.05,
                 label = paste0("estimate for w = ", w))
  }
})

test_that("without the structural controls the same data returns a biased answer", {
  ## This is the point of the controls: `k` is correlated with deciduous cover by construction,
  ## and an ecoregion-only fit credits that site effect to the species.
  skip_if_not_installed("terra")
  L <- makeDecidLandscape(0.85)
  pix <- L$pcd[, list(totalBiomass = totalBiomass[1L],
                      eco = sub("_.*", "", initialEcoregionCode[1L]),
                      D = sum(cover[speciesCode == "Popu_tre"]),
                      C = sum(cover[speciesCode == "Pice_mar"])),
               by = "pixelIndex"]
  y <- log(pix$totalBiomass)
  noControl <- stats::optimize(function(w) {
    r <- y - log(pix$C + w * pix$D)
    sum((r - stats::ave(r, pix$eco))^2)
  }, interval = c(0.05, 5))$minimum
  expect_gt(noControl, 1.5)     ## far above the true 0.85
})

test_that("it declines to answer when the landscape cannot identify it", {
  skip_if_not_installed("terra")
  ## almost no deciduous
  L <- makeDecidLandscape(0.85)
  L$pcd <- L$pcd[speciesCode == "Pice_mar" | pixelIndex <= 20L]
  expect_message(
    est <- deciduousCoverWeightFn(L$pcd, L$height, L$closure, minDeciduousPixels = 100L),
    "not enough deciduous cover"
  )
  expect_true(is.na(est))
})

test_that("it declines when every pixel has nearly the same deciduous share", {
  ## plenty of deciduous -- every pixel is 95-100% -- but no contrast between pixels, so the
  ## weight is absorbed by the intercept. ELF 10.1 (Manitoba parkland) is this case.
  skip_if_not_installed("terra")
  L <- makeDecidLandscape(0.85, dRange = c(0.95, 1))
  expect_message(
    est <- deciduousCoverWeightFn(L$pcd, L$height, L$closure, minDeciduousPixels = 100L),
    "barely varies between pixels"
  )
  expect_true(is.na(est))
})

test_that("an estimate on the search bound is reported as unidentified", {
  skip_if_not_installed("terra")
  L <- makeDecidLandscape(0.85)
  expect_message(
    est <- deciduousCoverWeightFn(L$pcd, L$height, L$closure,
                                  interval = c(1.5, 5), minDeciduousPixels = 100L),
    "hit the search bound"
  )
  expect_true(is.na(est))
})
