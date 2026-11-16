## `Cache()` does not key on RNG state. A call that draws random numbers and is cached therefore
## replays its first draw under another `.seed`. The two cached calls of `init` that draw
## (makeAndCleanInitialCohortData: LandR::subsetDT for the age model; convertUnwantedLCC with
## "nearestRandom") must have the seed of `init` in their cache key.

test_that("the cached random draws of init have the .seed of init in their key", {
  src <- parse(file.path(moduleRoot, paste0(moduleName, ".R")), keep.source = FALSE)
  cacheCalls <- list()
  walk <- function(e) {
    if (is.call(e)) {
      if (identical(e[[1]], as.name("Cache")) && length(e) > 1) cacheCalls[[length(cacheCalls) + 1]] <<- e
      for (a in as.list(e)) if (!missing(a)) walk(a)
    }
  }
  for (e in as.list(src)) walk(e)
  draws <- c("makeAndCleanInitialCohortData", "convertUnwantedLCC")
  for (fn in draws) {
    hit <- Filter(function(cl) is.call(cl[[2]]) && identical(cl[[2]][[1]], as.name(fn)), cacheCalls)
    expect_length(hit, 1L)
    extra <- hit[[1]][[".cacheExtra"]]
    expect_identical(extra, as.name("initSeed"), info = fn)
  }
  ## the biomass model is keyed on the total biomass of the whole data, not on its random subset
  stats <- Filter(function(cl) identical(cl[[2]], as.name("statsModel")) && !is.null(cl[[".cacheExtra"]]) &&
                    identical(cl[["modelFn"]], as.name("modelFn")), cacheCalls)
  expect_length(stats, 1L)
  expect_true("initSeed" %in% all.names(stats[[1]][[".cacheExtra"]]))
  ## and initSeed is the .seed of init, NULL (no key change) when none is set
  txt <- paste(deparse(src), collapse = "\n")
  expect_match(txt, "initSeed <- P\\(sim\\)\\$\\.seed\\[\\[\"init\"\\]\\]", fixed = FALSE)
})
