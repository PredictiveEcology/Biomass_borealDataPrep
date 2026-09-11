## ELF 3.1.2 has one tree species with cover 100 in every pixel, so every ecoregionGroup row had
## coverPres == coverNum. Those rows are left out of the cover model, which left nothing to fit,
## and glm() on the empty table stopped with "contrasts can be applied only to factors with 2 or
## more levels".
test_that("no cover model is fitted when every row has 100% presence", {
  cohortDataShort <- data.table::data.table(
    ecoregionGroup = factor(c("04_NA", "06_NA", "14_NA")), speciesCode = factor("Pice_mar"),
    coverPres = c(80L, 73L, 64L), coverNum = c(80L, 73L, 64L))
  called <- FALSE
  out <- estimateCoverModel(cohortDataShort, fitCover = function(cds) {
    called <<- TRUE
    stop("there is nothing to fit")
  })

  expect_false(called)
  expect_null(out$model)
  expect_identical(out$modelCover, c(1, 1, 1))
})

test_that("rows with 100% presence get probability 1 and the other rows are fitted", {
  cohortDataShort <- data.table::CJ(ecoregionGroup = sprintf("%02d_NA", 1:6),
                                    speciesCode = c("Pice_gla", "Pice_mar"))
  cohortDataShort[, `:=`(ecoregionGroup = factor(ecoregionGroup), speciesCode = factor(speciesCode),
                         coverNum = 50L, coverPres = c(50L, seq(10L, 32L, by = 2L))[seq_len(.N)])]
  fittedRows <- NULL
  out <- estimateCoverModel(cohortDataShort, fitCover = function(cds) {
    fittedRows <<- nrow(cds)
    list(mod = glm(cbind(coverPres, coverNum - coverPres) ~ speciesCode + ecoregionGroup,
                   family = binomial, data = cds))
  })

  expect_identical(fittedRows, nrow(cohortDataShort) - 1L)
  expect_false(is.null(out$model))
  expect_length(out$modelCover, nrow(cohortDataShort))
  expect_identical(out$modelCover[1], 1)
  expect_true(all(out$modelCover[-1] > 0 & out$modelCover[-1] < 1))
})

test_that("with no 100% presence rows the fitted model itself is returned", {
  cohortDataShort <- data.table::CJ(ecoregionGroup = sprintf("%02d_NA", 1:4),
                                    speciesCode = c("Pice_gla", "Pice_mar"))
  cohortDataShort[, `:=`(ecoregionGroup = factor(ecoregionGroup), speciesCode = factor(speciesCode),
                         coverNum = 40L, coverPres = seq(8L, 22L, by = 2L))]
  fit <- list(mod = "a fitted model")
  out <- estimateCoverModel(cohortDataShort, fitCover = function(cds) fit)

  expect_identical(out$model, fit)
  expect_identical(out$modelCover, fit)
})
