#' @importFrom crayon magenta
deciduousCoverDiscountFun <- function(pixelCohortData,
                                      coverPctToBiomassPctModel = quote(glm(I(log(B/100)) ~ logAge * I(log(totalBiomass/100)) * speciesCode * lcc)),
                                      subsetDataAgeModel) {

  message(magenta(paste0(format(coverPctToBiomassPctModel, appendLF = FALSE))))

  pixelCohortData[, lcc := as.factor(lcc)]

  plot.it <- FALSE
  sam <- subsetDT(pixelCohortData, by = c("speciesCode", "lcc"),
                  doSubset = subsetDataAgeModel,
                  indices = TRUE)
  pi <- unique(pixelCohortData[sam]$pixelIndex)
  sam <- which(pixelCohortData$pixelIndex %in% pi)

  system.time({
    out <- optimize(interval = c(0.1, 1), f = coverOptimFn, bm = coverPctToBiomassPctModel,
                    pixelCohortData = pixelCohortData, subset = sam, maximum = FALSE)
  })

  if (plot.it) {
    cover2BiomassModel <- coverOptimFn(out$minimum, pixelCohortData, subsetDataAgeModel,
                                       coverPctToBiomassPctModel, returnAIC = FALSE)
    sam1 <- sample(NROW(pixelCohortData), 1e5)
    dev()
    par(mfrow = c(1,2))
    plot(predict(cover2BiomassModel$modelBiomass1$mod,
                 newdata = cover2BiomassModel$pixelCohortData[sam1]),
         log(cover2BiomassModel$pixelCohortData$B / 100)[sam1], pch = ".")
    abline(a = 0, b = 1)

    cover2BiomassModel1 <- coverOptimFn(1, pixelCohortData, subsetDataAgeModel,
                                        coverPctToBiomassPctModel,
                                        returnAIC = FALSE)
    dev()
    plot(predict(cover2BiomassModel1$modelBiomass1$mod,
                 newdata = cover2BiomassModel1$pixelCohortData[sam1]),
         log(cover2BiomassModel1$pixelCohortData$B / 100)[sam1], pch = ".")
    abline(a = 0, b = 1)

    pcd <- pixelCohortData
    bb <- pcd[sample(sam)]
    cc <- bb[, cover3 := cover * c(1, out$minimum)[decid + 1]][
      , actualX := cover3 / sum(cover3) / (cover / 100), by = "pixelIndex"]
    setkey(cc, pixelIndex)
    mean(cc[speciesCode == "Popu_Tre"]$actualX)
  }
  return(out$minimum)
}



#' @importFrom crayon cyan
#' @importFrom LandR subsetDT asInteger
#' @importFrom data.table setDT
coverOptimFn <- function(x, pixelCohortData, subset, bm, returnAIC = TRUE) {
  pixelCohortData <- partitionBiomass(x, pixelCohortData)
  if (length(subset) > 1) {
    pixelCohortData2 <- pixelCohortData[subset]
  } else {
    pixelCohortData2 <- subsetDT(pixelCohortData, c("initialEcoregionCode", "speciesCode"), subset)
  }
  pixelCohortData2 <- pixelCohortData2[!is.infinite(pixelCohortData2$logAge)]
  pixelCohortData2 <- pixelCohortData2[pixelCohortData2$B > 0]

  modelBiomass1 <- statsModel(
    modelFn = bm,
    uniqueEcoregionGroups = .sortDotsUnderscoreFirst(as.character(unique(pixelCohortData2$initialEcoregionCode))),
    .specialData = pixelCohortData2#,
  )
  theAIC <- AIC(modelBiomass1$mod)
  message(cyan("#########################"))
  message(cyan(" -- deciduousDiscount:", round(x, 3), "; AIC=", round(theAIC, 3)))
  messageDF(modelBiomass1$rsq, round = 4, colour = "cyan")
  if (returnAIC)
    theAIC#unname(modelBiomass1$rsq[,2])
  else
    list(modelBiomass1 = modelBiomass1, pixelCohortData = pixelCohortData)
}
