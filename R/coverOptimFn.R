#' @importFrom crayon cyan
#' @importFrom LandR subsetDT asInteger
#' @importFrom data.table setDT

coverOptimFn <- function(x, pixelCohortData, subset, bm, returnRsq = TRUE) {

  pixelCohortData <- partitionBiomass(x, pixelCohortData)
  if (length(subset) > 1) {
    pixelCohortData2 <- pixelCohortData[subset]
  } else {
    pixelCohortData2 <- subsetDT(pixelCohortData, c("initialEcoregionCode", "speciesCode"),
                                 subset)
  }
  pixelCohortData2 <- pixelCohortData2[!is.infinite(pixelCohortData2$logAge)]
  pixelCohortData2 <- pixelCohortData2[pixelCohortData2$B > 0]

  modelBiomass1 <-
    statsModel(
      modelFn = bm,
      uniqueEcoregionGroup = .sortDotsUnderscoreFirst(as.character(unique(pixelCohortData2$initialEcoregionGroup))),
      .specialData = pixelCohortData2#,
    )
  theAIC <- AIC(modelBiomass1$mod)
  message(cyan("#########################"))
  message(cyan(" -- deciduousDiscount:", round(x, 3), "; AIC=", round(theAIC, 3)))
  messageDF(modelBiomass1$rsq, round = 4, colour = "cyan")
  if (returnRsq)
    theAIC#unname(modelBiomass1$rsq[,2])
  else
    list(modelBiomass1 = modelBiomass1, pixelCohortData = pixelCohortData)
}

partitionBiomass <- function(x, pixelCohortData) {
  if (!"decid" %in% colnames(pixelCohortData)) {
    pixelCohortData[, decid := speciesCode %in% c("Popu_Tre", "Betu_Pap")]
  }

  pixelCohortData[, cover2 := cover * c(1,x)[decid + 1]]
  pixelCohortData[, cover2 := cover2/sum(cover2), by = "pixelIndex"]
  pixelCohortData[, B := totalBiomass*cover2]
  pixelCohortData

}
