partitionBiomass <- function(x, pixelCohortData) {
  if (!"decid" %in% colnames(pixelCohortData)) {
    pixelCohortData[, decid := speciesCode %in% c("Popu_Tre", "Betu_Pap")]
  }

  pixelCohortData[, cover2 := cover * c(1,x)[decid + 1]]
  pixelCohortData[, cover2 := cover2/sum(cover2), by = "pixelIndex"]
  pixelCohortData[, B := totalBiomass*cover2]
  pixelCohortData

}
