#' @importFrom crayon green
#' @importFrom data.table setDT set setnames
#' @importFrom merTools predictInterval
#' @importFrom LandR asInteger

updateYoungBiomasses <- function(young, biomassModel) {
  if (is(biomassModel, "merMod")) {
    columns <- c("ecoregionGroup", "logAge", "speciesCode", "cover")
    young2 <- unique(young, by = columns)
    message(green("  -- Calculating bootstrap estimates around B; will replace B in young data if it is beyond 95% CI"))
    message(green("     This will take a bit"))
    PI.time <- system.time(
      PI <- predictInterval(merMod = biomassModel, newdata = young2,
                            level = 0.95, n.sims = 15,
                            stat = "median", type="linear.prediction",
                            include.resid.var = TRUE)
    )
    PI <- setDT(PI)
    young2 <- cbind(PI, young2)
    setnames(young2, old = "fit", new = "pred")
    set(young2, NULL, setdiff(colnames(young2), c("pred", "upr", "lwr", columns)), NULL)
    young <- young2[young, on = columns]
    young[, resid := B - pred]
    young[, beyond := resid > upr | resid < lwr]
    young[, tooLarge := resid > upr & beyond]
    young[, tooSmall := resid < lwr & beyond]
    young[tooLarge == TRUE, newB := upr]
    young[tooSmall == TRUE, newB := lwr]
  } else {
    pres <- predict(biomassModel, newdata = young, se = TRUE, type = "response")
    set(young, NULL, "pred", pres$fit)
    set(young, NULL, "se", pres$se.fit)
    young[, resid := B - pred]
    young[, beyond := abs(resid) > 2*se]
    young[, tooLarge := resid > 2*se & beyond]
    young[, tooSmall := resid < 2*se & beyond]
    young[tooLarge == TRUE, newB := pred + 2*se]
    young[tooSmall == TRUE, newB := pred - 2*se]
  }
  young[beyond == FALSE, newB := B]
  young[, B := asInteger(pmax(0, newB))]
  young[]
}
