#' @param young data.table not unlinke `cohortData`
#' @param modelBiomass named list with items "mod", "pred", "rsq", "scaledVarsModelB".
#' @param ... For anything, used for Cache. Not used internally here.
#'
#' @importFrom crayon green
#' @importFrom data.table setDT set setnames
#' @importFrom merTools predictInterval
#' @importFrom LandR asInteger
updateYoungBiomasses <- function(young, modelBiomass, ...) {
  young <- copy(young)

  useRescaled <- !is.null(modelBiomass$scaledVarsModelB)

  if (useRescaled) {
    if (!is(modelBiomass$scaledVarsModelB, "list"))
      stop("modelBiomass$scaledVarsModelB must be a list")

    if (!all(names(modelBiomass$scaledVarsModelB) %in% c("cover", "logAge")))
      stop("modelBiomass$scaledVarsModelB must be a list with 'cover' and 'logAge' entries")

    setnames(young, c("logAge", "cover"), c("logAge_orig", "cover_orig")) ## original, unscaled vars
    young[, `:=`(logAge = scale(logAge_orig,
                                center = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:center"),
                                scale = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:scale")),
                 cover = scale(cover_orig,
                               center = attr(modelBiomass$scaledVarsModelB$cover, "scaled:center"),
                               scale = attr(modelBiomass$scaledVarsModelB$cover, "scaled:scale")))]
  }

  if (is(modelBiomass$mod, "merMod")) {
    columns <- c("ecoregionGroup", "logAge", "speciesCode", "cover")
    young2 <- unique(young, by = columns)
    message(green("  -- Calculating bootstrap estimates around B; will replace B in young data if it is beyond 95% CI"))
    message(green("     This will take some time."))
    PI.time <- system.time({
      PI <- predictInterval(merMod = modelBiomass$mod, newdata = young2,
                            level = 0.95, n.sims = 15,
                            stat = "median", type = "linear.prediction",
                            include.resid.var = TRUE)
    })
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
    pres <- predict(modelBiomass$mod, newdata = young, se = TRUE, type = "response")
    set(young, NULL, "pred", pres$fit)
    set(young, NULL, "se", pres$se.fit)
    young[, resid := B - pred]
    young[, beyond := abs(resid) > 2*se]
    young[, tooLarge := resid > 2*se & beyond]
    young[, tooSmall := resid < 2*se & beyond]
    young[tooLarge == TRUE, newB := pred + 2*se]
    young[tooSmall == TRUE, newB := pred - 2*se]
  }
  if (sum(young$beyond) > 0)
  message("Within the cohorts aged ", max(young$age)," and younger, ",
          "there were ", sum(young$beyond), " cohorts whose biomass was way out of line for their ages. ",
          "Their biomasses have been adjusted down if too high (or up if too low) ",
          "to their predicted mean +1.96se (or - if too low) based on the fitted biomass model")
  young[beyond == FALSE, newB := B]
  young[, B := asInteger(pmax(0, newB))]
  if (useRescaled) {
    set(young, NULL, c("logAge", "cover"), NULL) ## remove scaled cols
    setnames(young, c("logAge_orig", "cover_orig"), c("logAge", "cover")) ## replace orig, unscaled
  }
  young[]
}
