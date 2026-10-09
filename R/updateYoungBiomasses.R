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
    if (!is(modelBiomass$scaledVarsModelB, "list")) {
      stop("modelBiomass$scaledVarsModelB must be a list")
    }

    if (!all(names(modelBiomass$scaledVarsModelB) %in% c("cover", "logAge"))) {
      stop("modelBiomass$scaledVarsModelB must be a list with 'cover' and 'logAge' entries")
    }

    setnames(young, c("logAge", "cover"), c("logAge_orig", "cover_orig")) ## original, unscaled vars
    young[, `:=`(logAge = scale(logAge_orig,                            # nolint: conflicting_fn_unqualified
                                center = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:center"),
                                scale = attr(modelBiomass$scaledVarsModelB$logAge, "scaled:scale")),
                 cover = scale(cover_orig,                              # nolint: conflicting_fn_unqualified
                               center = attr(modelBiomass$scaledVarsModelB$cover, "scaled:center"),
                               scale = attr(modelBiomass$scaledVarsModelB$cover, "scaled:scale")))]
  }

  if (is(modelBiomass$mod, "merMod")) {
    columns <- c("ecoregionGroup", "logAge", "speciesCode", "cover")
    young2 <- unique(young, by = columns)
    message(green("  -- Calculating bootstrap estimates around B; will replace B in young data if it is beyond 95% CI"))
    message(green("     This will take some time."))
    ## Not wrapped in system.time(): it defaults to gcFirst = TRUE, so this forced a
    ## full garbage collection every call, and `PI.time` was never read.
    PI <- predictInterval(merMod = modelBiomass$mod, newdata = young2,
                          level = 0.95, n.sims = 15,
                          stat = "median", type = "linear.prediction",
                          include.resid.var = TRUE)
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


#' Update Biomass using a spinup type approach
#'
#' This function takes a pixelCohortData, and updates all the B values according to
#' the succession and growth dynamics provided by \code{species} and \code{speciesEcoregion}
#' objects. It does this by setting B to zero in each pixelIndex, growing each pixelIndex
#' until it hits the age given by age, for each pixelIndex independently. This uses
#' Biomass_core module to do this. If this module does not exist in the modules(sim),
#' then it will download it and use the latest development branch version of Biomass_core.
#'
#' @param pixelCohortData A pixelCohortData object (must have columns pixelIndex,
#'                         ecoregionGroup, age, and B). Should only have rows where the
#'                         age <= maxAge.
#'                         It should have none where age = 0
#' @param speciesEcoregion A speciesEcoregion object (must have ecoregionGroup, maxB, maxANPP)
#' @param maxAge A numeric scalar, indicating the maximum age to simulation biomass.
#'                              This should correspond to the maximum age in pixelCohortData
#' @param minRelativeB A minRelativeB object
#' @param species A species table with species-level parameters
#' @param sppColorVect A sppColorVect object (required, but not used)
#' @param paths A list of spades paths e.g., from paths(sim)
#' @param currentModule A character string of the current module e.g., from currentModule(sim)
#' @param modules A list of character strings of the modules in the sim, e.g., from modules(sim)
#' @param initialB Biomass (g/m2) at which every cohort starts the spin-up. Also passed to
#'   Biomass_core as its `initialB`, i.e. the biomass a newly established cohort gets there
#'   (`LandR:::.initiateNewCohorts`). The default matches Biomass_core's default. If `NA`,
#'   cohorts start at `maxANPP` from `speciesEcoregion`, as Biomass_core does.
#' @export
spinUpPartial <- function(
  pixelCohortData,
  speciesEcoregion,
  maxAge,
  # rasterToMatch, speciesLayers,
  minRelativeB,
  species,
  sppEquiv,
  sppEquivCol,
  sppColorVect,
  paths,
  currentModule,
  modules,
  initialB = 10
) {
  rng <- range(pixelCohortData$age)
  if (rng[1] <= 0) {
    stop("This spinup is only tested with age > 0")
  }
  if (rng[2] > maxAge) {
    stop("This spinup is only tested with age <= maxAge")
  }
  cd <- copy(pixelCohortData)
  cd[, `:=`(pixelGroup = as.integer(factor(pixelIndex)))]
  pixelGroupMap <- pixelGroupMapGenerate(cd)
  # Maps
  studyArea <- vect(ext(pixelGroupMap), crs(pixelGroupMap))
  rasterToMatch <- pixelGroupMap
  ecoregionMap <- pixelGroupMap
  levels(ecoregionMap) <- data.frame( # nolint: conflicting_fn_unqualified
    ID = 1:max(cd$pixelGroup, na.rm = TRUE),
    ecoregion = 1,
    ecoregionGroup = 1,
    stringsAsFactors = TRUE
  )
  # minRelativeB <- sim$minRelativeB
  ecoregion <- makeEcoregionDT(cd, speciesEcoregion)
  parameters <- list(
    Biomass_core = list(
      .saveInitialTime = NA,
      .saveInterval = NA,
      .useParallel = 1,
      seedingAlgorithm = "noSeeding",
      calcSummaryBGM = NULL,
      .plots = NULL,
      .maxMemory = 1e9,
      sppEquivCol = sppEquivCol,
      .useCache = NULL,
      successionTimestep = 10,
      minCohortBiomass = 0,
      initialB = initialB,
      initialBiomassSource = "cohortData",
      vegLeadingProportion = 0
    )
  )
  #sppEquiv needed or module stops, but object unused, likewise with speciesLayers
  speciesLayers <- "species"

  times <- list(start = 0, end = as.numeric(maxAge))

  speciesEcoregion2 <- copy(speciesEcoregion)
  speciesEcoregion2[, year := times$start]
  cdZeroed <- spinUpStartCohorts(cd, initialB, speciesEcoregion)
  objectsForYoungSim <- list(
    studyArea = studyArea,
    rasterToMatch = rasterToMatch,
    cohortData = cdZeroed,
    species = species,
    speciesEcoregion = speciesEcoregion2,
    pixelGroupMap = pixelGroupMap,
    speciesLayers = speciesLayers,
    minRelativeB = minRelativeB,
    ecoregion = ecoregion,
    ecoregionMap = ecoregionMap,
    sppEquiv = sppEquiv,
    sppColorVect = sppColorVect
  )
  opts <- options(
    "LandR.assertions" = FALSE,
    "LandR.verbose" = 0,
    "spades.recoveryMode" = FALSE,
    "spades.moduleCodeChecks" = FALSE # Turn off all module's code checking
  )

  ## TODO: make the following into a function to use across modules (e.g. B_sppFactorial)
  curModPath <- file.path(paths$modulePath, currentModule)
  curModPath <- curModPath[dir.exists(curModPath)] # for multiple module paths
  submodulePath <- file.path(curModPath, "submodules") |> checkPath(create = TRUE)
  paths$outputPath <- file.path(submodulePath, "outputs", rndstr()) ## avoid race conditions
  on.exit(unlink(paths$outputPath, recursive = TRUE), add = TRUE)

  bcVersion <- "1.3.10"
  ## if Biomass_core doesn't exist in modulePath or is too old, then download it
  modulesInProject <- list.dirs(paths$modulePath, full.names = TRUE, recursive = FALSE) |> as.list()
  names(modulesInProject) <- modulesInProject
  modulesInProject <- lapply(modulesInProject, basename)
  Biomass_core_path <- paths$modulePath[dir.exists(file.path(paths$modulePath, "Biomass_core"))]
  BCore_missingOrOld <- TRUE
  if (length(Biomass_core_path) > 0) {
    if (moduleVersion("Biomass_core", Biomass_core_path) >= bcVersion) {
      ## trim unnecessary modules:
      BCore_missingOrOld <- FALSE
      # modules <- modules[modules == "Biomass_core"]
    }
  }
  if (BCore_missingOrOld) {
    ## NOTE: don't install pkgs mid-stream; use module metadata to declare pkgs for installation
    moduleNameAndBranch <- paste0("PredictiveEcology/Biomass_core@development (>= ", bcVersion, ")")
    modules <- Require::extractPkgName(moduleNameAndBranch)
    paths$modulePath <- file.path(submodulePath, "Biomass_core")
    getModule(moduleNameAndBranch, modulePath = paths$modulePath, overwrite = TRUE) ## will only overwrite if wrong version
  } else {
    modules <- "Biomass_core"
  }

  outputs <- data.frame(expand.grid(
    objectName = "cohortData",
    saveTime = unique(seq(times$start, times$end, by = 1)),
    eventPriority = 1,
    fun = "qs2::qs_save",
    stringsAsFactors = FALSE
  ))
  suppressMessages({
    ss <- simInit(paths = paths, outputs = outputs, times = times)
  })
  outputs <- outputs(ss)
  mySimOut <- simInitAndSpades(
    # .cacheExtra = list(knownDigest, paths$outputPath),
    # omitArgs = c("objects", "params", "debug", "paths"),
    times = times,
    params = parameters,
    modules = modules, # quick = "paths",
    paths = paths,
    objects = objectsForYoungSim,
    outputs = outputs,
    # outputObjects = "pixelGroupMap",
    debug = 1
  )
  cds <- ReadExperimentFiles(outputs)
  cd1 <- cds[cd[, -"B"], on = c("pixelGroup", "speciesCode", "age"), nomatch = NA]
  set(cd1, NULL, c("pixelGroup"), NULL)
  setcolorder(cd1, neworder = colnames(pixelCohortData))

  return(cd1[])
}

## Cohorts at the start of the spin-up: age 1 at the biomass a new cohort gets in Biomass_core.
## Starting at B = 0 left growthcurve-1 species (e.g. spruce) with only the 1 g/m2/yr growth
## floor, so they ended at B = age - 1 while pioneers took the pixel.
spinUpStartCohorts <- function(cd, initialB, speciesEcoregion) {
  cd <- copy(cd)
  if (is.null(initialB) || is.na(initialB)) {
    ## As LandR::.initiateNewCohorts (cohorts.R, `isTRUE(is.na(initialB))` branch) with sumB = 0:
    ## asInteger(pmin(maxANPP, asInteger(pmax(1, maxANPP * exp(-1.6 * 0 / maxB_eco)))))
    se <- unique(speciesEcoregion[, c("speciesCode", "ecoregionGroup", "maxANPP")],
                 by = c("speciesCode", "ecoregionGroup"))
    se[, startB := asInteger(pmin(maxANPP, asInteger(pmax(1, maxANPP))))]
    cd[, `:=`(age = 1L, B = se$startB[match(paste(speciesCode, ecoregionGroup),
                                           paste(se$speciesCode, se$ecoregionGroup))])]
  } else {
    cd[, `:=`(age = 1L, B = as.integer(initialB))]
  }
  cd[]
}

pixelGroupMapGenerate <- function(cohortData) {
  pixelGroupMap <- rast(res = c(1, 1))
  nrow(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0)
  ncol(pixelGroupMap) <- round(sqrt(max(cohortData$pixelGroup)), 0) + 1
  vals <- c(1:max(cohortData$pixelGroup), rep(NA, times = ncell(pixelGroupMap) - max(cohortData$pixelGroup)))
  pixelGroupMap[] <- vals
  pixelGroupMap
}

ReadExperimentFiles <- function(outputs) {
  outputs <- as.data.table(outputs)[objectName == "cohortData"]
  fEs <- .fileExtensions()
  cdsList <- by(outputs, outputs[, "saveTime"], function(x) {
    fE <- reproducible:::fileExt(x$file)
    wh <- fEs[fEs$exts %in% fE,]
    message(crayon::green("reading: "))
    cat(crayon::green(x$file, "..."))
    cd <- getFromNamespace(wh$fun, ns = asNamespace(wh$package))(x$file)[, .(speciesCode, age, B, pixelGroup)]
    cat(crayon::green(" Done!\n"))
    return(cd)
  })
  message("rbindlisting the cohortData objects")
  cds <- rbindlist(cdsList, use.names = TRUE, fill = TRUE)

  return(invisible(cds))
}

#' Upper limit for the biomass of young cohorts re-estimated inside fire perimeters
#'
#' The age-based limit (2.8 x the raw map maximum, scaled by how much of a species' longevity
#' the young cohort has lived) is widened by 1.2 and capped at the raw map maximum.
#' The 1.2 is there because fast-growing pioneers legitimately exceed the age-based limit while
#' still below their own `maxB`: a 59-year Lari_lar cohort was 17% over (3,290 vs 2,817 g/m2)
#' and Pinu_con cohorts up to 2.9% over (ELF 3.2.3 and 10.3.1, 2026-10-09).
#' A single scalar limit is kept on purpose, not a per-species x ecoregion `maxB`.
#'
#' @param maxRawB Maximum of the raw biomass map, in g/m2.
#' @param longevity Species longevities (`sim$species$longevity`).
#' @param maxAgeHighQualityData Age up to which ages (and so biomass) are considered high quality.
#' @return A numeric scalar.
#' @noRd
youngBiomassLimit <- function(maxRawB, longevity, maxAgeHighQualityData) {
  min(maxRawB, 1.2 * 2.8 * maxRawB / min(longevity / maxAgeHighQualityData))
}
