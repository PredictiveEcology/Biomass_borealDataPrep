defineModule(sim, list(
  name = "Boreal_LBMRDataPrep",
  description = "A data preparation module for parameterizing LBMR from open data sources, within the Boreal forest of Canada",
  keywords = c("LandWeb", "LBMR"),
  authors = c(
    person("Yong", "Luo", email = "yong.luo@canada.ca", role = c("aut")),
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person(c("Ceres"), "Barros", email = "cbarros@mail.ubc.ca", role = c("ctb")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@friresearch.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(SpaDES.core = "0.2.3.9009", Boreal_LBMRDataPrep = numeric_version("1.3.3")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Boreal_LBMRDataPrep.Rmd"),
  reqdPkgs = list("data.table", "dplyr", "fasterize", "gdalUtils", "raster", "rgeos", "sp",
                  "PredictiveEcology/LandR@development", "lme4",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    #defineParameter("speciesPresence", "numeric", 50, NA, NA,
    #                "minimum percent cover required to classify a species as present"),
    # defineParameter("minNumPixelsToEstMaxBiomass", "integer", 100, NA, NA,
    #                 "When estimating maximum biomass by species and ecoregion, this number indicates the minimum number of pixels with data required before a maximum is estimated."),
    # defineParameter("quantileForMaxBiomass", "numeric", 0.99, NA, NA,
    #                 "When estimating maximum biomass by species and ecoregion, rather than take the absolute max(biomass), the quantile is taken. This gives the capacity to remove outliers."),
    defineParameter("biomassQuotedFormula", "name",
                    quote(biomass ~ logAge * speciesCode + (speciesCode | ecoregionCode : lcc) + cover),
                    NA, NA,
                    "This formula is for estimating biomass from landcover, ecoregionCode, speciesCode, age, and cover"),
    defineParameter("coverQuotedFormula", "name",
                    quote(cbind(coverPres, coverNum) ~ speciesCode + (1 | ecoregionCode : lcc)),
                    NA, NA,
                    "This formula is for estimating biomass from landcover, ecoregion, and speciesCode"),
    defineParameter("pixelGroupAgeClass", "numeric", P(sim)$successionTimestep, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same"),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same"),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events")
  ),
  inputObjects = bind_rows(
    expectsInput("biomassMap", "RasterLayer",
                 desc = "total biomass raster layer in study area, default is Canada national biomass map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureBiomass.tar"),
    expectsInput("ecoDistrict", "SpatialPolygonsDataFrame",
                 desc = "ecodistricts in study area, default is Canada national ecodistricts",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("ecoRegion", "SpatialPolygonsDataFrame",
                 desc = "ecoregions in study area, default is Canada national ecoregions",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/region/ecoregion_shp.zip"),
    expectsInput("ecoZone", "SpatialPolygonsDataFrame",
                 desc = "ecozones in study area, default is Canada national ecozones",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/zone/ecozone_shp.zip"),
    expectsInput("LCC2005", "RasterLayer",
                 desc = "2005 land classification map in study area, default is Canada national land classification in 2005",
                 #sourceURL = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip"),
                 sourceURL = "https://drive.google.com/file/d/1g9jr0VrQxqxGjZ4ckF6ZkSMP-zuYzHQC/view?usp=sharing"),
    expectsInput("rasterToMatch", "RasterLayer",
                 #desc = "this raster contains two pieces of information: Full study area with fire return interval attribute",
                 desc = "DESCRIPTION NEEDED",
                 sourceURL = NA),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "biomass percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("speciesTable", "data.table",
                 desc = "species attributes table, default is from Dominic and Yan's project",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon to use as the study area,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (larger area than studyArea) to use for parameter estimation,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput("standAgeMap", "RasterLayer",
                 desc = "stand age map in study area, default is Canada national stand age map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-StructureStandVolume.tar"),
    expectsInput("sufficientLight", "data.frame",
                 desc = "define how the species with different shade tolerance respond to stand shadeness")
  ),
  outputObjects = bind_rows(
    createsOutput("ecoDistrict", "", desc = ""), ## TODO: description and type needed
    createsOutput("ecoRegion", "", desc = ""),   ## TODO:description and type needed
    createsOutput("ecoregion", "data.table",
                  desc = "ecoregion look up table"),
    createsOutput("ecoregionMap", "RasterLayer",
                  desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table"),
    createsOutput("ecoZone", "", desc = ""),
    createsOutput("initialCommunities", "data.table",
                  desc = "initial community table"),
    createsOutput("initialCommunitiesMap", "RasterLayer",
                  desc = "initial community map that has mapcodes match initial community table"),
    createsOutput("minRelativeB", "data.frame",
                  desc = "define the cut points to classify stand shadeness"),
    createsOutput("notEnoughDataMaxBiomass", "data.table",
                  desc = "The collection of ecoregion-species combinations that don't have values for maxBiomass or maxANPP"),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
    createsOutput("speciesEcoregion", "data.table",
                  desc = "define the maxANPP, maxB and SEP change with both ecoregion and simulation time"),
    createsOutput("studyArea", "", desc = ""),
    createsOutput("speciesEstablishmentProbMap", "RasterStack",
                  paste("Species establishment probability as a map, ",
                        "by species. This is written to disk to save RAM space")),
    createsOutput("useCache", "logic",
                  desc = "define which the caching for spinup simulation should be used, default is TRUE")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Boreal_LBMRDataPrep <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    # names(sim$speciesLayers) <- equivalentName(names(sim$speciesLayers), sim$sppEquiv, "Latin_full")
    sim <- estimateParameters(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Boreal_LBMRDataPrep", "save")
  } else if (eventType == "save") {
    sim <- Save(sim)
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

## event functions
#   - follow the naming convention `modulenameEventtype()`;
#   - `modulenameInit()` function is required for initiliazation;
#   - keep event functions short and clean, modularize by calling subroutines from section below.

estimateParameters <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  message("1: Starting to estimate parameters in Boreal_LBMRDataPrep: ", Sys.time())
  cPath <- cachePath(sim)
  sim$ecoDistrict <- spTransform(sim$ecoDistrict, crs(sim$speciesLayers))
  sim$ecoRegion <- spTransform(sim$ecoRegion, crs(sim$speciesLayers))
  sim$ecoZone <- spTransform(sim$ecoZone, crs(sim$speciesLayers))

  sim$standAgeMap <- round(sim$standAgeMap / 20, 0) * 20 # use 20-year bins (#103)
  sim$standAgeMap[] <- asInteger(sim$standAgeMap[])

  ################################################################
  ## species traits inputs
  ################################################################
  message("1: Prepare 'species' table, i.e., species level traits", Sys.time())
  sim$species <- prepSpeciesTable(speciesTable = sim$speciesTable,
                                  speciesLayers = sim$speciesLayers,
                                  sppEquiv = sim$sppEquiv,
                                  sppEquivCol = P(sim)$sppEquivCol)

  # rasterToMatchBinary <- raster(sim$rasterToMatch)
  # rasterToMatchBinary[] <- NA
  # rasterToMatchBinary[!is.na(sim$rasterToMatch[])] <- 1

  ecoDistrictMap <- Cache(postProcess, sim$ecoDistrict, studyArea = sim$studyArea, filename2 = NULL)
  ecoRegionMap <- Cache(postProcess, sim$ecoRegion, studyArea = sim$studyArea, filename2 = NULL)

  ecoregionMap <- Cache(postProcess, sim$ecoDistrict, studyArea = sim$studyArea, filename2 = NULL)
  ecoregionstatus <- data.table(active = "yes",
                                ecoregion = 1:1031)
  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(ecoregionMap, sim$LCC2005),
                          ecoregionName = "ECODISTRIC",
                          ecoregionActiveStatus = ecoregionstatus,
                          rasterToMatch = sim$rasterToMatch,
                          userTags = "stable")

  # put together cohortData object --
  #  Round age to successionTimestep
  message(crayon::blue("Step 1: Round biomass to nearest P(sim)$pixelGroupAgeClass, which is", P(sim)$pixelGroupAgeClass))
  dt <- data.table(age = asInteger(ceiling(asInteger(sim$standAgeMap[]) / P(sim)$pixelGroupAgeClass) *
                                     P(sim)$pixelGroupAgeClass),
                   #asInteger(sim$standAgeMap[]),
                   logAge = log(sim$standAgeMap[]),
                   ecoregionCode = factor(factorValues2(ecoregionFiles$ecoregionMap,
                                                        ecoregionFiles$ecoregionMap[],
                                                        att = 5)),
                   #ecoRegion = factor(fasterize::fasterize(sf::st_as_sf(ecoRegionMap),
                   #                                         raster = sim$rasterToMatch, "ECOREGION")[]),
                   totalBiomass = asInteger(sim$biomassMap[]) * 100, # change units
                   cover = sim$speciesLayers[],
                   pixelIndex = seq(ncell(sim$standAgeMap)),
                   lcc = factor(sim$LCC2005[])
  )

  message(crayon::blue("This is the summary of the input data for age, ecoregionCode, biomass, speciesLayers"))
  print(summary(dt))
  message(crayon::blue("Step 2: rm NAs, i.e., no input data in some pixels, including outside studyArea"))
  dt <- na.omit(dt)
  dt <- dt[age > 0]
  message(crayon::blue("Step 3: rm age == 0, leaving", NROW(dt), "pixels with data"))

  ### Create groupings
  coverColNames <- grep(colnames(dt), pattern = "cover", value = TRUE)
  newCoverColNames <- gsub("cover\\.", "", coverColNames)
  setnames(dt, old = coverColNames, new = newCoverColNames)
  message(crayon::blue("Step 4: Create initial cohortData object, with no pixelGroups yet"))

  cohortData <- data.table::melt(dt,
                                 value.name = "cover",
                                 measure.vars = newCoverColNames,
                                 variable.name = "speciesCode")
  cohortData[, coverOrig := cover]

  message(crayon::blue("Step 5: POSSIBLE ALERT -- assume deciduous cover is 1/2 the conversion to biomass as conifer"))
  cohortData[speciesCode == "Popu_sp", cover := asInteger(cover / 2)] # CRAZY TODO -- DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  cohortData[ , cover := asInteger(cover/(sum(cover) + 0.0001) * 100L), by = "pixelIndex"]

  # Biomass -- by cohort
  message(crayon::blue("Step 6: Divide total biomass of each pixel by the relative cover of the cohorts"))
  cohortData[ , biomass := asInteger(mean(totalBiomass) * cover / 100), by = "pixelIndex"] # /100 because cover is percent
  message(crayon::blue("Step 7: Round biomass to nearest P(sim)$pixelGroupBiomassClass"))
  cohortData[ , biomass := asInteger(ceiling(biomass / P(sim)$pixelGroupBiomassClass) *
                                       P(sim)$pixelGroupBiomassClass)] # /100 because cover is percent

  # https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
  # cohortData[ , coverProp := (cover/100 * (NROW(cohortData) - 1) + 0.5) / NROW(cohortData)]


  dtShort <- cohortData[, list(coverNum = .N,
                               coverPres = sum(cover > 0)),
                        by = c("ecoregionCode", "speciesCode", "lcc")]

  message(crayon::blue("Step 8: Estimaing Species Establishment Probability using P(sim)$coverQuotedFormula, which is\n",
                       format(P(sim)$coverQuotedFormula)))
  coverModel <- function(form, .specialData) {

    coverMod1 <- lme4::glmer(eval(form),
                             data = .specialData,
                             family = binomial)
    # MuMIn package MUST have the object, by name, in the global environment
    assign(".specialData", .specialData, envir = .GlobalEnv)
    on.exit(rm(list = ".specialData", envir = .GlobalEnv))
    return(list(mod = coverMod1, pred = fitted(coverMod1),
                rsq = MuMIn::r.squaredGLMM(coverMod1)))
  }
  #coverQuotedFormula <- quote(cbind(coverPres, coverNum) ~ speciesCode + (1 | ecoregionCode : lcc))
  outCover <- Cache(coverModel, form = P(sim)$coverQuotedFormula, .specialData = dtShort)
  message(crayon::blue("  The rsquared is: "))
  print(outCover$rsq)


  # For biomass
  cohortData <- cohortData[biomass > 0]
  message(crayon::blue("Step 9: rm pixels with no biomass, leaving", length(unique(cohortData$pixelIndex)), "pixels"))

  # For Cache -- doesn't need to cache all columns in the data.table -- only the ones in the model
  dtLongForLMER <- cohortData[, .(biomass, logAge, speciesCode, ecoregionCode, lcc, cover)]

  # biomassQuotedFormula <- quote(biomass ~ logAge * speciesCode + (speciesCode | ecoregionCode : lcc) + cover)
  message("Step 10: Estimaing maxBiomass with P(sim)$biomassQuotedFormula, which is:",
          paste0(format(P(sim)$biomassQuotedFormula, appendLF = FALSE), collapse = ""))
  biomassModel <- function(form, .specialData) {
    ba4 <- lmer(
      formula = form,
      data = .specialData)
    assign(".specialData", .specialData, envir = .GlobalEnv)
    on.exit(rm(list = ".specialData", envir = .GlobalEnv))
    # MuMIn package MUST have the object, by name, in the global environment
    list(mod = ba4, pred = fitted(ba4),
         rsq = MuMIn::r.squaredGLMM(ba4))
  }
  system.time(outBiomass <- Cache(biomassModel, form = P(sim)$biomassQuotedFormula, .specialData = dtLongForLMER, showSimilar = TRUE))

  # Add SEP to dtShort
  dtShort[, SEP := outCover$pred]

  # Create initial communities, i.e., pixelGroups
  columnsForPG <- c("ecoregionCode", "speciesCode", "age", "biomass")
  cd <- cohortData[,c("pixelIndex", columnsForPG), with = FALSE]
  cohortData[, pixelGroup :=
               Cache(generatePixelGroups, cd, maxPixelGroup = 0,
                     columns = columnsForPG)]
  message(crayon::blue("Step 11: Create pixelGroups based on: ", paste(columnsForPG, collapse = ", "),
                       " Resulted in", length(unique(cohortData$pixelGroup)), "uniqe pixelGroup values"))

  if (FALSE) {

    library(glmmTMB)
    a <- glmmTMB::glmmTMB(coverProp ~ logAge * speciesCode + (1 | ecoregionCode) + (1 | lcc),
                          data = cohortData,
                          family=beta_family,
                          verbose = TRUE)

    library(gamlss)
    cohortData[ , coverProp := (cover/100)]
    dtLong2 <- setDF(cohortData)
    assign("cohortData", cohortData, envir = .GlobalEnv)
    form <- formula(coverProp ~ logAge * speciesCode + re(random=~1|ecoregionCode) + re(random=~1|lcc))
    m.gamlss <- gamlss(formula = form,
                       # sigma.formula = ~logAge * speciesCode,
                       family=BEOI,
                       data=cohortData)
    setDT(cohortData)
    cohortData[ , fittedCover := fitted(m.gamlss, type = "response")]


    library(googledrive)

    folder <- drive_mkdir("upload-into-me-article-demo")
    #> Auto-refreshing stale OAuth token.
    #> Folder created:
    #>   * upload-into-me-article-demo
    files <- map(local_files, drive_upload, path = folder, verbose = FALSE)

    #b <- lme4::lmer(biomass ~ logAge + (1 | ecoregionCode) + (1 | lcc) + speciesCode + cover, data = cohortData)
    b <- lme4::lmer(biomass ~ logAge + (1 | ecoregionCode) + (1 | lcc) + speciesCode + cover, data = cohortData)
    b3 <- lme4::lmer(biomass ~ 0 + logAge + (1 | ecoregionCode) + (1 | lcc) + speciesCode + cover, data = cohortData)
    ba <- lme4::lmer(biomass ~ logAge + (1 | ecoregionCode : lcc) + speciesCode + cover, data = cohortData)
    ba2 <- lme4::lmer(biomass ~ logAge * speciesCode + (1 | ecoregionCode : lcc) + cover, data = cohortData)
    ba3 <- lme4::lmer(biomass ~ logAge * speciesCode + (logAge | ecoregionCode : lcc) + cover, data = cohortData, verbose = 100)
    ba5 <- lme4::lmer(biomass ~ (logAge * speciesCode | ecoregionCode : lcc) + cover, data = cohortData,
                      verbose = 100)
    ba6 <- lme4::lmer(biomass ~ logAge * speciesCode + cover + (logAge * speciesCode | ecoregionCode : lcc),
                      data = cohortData,
                      verbose = 100)

    # dtLong2 <- copy(cohortData)
    # dtLong3 <- copy(cohortData)
    # cohortData <- copy(dtLong2)
    library(nlme)
    data <- groupedData(y ~ t | UID, data=data) ## not strictly necessary
    initVals <- getInitial(y ~ SSlogis(t, Asym, xmid, scal), data = data)
    baseModel<- nlme(y ~ SSlogis(t, Asym, xmid, scal),
                     data = data,
                     fixed = list(Asym ~ 1, xmid ~ 1, scal ~ 1),
                     random = Asym + xmid + scal ~ 1|UID,
                     start = initVals
    )
    library(nlme)
    data <- cohortData
    plot(biomass ~ age, data = data)

    fit <- nls(biomass ~ SSlogis(age, Asym, xmid, scal), data=data)
    summary(fit)
    curve(predict(fit, newdata = data.frame(age=x)), add=TRUE)

    #data <- groupedData(biomass ~ age | ecoregionCode, data=data) ## not strictly necessary
    #initVals <- getInitial(biomass ~ SSlogis(age, Asym, xmid, scal), data = data)
    initVals <- coef(fit)
    initVals <- append(initVals, 30, 2)
    initVals <- append(initVals, 10, 1)

    baseModel<- nlme(biomass ~ SSlogis(age, Asym, xmid, scal),
                     data = data,
                     fixed = list(Asym ~ age, xmid ~ age, scal ~ 1),
                     random = Asym + xmid + scal ~ 1|ecoregionCode,
                     start = initVals,
                     control = lmeControl(msMaxIter = 150),
                     verbose = 100
    )

    nestedModel <- update(baseModel,
                          fixed=list(Asym ~ var1, xmid ~ 1, scal ~ 1),
                          start = c(fixef(baseModel)[1], 0, fixef(baseModel)[2], fixef(baseModel)[3]))


    ba1 <- lme4::glmer(biomass ~ logAge + (1 | ecoregionCode : lcc) + speciesCode + cover, data = cohortData,
                       family = gaussian("log"))
    b1 <- lme4::glmer(biomass ~ logAge + (1 | ecoregionCode) + (1 | lcc) + speciesCode + cover, data = cohortData,
                      family = gaussian(link = "log"))
    b1b <- lme4::glmer(biomass ~ logAge + (1 | ecoregionCode + lcc) + speciesCode + cover, data = cohortData,
                       family = gaussian(link = "log"))
    #b1 <- lme4::lmer(biomass ~ logAge + (speciesCode | ecoregionCode) + (speciesCode | lcc) + speciesCode + cover, data = cohortData)
    b2 <- lme4::lmer(biomass ~ logAge * speciesCode + (speciesCode | ecoregionCode) + (speciesCode | lcc) + cover, data = cohortData)
    b1a <- lme4::glmer(biomass ~ logAge * speciesCode + (speciesCode | ecoregionCode : lcc)  + cover, data = cohortData,
                       family = gaussian(link = "log"))

    MuMIn::r.squaredGLMM(b)
    MuMIn::r.squaredGLMM(b1)
    MuMIn::r.squaredGLMM(ba)
    MuMIn::r.squaredGLMM(ba2)
    MuMIn::r.squaredGLMM(ba3)
    #piecewiseSEM::rsquared(b1)
  }


  ######### CHANGE 34 and 35 values -- these are burns
  burnedLCC <- raster(sim$LCC2005);
  burnedLCC[] <- NA;
  burnedLCC[sim$LCC2005[] %in% 34:35] <- 1
  theBurnedCells <- which(burnedLCC[] == 1)
  iterations <- 1
  try(rm(list = "out3"), silent = TRUE)
  while (length(theBurnedCells) > 0) {
    iterations <- iterations + 1
    print(iterations)
    out <- spread2(sim$LCC2005, start = theBurnedCells, asRaster = FALSE,
                   iterations = iterations, allowOverlap = TRUE, spreadProb = 1)
    out[, lcc := sim$LCC2005[][pixels]]
    out[lcc %in% c(34:35, 37:39), lcc:=NA]
    out <- na.omit(out)
    out2 <- out[, list(newPossLCC = sample(lcc, 1)), by = "initialPixels"]
    if (!exists("out3")) {
      out3 <- out2
    } else {
      out3 <- rbindlist(list(out2, out3))
    }
    theBurnedCells <- theBurnedCells[!theBurnedCells %in% out2$initialPixels]
  }
  setnames(out3, "initialPixels", "pixelIndex")
  cohortData <- out3[cohortData, on = "pixelIndex"]
  cohortData[, lcc2 := as.integer(as.character(lcc))]
  cohortData[lcc2 %in% 34:35, lcc := newPossLCC]
  cohortData[, lcc := NULL]

  cohortData[, lcc3 := factor(lcc2)]
  if (getOption("LandR.assertions")) {
    assert1 <- all(as.integer(as.character(cohortData$lcc3)) == cohortData$lcc)
    if (!assert1)
      stop("lcc classes were mismanaged; contact developers")
  }

  cohortData[, lcc := lcc3]
  cohortData[, `:=`(lcc3 = NULL, lcc2 = NULL, newPossLCC = NULL,
                    logAge = NULL, cover = NULL, coverOrig = NULL,
                    #coverProp = NULL,
                    totalBiomass = NULL)]



  ################### Done
  joinOn <- c("ecoregionCode", "lcc", "speciesCode")
  speciesEcoregionTable <- unique(cohortData, by = joinOn)
  speciesEcoregionTable[, c("pixelIndex", "age", "biomass") := NULL]
  speciesEcoregionTable[lcc %in% unique(dtLongForLMER$lcc)]
  sim$species[, speciesCode := as.factor(species)]
  browser()
  speciesEcoregionTable <- sim$species[, .(speciesCode, longevity)][speciesEcoregionTable, on = "speciesCode"]

  # Set age to the age of longevity and cover to 100%
  speciesEcoregionTable[, `:=`(logAge = log(longevity), cover = 100)]

  speciesEcoregionTable[ , maxBiomass := asInteger(predict(outBiomass$mod,
                                                           newdata = speciesEcoregionTable,
                                                           type = "response"))]
  speciesEcoregionTable[maxBiomass < 0, maxBiomass := 0]
  speciesEcoregionTable[ , maxANPP := maxBiomass / 30]

  # Get SEP
  speciesEcoregionTable <- dtShort[, .(ecoregionCode, lcc, speciesCode, SEP)][speciesEcoregionTable, on = joinOn]

  speciesEcoregionTable[ , `:=`(logAge = NULL, cover = NULL, longevity = NULL, #pixelIndex = NULL,
                                lcc = NULL)]
  # speciesEcoregionTable[, lccChar := as.character(lcc)]
  # speciesEcoregionTable[, lccCode := paddedFloatToChar(as.integer(lcc),
  #                                                      padL = max(nchar(lccChar)),
  #                                                      padR = 0)]
  # speciesEcoregionTable[ , ecoregionCode := factor(paste(ecoregionCode, lccCode, sep = "_"))]#,
  #by = seq(1:NROW(speciesEcoregionTable))]

  if (!is.na(P(sim)$.plotInitialTime)) {
    uniqueSpeciesNames <- as.character(unique(speciesEcoregionTable$speciesCode))
    names(uniqueSpeciesNames) <- uniqueSpeciesNames
    speciesEcoregionTable2 <- copy(speciesEcoregionTable)
    speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionCode)]
    maxBiomass <- stack(lapply(uniqueSpeciesNames, function(sp) {
      rasterizeReduced(speciesEcoregionTable2[speciesCode == sp], ecoregionFiles$ecoregionMap,
                       "maxBiomass", "ecoregionInt")
    }))
    Plot(maxBiomass, legendRange = c(0, max(maxValue(maxBiomass))))
  }


  if (FALSE) {
    maxAge <- 1:100
    maxAgeDT <- data.table(logAge = rep(log(maxAge), length(unique(speciesEcoregionTable$speciesCode))),
                           speciesCode = rep(unique(speciesEcoregionTable$speciesCode), each = length(maxAge)))
    newdat <- maxAgeDT[speciesEcoregionTable, on = c("speciesCode"), allow.cartesian = TRUE]
    plot(maxAge, predict(bestModel, newdat[pixelIndex == 651155]))

    # Checks
    dd3 <- speciesEcoregionTable[, mean(maxBiomass), by = c("lcc")]
    setkey(dd3, lcc)
    # slope is negative meaning the higher lcc classes have lower maxBiomass -- TODO: should be a better test
    lm(V1 ~ as.numeric(lcc), data = dd3)



    d <- speciesEcoregionTable[, list(meanMaxBiomass = mean(maxBiomass)), by = c("ecoregionCode", "lcc", "speciesCode")]
    dd <- d[, list(mostAbundantSpecies = speciesCode[which.max(meanMaxBiomass)]), by = c("ecoregionCode", "lcc")]
    table(dd$mostAbundantSpecies)
    # plot(cohortData$biomass, cohortData$fittedBiomass)
    # coef(b1)
    maxAge <- 1:100
    newdat <- data.frame(ecoregionCode = "643", lcc = "1", logAge = log(maxAge), speciesCode = "Pice_gla", cover = 100)
    plot(maxAge, predict(bestModel, newdat))
  }

  if (FALSE) { # old way -- superceded by Eliot Dec 30, 2018 -- was giving bad values
    message("2: ecoregionProducer: ", Sys.time())
    # Note: this ecoregionMap is NOT the Canadian EcoRegion -- it is for LBMR, which uses "ecoregion"
    ecoregionMap <- Cache(postProcess, sim$ecoDistrict, studyArea = sim$studyArea, filename2 = NULL)
    ecoregionstatus <- data.table(active = "yes",
                                  ecoregion = 1:1031)
    ecoregionFiles <- Cache(ecoregionProducer,
                            ecoregionMaps = list(ecoregionMap, sim$LCC2005),
                            ecoregionName = "ECODISTRIC",
                            ecoregionActiveStatus = ecoregionstatus,
                            rasterToMatch = sim$rasterToMatch,
                            userTags = "stable")

    ecoRegionMap <- Cache(postProcess, sim$ecoRegion, studyArea = sim$studyArea, filename2 = NULL)
    ecoRegionFiles <- Cache(ecoregionProducer,
                            ecoregionMaps = list(ecoRegionMap, sim$LCC2005),
                            ecoregionName = "ECOREGION",
                            ecoregionActiveStatus = ecoregionstatus,
                            rasterToMatch = sim$rasterToMatch,
                            userTags = "stable")

    ############################################################
    # Create initialCommunities
    ############################################################
    message("3: Create initial communities map and data.table ", Sys.time())
    initialCommunities <- Cache(initialCommunityProducer,
                                speciesLayers = sim$speciesLayers,
                                # This cuts up speciesCode pct into x groups, and age into x-year groups
                                percentileGrps = 10,
                                ecoregionMap = ecoregionFiles$ecoregionMap,
                                standAgeMap = sim$standAgeMap,
                                userTags = "stable")
    # remove unneeded columns of initialCommunities based on sim$speciesCode
    # initialCommunities1 <- initialCommunities[, .(mapcode,
    #                                              description = NA,
    #                                              species,
    #                                              speciesPresence = as.integer(speciesPresence),
    #                                              age = as.integer(age1),
    #                                              pixelIndex = as.integer(pixelIndex))]

    ############################################################
    # Create initialCommunitiesMap
    ############################################################
    browser()
    initialCommunitiesMap <- raster(sim$speciesLayers[[1]])
    initialCommunitiesMap[] <- NA_integer_
    initialCommunitiesByPixel <- unique(initialCommunities, by = c("mapcode", "pixelIndex"))
    initialCommunitiesMap[initialCommunitiesByPixel$pixelIndex] <- as.integer(initialCommunitiesByPixel$mapcode)
    datatype <- if (max(initialCommunities$mapcode) > 64e3)
      "INT4U" else "INT2U"
    sim$initialCommunitiesMap <- Cache(raster::writeRaster,
                                       initialCommunitiesMap,
                                       filename = file.path(outputPath(sim), "initialCommunitiesMap.tif"),
                                       overwrite = TRUE,
                                       datatype = datatype)

    if (ncell(sim$rasterToMatch) > 3e6)  .gc()

    ############################################################
    # Create speciesEcoregion Table
    ############################################################
    message("4: Prepare speciesEcoregion Table -- running obtainMaxBandANPP", Sys.time())

    possibleEcoregionSrcs <- list(ecoDistrict = ecoregionFiles$ecoregionMap,
                                  ecoRegion = ecoRegionFiles$ecoregionMap)

    speciesEcoregionTable <- Cache(createSpeciesEcoregion,
                                   possibleEcoregionSrcs = possibleEcoregionSrcs,
                                   rasterToMatch = sim$rasterToMatch,
                                   speciesLayers = sim$speciesLayers,
                                   biomassMap = sim$biomassMap,
                                   minNumPixelsToEstMaxBiomass = P(sim)$minNumPixelsToEstMaxBiomass,
                                   quantileForMaxBiomass = P(sim)$quantileForMaxBiomass)
    if (ncell(sim$rasterToMatch) > 3e6) .gc()

    # This will create stacks that may be useful
    if (!is.na(P(sim)$.plotInitialTime)) {
      uniqueSpeciesNames <- as.character(unique(speciesEcoregionTable$species))
      names(uniqueSpeciesNames) <- uniqueSpeciesNames
      speciesEcoregionTable2 <- copy(speciesEcoregionTable)
      speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionCode)]
      maxBiomass <- stack(lapply(uniqueSpeciesNames, function(sp) {
        rasterizeReduced(speciesEcoregionTable2[species == sp], ecoregionFiles$ecoregionMap,
                         "maxBiomass", "ecoregionInt")
      }))
      Plot(maxBiomass, legendRange = c(0, max(maxValue(maxBiomass))))
    }

    ################################################################
    # SEP
    ################################################################
    message("5: Derive Species Establishment Probability (SEP) from sim$speciesLayers: ", Sys.time())
    sepTable <- Cache(obtainSEP, possibleEcoregionSrcs = possibleEcoregionSrcs,
                      speciesEcoregionTable = speciesEcoregionTable,
                      speciesLayers = sim$speciesLayers)

    sim$speciesEstablishmentProbMap <- createSEPStack(speciesEcoregionTable = sepTable,
                                                      ecoregionMap = ecoregionFiles$ecoregionMap,
                                                      destinationPath = Paths$inputPath)

    speciesEcoregionTable[, active := "yes"]#ecoregionFiles$ecoregion, on = c(ecoregionCode = "ecoregion")]

    if ("maxBiomass" %in% names(sepTable))
      setnames(sepTable, c("maxBiomass"), c("maxB"))
    speciesEcoregion <- speciesEcoregionTable[, .(ecoregionCode, species, active, maxB = maxBiomass, maxANPP)][
      sepTable, on = c("ecoregionCode", "species")]
    if ("SEP" %in% names(speciesEcoregion))
      setnames(speciesEcoregion, c("SEP"), c("establishprob"))
    set(speciesEcoregion, NULL, "year", 0L)


    lvs <- as.data.table(raster::levels(sim$ecoregionMap)[[1]] )

    sim$notEnoughDataMaxBiomass <- speciesEcoregion[is.na(maxB)]

    if (!is.na(P(sim)$.plotInitialTime)) {
      uniqueSpeciesNames <- as.character(unique(speciesEcoregion$species))
      names(uniqueSpeciesNames) <- uniqueSpeciesNames
      sim$notEnoughDataMaxBiomass[, ecoregionInt := as.integer(ecoregionCode)] # need for plotting

      noMaxBiomass <- stack(lapply(uniqueSpeciesNames, function(sp) {
        r <- rasterizeReduced(sim$notEnoughDataMaxBiomass[species == sp], ecoregionFiles$ecoregionMap,
                              "year", "ecoregionInt")
        r[!is.na(r)] <- 1
        r
      }))
      Plot(noMaxBiomass, cols = "blue")

    }


    # Clean it up for return to LBMR
    speciesEcoregion <- speciesEcoregion[, .(year, ecoregion = ecoregionCode,
                                             species, maxB, maxANPP, establishprob)]
    speciesEcoregion[is.na(maxB)]

    sim$speciesEcoregion <- speciesEcoregion
    sim$ecoregion <- ecoregionFiles$ecoregion
    sim$ecoregionMap <- ecoregionFiles$ecoregionMap

  }
  if (ncell(sim$rasterToMatch) > 3e6) .gc()



  browser()
  ## filter communities to species that have traits
  sim$initialCommunities <-  cohortData



  sim$minRelativeB <- data.frame(ecoregion = sim$ecoregion[active == "yes",]$ecoregion,
                                 X1 = 0.2, X2 = 0.4, X3 = 0.5,
                                 X4 = 0.7, X5 = 0.9)

  message("Done Boreal_LBMRDataPrep: ", Sys.time())

  sim$speciesLayers <- lapply(seq(numLayers(sim$speciesLayers)), function(x) {
    writeRaster(sim$speciesLayers[[x]],
                file.path(outputPath(sim), paste0(names(sim$speciesLayers)[x], ".tif")),
                datatype = "INT2U", overwrite = TRUE)
  }) %>% raster::stack()

  return(invisible(sim))
}

Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

## see other helper functions in R/ subdirectory

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "function:.inputObjects")
  dPath <- asPath(getOption("reproducible.destinationPath", dataPath(sim)), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # 1. test if all input objects are already present (e.g., from inputs, objects or another module)
  a <- depends(sim)
  whThisMod <- which(unlist(lapply(a@dependencies, function(x) x@name)) == "Boreal_LBMRDataPrep")
  objNames <- a@dependencies[[whThisMod]]@inputObjects$objectName
  objExists <- !unlist(lapply(objNames, function(x) is.null(sim[[x]])))
  names(objExists) <- objNames

  # Filenames
  ecoregionFilename <-   file.path(dPath, "ecoregions.shp")
  ecodistrictFilename <- file.path(dPath, "ecodistricts.shp")
  ecozoneFilename <-   file.path(dPath, "ecozones.shp")
  biomassMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.tif")
  lcc2005Filename <- file.path(dPath, "LCC2005_V1_4a.tif")
  standAgeMapFilename <- file.path(dPath, "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.tif")

  # Also extract
  fexts <- c("dbf", "prj", "sbn", "sbx", "shx")
  ecoregionAE <- basename(paste0(tools::file_path_sans_ext(ecoregionFilename), ".", fexts))
  ecodistrictAE <- basename(paste0(tools::file_path_sans_ext(ecodistrictFilename), ".", fexts))
  ecozoneAE <- basename(paste0(tools::file_path_sans_ext(ecozoneFilename), ".", fexts))

  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using a polygon in southwestern Alberta, Canada,")

    sim$studyArea <- randomStudyArea(seed = 1234)
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
    sim$studyAreaLarge <- sim$studyArea
  }

  needRTM <- FALSE
  if (is.null(sim$rasterToMatch)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use biomassMap")
    } else {
      stop("rasterToMatch is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("biomassMap", sim) || needRTM) {
    sim$biomassMap <- Cache(prepInputs,
                            targetFile = asPath(basename(biomassMapFilename)),
                            archive = asPath(c("kNN-StructureBiomass.tar",
                                               "NFI_MODIS250m_kNN_Structure_Biomass_TotalLiveAboveGround_v0.zip")),
                            url = extractURL("biomassMap"),
                            destinationPath = dPath,
                            studyArea = sim$studyArea,
                            rasterToMatch = sim$rasterToMatch, ## TODO: biomass map needs rasterToMatch but it _is_ the rasterToMatch!!
                            maskWithRTM = TRUE,
                            useSAcrs = TRUE,
                            method = "bilinear",
                            datatype = "INT2U",
                            filename2 = TRUE, overwrite = TRUE,
                            userTags = c(cacheTags, "stable"))
  }

  if (needRTM) {
    # if we need rasterToMatch, that means a) we don't have it, but b) we will have biomassMap
    sim$rasterToMatch <- sim$biomassMap
    message("  Rasterizing the studyArea polygon map")
    if (!is(sim$studyArea, "SpatialPolygonsDataFrame")) {
      dfData <- if (is.null(rownames(sim$studyArea))) {
        polyID <- sapply(slot(sim$studyArea, "polygons"), function(x) slot(x, "ID"))
        data.frame("field" = as.character(seq_along(length(sim$studyArea))), row.names = polyID)
      } else {
        polyID <- sapply(slot(sim$studyArea, "polygons"), function(x) slot(x, "ID"))
        data.frame("field" = rownames(sim$studyArea), row.names = polyID)
      }
      sim$studyArea <- SpatialPolygonsDataFrame(sim$studyArea, data = dfData)
    }
    #TODO: review whether this is necessary (or will break LandWeb if removed) see Git Issue #22
    # layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
    LTHxC <- grep("(LTH.+C)", names(sim$studyArea), value = TRUE)
    fieldName <- if (length(LTHxC)) {
      LTHxC
    } else {
      if (length(names(sim$studyArea)) > 1) {
        ## study region may be a simple polygon
        names(sim$studyArea)[1]
      } else NULL
    }

    sim$rasterToMatch <- crop(fasterizeFromSp(sim$studyArea, sim$rasterToMatch, fieldName),
                              sim$studyArea)
    sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                               filename = file.path(dataPath(sim), "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE)
  }

  # LCC2005
  if (!suppliedElsewhere("LCC2005", sim)) {
    sim$LCC2005 <- Cache(prepInputs,
                         targetFile = lcc2005Filename,
                         archive = asPath("LandCoverOfCanada2005_V1_4.zip"),
                         url = extractURL("LCC2005"),
                         destinationPath = dPath,
                         studyArea = sim$studyArea,
                         rasterToMatch = sim$rasterToMatch,
                         method = "bilinear",
                         datatype = "INT2U",
                         filename2 = TRUE, overwrite = TRUE,
                         userTags = currentModule(sim))

    projection(sim$LCC2005) <- projection(sim$rasterToMatch)
  }
  if (!suppliedElsewhere("ecoDistrict", sim)) {
    sim$ecoDistrict <- Cache(prepInputs,
                             targetFile = asPath(ecodistrictFilename),
                             archive = asPath("ecodistrict_shp.zip"),
                             url = extractURL("ecoDistrict"),
                             alsoExtract = ecodistrictAE,
                             destinationPath = dPath,
                             studyArea = sim$studyAreaLarge,
                             overwrite = TRUE,
                             useSAcrs = TRUE, # this is required to make ecoZone be in CRS of studyArea
                             fun = "raster::shapefile",
                             filename2 = TRUE,
                             userTags = cacheTags)
  }

  if (!suppliedElsewhere("ecoRegion", sim)) {
    sim$ecoRegion <- Cache(prepInputs,
                           targetFile = asPath(ecoregionFilename),
                           archive = asPath("ecoregion_shp.zip"),
                           alsoExtract = ecoregionAE,
                           url = extractURL("ecoRegion"),
                           destinationPath = dPath,
                           studyArea = sim$studyAreaLarge,
                           overwrite = TRUE,
                           useSAcrs = TRUE, # this is required to make ecoZone be in CRS of studyArea
                           fun = "raster::shapefile",
                           filename2 = TRUE,
                           userTags = cacheTags)
  }

  if (!suppliedElsewhere("ecoZone", sim)) {
    sim$ecoZone <- Cache(prepInputs, #notOlderThan = Sys.time(),
                         targetFile = asPath(ecozoneFilename),
                         archive = asPath("ecozone_shp.zip"),
                         url = extractURL("ecoZone"),
                         alsoExtract = ecozoneAE,
                         destinationPath = dPath,
                         studyArea = sim$studyAreaLarge,
                         overwrite = TRUE,
                         useSAcrs = TRUE, # this is required to make ecoZone be in CRS of studyArea
                         fun = "raster::shapefile",
                         filename2 = TRUE,
                         userTags = cacheTags)
  }

  # stand age map
  if (!suppliedElsewhere("standAgeMap", sim)) {
    sim$standAgeMap <- Cache(prepInputs,
                             targetFile = basename(standAgeMapFilename), ## TODO: undefined filename
                             archive = asPath(c("kNN-StructureStandVolume.tar",
                                                "NFI_MODIS250m_kNN_Structure_Stand_Age_v0.zip")),
                             destinationPath = dPath,
                             url = extractURL("standAgeMap"),
                             fun = "raster::raster",
                             studyArea = sim$studyAreaLarge,
                             rasterToMatch = sim$rasterToMatch,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = TRUE, overwrite = TRUE,
                             userTags = c("stable", currentModule(sim)))
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)

    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]
  }

  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    speciesLayersList <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppEquivCol = P(sim)$sppEquivCol,
                               # thresh = 10,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"))

    #options(opts)
    speciesLayersList$speciesLayers <- writeRaster(speciesLayersList$speciesLayers,
                                                   file.path(outputPath(sim), "speciesLayers.grd"),
                                                   overwrite = TRUE)
    sim$speciesLayers <- speciesLayersList$speciesLayers
  }

  # 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath, cacheTags)
  }

  if (!suppliedElsewhere("sufficientLight", sim)) {
    sim$sufficientLight <- data.frame(speciesshadetolerance = 1:5,
                                      X0 = 1,
                                      X1 = c(0.5, rep(1, 4)),
                                      X2 = c(0, 0.5, rep(1, 3)),
                                      X3 = c(rep(0, 2), 0.5, rep(1, 2)),
                                      X4 = c(rep(0, 3), 0.5, 1),
                                      X5 = c(rep(0, 4), 1))
  }

  return(invisible(sim))
}
