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
  version = list(SpaDES.core = "0.2.3.9009", Boreal_LBMRDataPrep = numeric_version("1.3.4")),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Boreal_LBMRDataPrep.Rmd"),
  reqdPkgs = list("crayon", "data.table", "dplyr", "fasterize", "gdalUtils", "lme4", "plyr",
                  "raster", "rgeos", "sp",
                  "achubaty/amc@development",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    defineParameter("biomassModel", "call",
                    quote(lme4::lmer(B ~ logAge * speciesCode + cover * speciesCode + (logAge + cover + speciesCode | ecoregionGroup))),
                    NA, NA,
                    paste("Model and formula for estimating biomass (B) from ecoregionGroup (currently ecoDistrict * LandCoverClass),",
                          "speciesCode, logAge (gives a downward curving relationship), and cover. Defaults to a LMEM, which",
                          "can be slow if dealing with very large datasets (e.g. 36 000 points take 20min).",
                          "For faster fitting try P(sim)$subsetDataBiomassModel == TRUE, or",
                          "quote(RcppArmadillo::fastLm(formula = B ~ logAge * speciesCode * ecoregionGroup + cover",
                          "* speciesCode * ecoregionGroup)). A custom model call can also be provided,",
                          "as long as the 'data' argument is NOT included")),
    defineParameter("coverModel", "call",
                    quote(lme4::glmer(cbind(coverPres, coverNum) ~ speciesCode + (1 | ecoregionGroup),
                                      family = binomial)),
                    NA, NA,
                    paste("Model and formula used for estimating cover from ecoregion and speciesCode",
                          "and potentially others. Defaults to a GLMEM if there are > 1 grouping levels.",
                          "A custom model call can also be provided, as long as the 'data' argument is NOT included")),
    defineParameter("forestedLCCClasses", "numeric", c(1:15, 20, 32, 34:35), 0, 39,
                    "The classes in the rstLCC layer that are 'treed' and will therefore be run in LBMR.
                    Defaults to forested classes in LCC2005 map."),
    defineParameter("growthCurveDecid", "numeric", 0, 0, 1,
                    "growth curve shape for deciduous species (i.e., aspen). LANDIS-II uses 0 for aspen."),
    defineParameter("growthCurveNonDecid", "numeric", 1, 0, 1,
                    "growth curve shape for non-deciduous species. LANDIS-II uses 1 for spruce, 0 for other non deciduous."),
    defineParameter("LCCClassesToReplaceNN", "numeric", 34:35, NA, NA,
                    paste("This will replace these classes on the landscape with the closest forest class P(sim)$forestedLCCClasses.",
                          "If the user is using the default 2005 data product for rstLCC, then users may wish to",
                          "include 36 (cities -- if running a historic range of variation project), and 34:35 (burns)",
                          "Since this is about estimating parameters for growth, it doesn't make any sense to have",
                          "unique estimates for transient classes in most cases")),
    defineParameter("mortalityShapeDecid", "numeric", 25, 12, 27,
                    "mortality curve shape for deciduous species (i.e., aspen). LANDIS-II uses 25 by default."),
    defineParameter("mortalityShapeNonDecid", "numeric", 15, 12, 27,
                    "mortality curve shape for non-deciduous species. LANDIS-II uses 15 by default."),
    defineParameter("omitNonTreedPixels", "logical", TRUE, FALSE, TRUE,
                    "Should this module use only treed pixels, as identified by P(sim)$forestedLCCClasses?"),
    defineParameter("pixelGroupAgeClass", "numeric", params(sim)$Boreal_LBMRDataPrep$successionTimestep, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same"),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same"),
    defineParameter("runName", "character", "", NA, NA,
                    paste("A description for run.",
                          "This will form the basis of cache path and output path, and affect dispersal parameterization.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("subsetDataAgeModel", "numeric", NULL, NA, NA,
                    "the number of samples to use when subsampling the biomass data model; if TRUE, uses 50"),
    defineParameter("subsetDataBiomassModel", "numeric", NULL, NA, NA,
                    "the number of samples to use when subsampling the biomass data model; if TRUE, uses 50"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter("useCloudCacheForStats", "logical", TRUE, NA, NA,
                    paste("Some of the statistical models take long (at least 30 minutes, likely longer).",
                          "If this is TRUE, then it will try to get previous cached runs from googledrive")),
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
    expectsInput("cloudFolderID", "character",
                 "The google drive location where cloudCache will store large statistical objects"),
    expectsInput("columnsForPixelGroups", "character",
                 "The names of the columns in cohortData that define unique pixelGroups. Default is c('ecoregionGroup', 'speciesCode', 'age', 'B') "),
    expectsInput("ecoDistrict", "SpatialPolygonsDataFrame",
                 desc = "ecodistricts in study area, default is Canada national ecodistricts",
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("rstLCC", "RasterLayer",
                 desc = paste("A land classification map in study area. It must be 'corrected', in the sense that:\n",
                              "1) Every class must not conflict with any other map in this module\n",
                              "    (e.g., speciesLayers should not have data in LCC classes that are non-treed);\n",
                              "2) It can have treed and non-treed classes. The non-treed will be removed within this\n",
                              "    module if P(sim)$omitNonTreedPixels is TRUE;\n",
                              "3) It can have transient pixels, such as 'young fire'. These will be converted to a\n",
                              "    the nearest non-transient class, probabilistically if there is more than 1 nearest\n",
                              "    neighbour class, based on P(sim)$LCCClassesToReplaceNN.\n",
                              "The default layer used, if not supplied, is Canada national land classification in 2005"),
                 sourceURL = "https://drive.google.com/file/d/1g9jr0VrQxqxGjZ4ckF6ZkSMP-zuYzHQC/view?usp=sharing"),
    expectsInput("rasterToMatch", "RasterLayer",
                 #desc = "this raster contains two pieces of information: Full study area with fire return interval attribute",
                 desc = "DESCRIPTION NEEDED",
                 sourceURL = NA),
    expectsInput("speciesLayers", "RasterStack",
                 desc = "cover percentage raster layers by species in Canada species map",
                 sourceURL = "http://tree.pfc.forestry.ca/kNN-Species.tar"),
    expectsInput("speciesTable", "data.table",
                 desc = "species attributes table, default is from Dominic Cyr and Yan Boulanger's project",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species",
                 sourceURL = ""),
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
    createsOutput("ecoregion", "data.table",
                  desc = "ecoregion look up table"),
    createsOutput("ecoregionMap", "RasterLayer",
                  desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table"),
    createsOutput("cohortData", "data.table",
                  desc = "initial community table"),
    createsOutput("pixelGroupMap", "RasterLayer",
                  desc = "initial community map that has mapcodes match initial community table"),
    createsOutput("minRelativeB", "data.frame",
                  desc = "define the cut points to classify stand shadeness"),
    createsOutput("notEnoughDataMaxBiomass", "data.table",
                  desc = "The collection of ecoregion-species combinations that don't have values for maxB or maxANPP"),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
    createsOutput("speciesEcoregion", "data.table",
                  desc = "define the maxANPP, maxB and establishprob change with both ecoregion and simulation time"),
    createsOutput("studyArea", "", desc = ""),
    # createsOutput("speciesEstablishmentProbMap", "RasterStack",
    #               paste("Species establishment probability as a map, ",
    #                     "by species. This is written to disk to save RAM space")),
    createsOutput("useCache", "logic",
                  desc = "define which the caching for spinup simulation should be used, default is TRUE")
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Boreal_LBMRDataPrep <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    # names(sim$speciesLayers) <- equivalentName(names(sim$speciesLayers), sim$sppEquiv, "Latin_full")
    sim <- createLBMRInputs(sim)

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

createLBMRInputs <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  if (is.null(P(sim)$pixelGroupAgeClass))
    params(sim)[[currentModule(sim)]]$pixelGroupAgeClass <- P(sim)$successionTimestep

  message(blue("Starting to createLBMRInputs in Boreal_LBMRDataPrep: ", Sys.time()))
  cPath <- cachePath(sim)
  sim$ecoDistrict <- spTransform(sim$ecoDistrict, crs(sim$speciesLayers))

  sim$standAgeMap <- round(sim$standAgeMap / 20, 0) * 20 # use 20-year bins (#103)
  sim$standAgeMap[] <- asInteger(sim$standAgeMap[])

  ################################################################
  ## species traits inputs
  ################################################################
  message(blue("Prepare 'species' table, i.e., species level traits", Sys.time()))
  sim$species <- prepSpeciesTable(speciesTable = sim$speciesTable,
                                  speciesLayers = sim$speciesLayers,
                                  sppEquiv = sim$sppEquiv[get(P(sim)$sppEquivCol) %in%
                                                            names(sim$speciesLayers)],
                                  sppEquivCol = P(sim)$sppEquivCol)

  ### override species table values ##############################

  if (grepl("aspenDispersal", P(sim)$runName)) {
    ## seed dispersal (see LandWeb#96, LandWeb#112)
    sim$species[species == "Abie_sp", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 25, 160
    sim$species[species == "Pice_gla", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 100, 303
    sim$species[species == "Pice_mar", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 80, 200
    sim$species[species == "Pinu_ban", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 30, 100
    sim$species[species == "Pinu_con", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 30, 100
    sim$species[species == "Pinu_sp", `:=`(seeddistance_eff = 0, seeddistance_max = 125)] # defaults 30, 100
    sim$species[species == "Popu_sp", `:=`(seeddistance_eff = 100, seeddistance_max = 235)] # defaults 200, 5000
  } else {
    ## seed dispersal (see LandWeb#96, LandWeb#112)
    sim$species[species == "Abie_sp", `:=`(seeddistance_eff = 100, seeddistance_max = 750)] # defaults 25, 160
    sim$species[species == "Pice_gla", `:=`(seeddistance_eff = 300, seeddistance_max = 750)] # defaults 100, 303
    sim$species[species == "Pice_mar", `:=`(seeddistance_eff = 300, seeddistance_max = 750)] # defaults 80, 200
    sim$species[species == "Pinu_ban", `:=`(seeddistance_eff = 300, seeddistance_max = 750)] # defaults 30, 100
    sim$species[species == "Pinu_con", `:=`(seeddistance_eff = 300, seeddistance_max = 750)] # defaults 30, 100
    sim$species[species == "Pinu_sp", `:=`(seeddistance_eff = 300, seeddistance_max = 750)] # defaults 30, 100
    #sim$species[species == "Popu_sp", `:=`(seeddistance_eff = 300, seeddistance_max = 3000)] # defaults 200, 500
  }

  if (grepl("noDispersal|aspenDispersal", P(sim)$runName)) {
    sim$species[, postfireregen := "none"]
  }

  ## resprouting (only aspen resprouts)
  sim$species[species == "Popu_sp", resproutage_min := 25] # default 10
  #speciesTable[species == "Popu_sp", resproutprob := 0.1] # default 0.5

  ## growth curves:
  #   Biomass Succession User Guide p17, 0 is faster growth, 1 was the prev assumption
  sim$species[species == "Abie_sp", growthcurve := P(sim)$growthCurveNonDecid] # original default 0
  sim$species[species == "Pice_gla", growthcurve := P(sim)$growthCurveNonDecid] # original default 1
  sim$species[species == "Pice_mar", growthcurve := P(sim)$growthCurveNonDecid] # original default 1
  sim$species[species == "Pinu_ban", growthcurve := P(sim)$growthCurveNonDecid] # original default 0
  sim$species[species == "Pinu_con", growthcurve := P(sim)$growthCurveNonDecid] # original default 0
  sim$species[species == "Pinu_sp", growthcurve := P(sim)$growthCurveNonDecid] # original default 0
  sim$species[species == "Popu_sp", growthcurve := P(sim)$growthCurveDecid] # original default 0

  ## mortality
  sim$species[species == "Abie_sp", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Pice_gla", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Pice_mar", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Pinu_ban", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Pinu_con", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Pinu_sp", mortalityshape := P(sim)$mortalityShapeNonDecid] # default 15
  sim$species[species == "Popu_sp", mortalityshape := P(sim)$mortalityShapeDecid] # default 25

  if (grepl("aspen80", P(sim)$runName)) {
    sim$species[species == "Popu_sp", longevity := 80] # default 150
  }

  if (getOption("LandR.verbose") > 0) {
    message("Adjusting species-level traits, part 2, for LandWeb")
    print(sim$species)
  }

  ################################################################
  ## initialEcoregionMap
  ################################################################
  if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  sim$ecoDistrict <- fixErrors(sim$ecoDistrict)

  ecoregionMap <- Cache(postProcess, sim$ecoDistrict, studyArea = sim$studyArea, filename2 = NULL)
  rstEcoregionMap <- fasterize::fasterize(sf::st_as_sf(ecoregionMap), raster = sim$rasterToMatch,
                                          field = "ECODISTRIC")
  ecoregionstatus <- data.table(active = "yes", ecoregion = 1:1031)
  rstLCCAdj <- sim$rstLCC

  ## Clean pixels for veg. succession model
  ## remove pixes with no spp data
  pixelsToRm <- is.na(sim$speciesLayers[[1]][])

  ## remove non-forested if asked by user
  if (P(sim)$omitNonTreedPixels) {
    if (is.null(P(sim)$forestedLCCClasses))
      stop("No P(sim)$forestedLCCClasses provided, but P(sim)$omitNonTreedPixels is TRUE.
           \nPlease provide a vector of forested classes in P(sim)$forestedLCCClasses")
    lccPixelsRemoveTF <- !(sim$rstLCC[] %in% P(sim)$forestedLCCClasses)
    pixelsToRm <- lccPixelsRemoveTF | pixelsToRm
  }

  rstLCCAdj[pixelsToRm] <- NA
  rstEcoregionMap[pixelsToRm] <- NA

  ## TODO: clean up - not the most effient function (maybe contains redundancies). Producing a non-used object
  message(blue("Make initial ecoregionGroups ", Sys.time()))
  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(rstEcoregionMap, rstLCCAdj),
                          ecoregionName = "ECODISTRIC",
                          ecoregionActiveStatus = ecoregionstatus,
                          rasterToMatch = sim$rasterToMatch,
                          userTags = "stable")

  ################################################################
  ## put together pixelTable object
  ################################################################
  #  Round age to pixelGroupAgeClass
  pixelTable <- makePixelTable(speciesLayers = sim$speciesLayers, species = sim$species,
                               standAgeMap = standAgeMap, ecoregionFiles = ecoregionFiles,
                               biomassMap = biomassMap, rasterToMatch = sim$rasterToMatch,
                               LCC2005 = LCC2005, pixelGroupAgeClass = P(sim)$pixelGroupAgeClass)
  # message(blue("Round age to nearest P(sim)$pixelGroupAgeClass, which is",
  #              P(sim)$pixelGroupAgeClass))
  # coverMatrix <- matrix(asInteger(sim$speciesLayers[]),
  #                       ncol = length(names(sim$speciesLayers)))
  # colnames(coverMatrix) <- names(sim$speciesLayers)
  #
  # pixelTable <- data.table(age = asInteger(ceiling(asInteger(sim$standAgeMap[]) /
  #                                                    P(sim)$pixelGroupAgeClass) *
  #                                            P(sim)$pixelGroupAgeClass),
  #                          logAge = log(sim$standAgeMap[]),
  #                          initialEcoregionCode = factor(factorValues2(ecoregionFiles$ecoregionMap,
  #                                                                      ecoregionFiles$ecoregionMap[],
  #                                                                      att = 5)),
  #                          totalBiomass = asInteger(sim$biomassMap[]) * 100, # change units
  #                          cover = coverMatrix,
  #                          pixelIndex = seq(ncell(sim$standAgeMap)),
  #                          lcc = rstLCCAdj[],
  #                          rasterToMatch = sim$rasterToMatch[]
  # )
  #
  # coverColNames <- paste0("cover.", sim$species$species)
  # pixelTable1 <- na.omit(pixelTable, cols = c("rasterToMatch"))
  # pixelTable2 <- na.omit(pixelTable, cols = c("rasterToMatch", "initialEcoregionCode"))
  # pixelTable <- na.omit(pixelTable2, cols = c(coverColNames))
  #
  # if (NROW(pixelTable1) != NROW(pixelTable))
  #   warning("Setting pixels to NA where there is NA in sim$speciesLayers'. Vegetation succession",
  #           "\n  parameters will only be calculated where there is data for species cover.",
  #           "\n  Check if these pixels should also be excluded in other layers,",
  #           "\n  as this may affect other modules.")
  # if (NROW(pixelTable2) != NROW(pixelTable))
  #   warning("Setting pixels to NA where there is NA in sim$ecoDistrict")
  #
  # message(blue("rm NAs, leaving", magenta(NROW(pixelTable)), "pixels with data"))
  # message(blue("This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"))
  # print(summary(pixelTable))

  #######################################################
  # Make the initial pixelCohortData table
  #######################################################
  pixelCohortData <- Cache(makeAndCleanInitialCohortData, pixelTable,
                           sppColumns = coverColNames,
                           pixelGroupBiomassClass = P(sim)$pixelGroupBiomassClass,
                           doSubset = P(sim)$subsetDataAgeModel)

  #######################################################
  # replace 34 and 35 and 36 values -- burns and cities -- to a neighbour class *that exists*
  #######################################################
  uwc <- P(sim)$LCCClassesToReplaceNN

  message("Replace ", paste(uwc, collapse = ", "),
          " values -- ", "burns"[any(uwc %in% 34:35)], " and cities"[any(uwc %in% 36)],
          " -- to a neighbour class *that exists*")

  rmZeroBiomassQuote <- quote(B > 0)
  # availableCombinations <- unique(pixelCohortData[eval(rmZeroBiomassQuote),
  #                                                 .(speciesCode, initialEcoregionCode, pixelIndex)])
  availableCombinations <- unique(pixelCohortData[, .(speciesCode, initialEcoregionCode, pixelIndex)])
  pseudoSpeciesEcoregion <- unique(availableCombinations[, .(speciesCode, initialEcoregionCode)])
  newLCCClasses <- Cache(convertUnwantedLCC, classesToReplace = P(sim)$LCCClassesToReplaceNN,
                         rstLCC = rstLCCAdj, availableERC_by_Sp = availableCombinations)

  ## split pixelCohortData into 2 parts -- one with the former 34:36 pixels, one without
  #    The one without 34:36 can be used for statistical estimation, but not the one with
  cohortData34to36 <- pixelCohortData[pixelIndex %in% newLCCClasses$pixelIndex]
  cohortData34to36 <- merge(newLCCClasses, cohortData34to36, all.x = TRUE,
                            all.y = FALSE, by = "pixelIndex")
  cohortDataNo34to36 <- pixelCohortData[!pixelIndex %in% newLCCClasses$pixelIndex]
  setnames(cohortDataNo34to36, "initialEcoregionCode", "ecoregionGroup")
  #cohortDataNo34to36[, ecoregionGroup := initialEcoregionCode]
  cohortDataNo34to36NoBiomass <- cohortDataNo34to36[eval(rmZeroBiomassQuote),
                                                    .(B, logAge, speciesCode, ecoregionGroup, lcc, cover)]

  assert1(cohortData34to36, pixelCohortData[eval(rmZeroBiomassQuote)])

  ##############################################################
  # Statistical estimation of establishprob, maxB and maxANPP
  ##############################################################
  cohortDataShort <- cohortDataNo34to36[, list(coverNum = .N,
                                               coverPres = sum(cover > 0)),
                                        by = c("ecoregionGroup", "speciesCode")]
  cohortDataShortNoCover <- cohortDataShort[coverPres == 0] #
  cohortDataShort <- cohortDataShort[coverPres > 0] # remove places where there is 0 cover
  # will be added back as establishprob = 0
  message(blue("Estimating Species Establishment Probability using P(sim)$coverModel, which is\n",
               magenta(paste0(format(P(sim)$coverModel, appendLF = FALSE), collapse = ""))))

  # for backwards compatibility -- change from parameter to object
  if (is.null(sim$cloudFolderID))
    if (!is.null(P(sim)$cloudFolderID))
      sim$cloudFolderID <- P(sim)$cloudFolderID

  useCloud <- if (!is.null(sim$cloudFolderID)) {
    (getOption("reproducible.useCache", FALSE) && P(sim)$useCloudCacheForStats)
  } else {
    FALSE
  }

  modelCover <- cloudCache(statsModel,
                           modelFn = P(sim)$coverModel,
                           uniqueEcoregionGroup = .sortDotsUnderscoreFirst(unique(cohortDataShort$ecoregionGroup)),
                           .specialData = cohortDataShort,
                           useCloud = useCloud,
                           cloudFolderID = sim$cloudFolderID,
                           showSimilar = getOption("reproducible.showSimilar", FALSE),
                           omitArgs = c("showSimilar", ".specialData",
                                        "useCloud", "cloudFolderID"))
  message(blue("  The rsquared is: "))
  print(modelCover$rsq)

  ## For biomass
  ### Subsample cases where there are more than 50 points in an ecoregionGroup * speciesCode
  cohortDataNo34to36NoBiomass <- Cache(subsetDT,
                                       DT = cohortDataNo34to36NoBiomass,
                                       by = c("ecoregionGroup", "speciesCode"),
                                       doSubset = P(sim)$subsetDataBiomassModel)

  ### For Cache -- doesn't need to cache all columns in the data.table -- only the ones in the model
  ### force parameter values to avoid more checks
  message(blue("Estimating biomass using P(sim)$biomassModel as:\n"),
          magenta(paste0(format(P(sim)$biomassModel, appendLF = FALSE), collapse = "")))
  modelBiomass <- cloudCache(statsModel,
                             modelFn = P(sim)$biomassModel,
                             uniqueEcoregionGroup = .sortDotsUnderscoreFirst(unique(cohortDataNo34to36NoBiomass$ecoregionGroup)),
                             .specialData = cohortDataNo34to36NoBiomass,
                             useCloud = useCloud,
                             cloudFolderID = sim$cloudFolderID,
                             showSimilar = getOption("reproducible.showSimilar", FALSE),
                             omitArgs = c("showSimilar", ".specialData",
                                          "useCloud", "cloudFolderID"))

  message(blue("  The rsquared is: "))
  print(modelBiomass$rsq)

  ########################################################################
  # create speciesEcoregion -- a single line for each combination of ecoregionGroup & speciesCode
  #   doesn't include combinations with B = 0 because those places can't have the species/ecoregion combo
  ########################################################################
  message(blue("Create speciesEcoregion"))
  speciesEcoregion <- makeSpeciesEcoregion(cohortData = cohortDataNo34to36NoBiomass,
                                           species = sim$species)
  # joinOn <- c("ecoregionGroup", "speciesCode")
  # speciesEcoregion <- unique(cohortDataNo34to36NoBiomass, by = joinOn)
  # speciesEcoregion[, c("B", "logAge", "cover") := NULL]
  # speciesEcoregion[lcc %in% unique(cohortDataNo34to36NoBiomass$lcc)] # shouldn't do anything because already correct
  # sim$species[, speciesCode := as.factor(species)]
  # speciesEcoregion <- sim$species[, .(speciesCode, longevity)][speciesEcoregion, on = "speciesCode"]
  # speciesEcoregion[ , ecoregionGroup := factor(as.character(ecoregionGroup))]

  ########################################################################
  # Make predictions from statistical models for
  ########################################################################
  # establishprob -- already is on the short dataset -- need to add back the zeros too
  establishprobBySuccessionTimestep <- 1 - (1 - modelCover$pred)^P(sim)$successionTimestep
  cohortDataShort[, establishprob := establishprobBySuccessionTimestep]

  ############################################
  # Calc. establishProb
  ############################################
  ## for resprouters, establishProb is calculated as the fraction of predicted cover (establishprobBySuccessionTimestep)
  ## that did not result from resprouting. Both reprouters and non-resprouters can be dealt with at the same time
  ## because resproutprob = 0 for non-resprouters
  cohortDataShort <- sim$species[, .(resproutprob, postfireregen, speciesCode)][cohortDataShort, on = "speciesCode"]
  cohortDataShort[, establishprob := pmax(0, pmin(1, (establishprob * (1 - resproutprob))))]

  cohortDataShort <- rbindlist(list(cohortDataShort, cohortDataShortNoCover),
                               use.names = TRUE, fill = TRUE)
  cohortDataShort[is.na(establishprob), establishprob := 0]

  # Join cohortDataShort with establishprob predictions to speciesEcoregion
  speciesEcoregion <- cohortDataShort[, .(ecoregionGroup, speciesCode, establishprob)][
    speciesEcoregion, on = joinOn]

  ########################################################################
  # maxB
  # Set age to the age of longevity and cover to 100%
  speciesEcoregion[, `:=`(logAge = log(longevity), cover = 100)]
  speciesEcoregion[ , maxB := asInteger(predict(modelBiomass$mod,
                                                newdata = speciesEcoregion,
                                                type = "response"))]
  speciesEcoregion[maxB < 0, maxB := 0L] # fix negative predictions

  ########################################################################
  # maxANPP
  message(blue("Add maxANPP to speciesEcoregion -- currently --> maxB/30"))
  speciesEcoregion[ , maxANPP := asInteger(maxB / 30)]

  ########################################################################
  # Clean up unneeded columns
  ########################################################################
  speciesEcoregion[ , `:=`(logAge = NULL, cover = NULL, longevity = NULL, lcc = NULL)]

  speciesEcoregion[ , year := time(sim)]

  #######################################
  if (!is.na(P(sim)$.plotInitialTime)) {
    uniqueSpeciesNames <- as.character(unique(speciesEcoregion$speciesCode))
    names(uniqueSpeciesNames) <- uniqueSpeciesNames
    speciesEcoregionTable2 <- copy(speciesEcoregion)
    speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionGroup)]
    maxB <- stack(lapply(uniqueSpeciesNames, function(sp) {
      rasterizeReduced(speciesEcoregionTable2[speciesCode == sp], ecoregionFiles$ecoregionMap,
                       "maxB", "ecoregionInt")
    }))
    curDev <- dev.cur()
    quickPlot::dev(6, width = 18, height = 10)
    Plot(maxB, legendRange = c(0, max(maxValue(maxB))))
    quickPlot::dev(curDev)
  }

  if (ncell(sim$rasterToMatch) > 3e6) .gc()

  ########################################################################
  # Create initial communities, i.e., pixelGroups
  ########################################################################
  # Rejoin back the pixels that were 34 and 35
  pixelCohortData <- rbindlist(list(cohortData34to36, cohortDataNo34to36),
                               use.names = TRUE, fill = TRUE)
  pixelCohortData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # refactor because the "_34" and "_35" ones are still levels

  pixelCohortData[ , `:=`(logAge = NULL, coverOrig = NULL, totalBiomass = NULL,
                          initialEcoregionCode = NULL, cover = NULL, lcc = NULL)]
  pixelCohortData <- pixelCohortData[B > 0]
  cd <- pixelCohortData[, .SD, .SDcols = c("pixelIndex", sim$columnsForPixelGroups)]
  pixelCohortData[, pixelGroup := Cache(generatePixelGroups, cd, maxPixelGroup = 0,
                                        columns = sim$columnsForPixelGroups)]

  ########################################################################
  ## rebuild ecoregion, ecoregionMap objects -- some initial ecoregions disappeared (e.g., 34, 35, 36)
  ## rebuild biomassMap object -- biomasses have been adjusted
  ecoregionsWeHaveParametersFor <- levels(speciesEcoregion$ecoregionGroup)

  pixelCohortData <- pixelCohortData[ecoregionGroup %in% ecoregionsWeHaveParametersFor] # keep only ones we have params for
  pixelCohortData[ , ecoregionGroup := factor(as.character(ecoregionGroup))]
  pixelCohortData[, totalBiomass := asInteger(sum(B)), by = "pixelIndex"]
  sim$ecoregion <- data.table(active = "yes",
                              ecoregionGroup = factor(as.character(unique(pixelCohortData$ecoregionGroup))))

  # Some ecoregions have NO BIOMASS -- so they are not active
  sim$ecoregion[!ecoregionGroup %in% unique(speciesEcoregion$ecoregionGroup), active := "no"]

  pixelData <- unique(pixelCohortData, by = "pixelIndex")
  pixelData[, ecoregionGroup := factor(as.character(ecoregionGroup))] # resorts them in order

  sim$biomassMap <- raster(sim$rasterToMatch)
  sim$biomassMap[pixelData$pixelIndex] <- pixelData$totalBiomass

  sim$ecoregionMap <- raster(ecoregionFiles$ecoregionMap)
  sim$ecoregionMap[pixelData$pixelIndex] <- as.integer(pixelData$ecoregionGroup)
  levels(sim$ecoregionMap) <- data.frame(ID = seq(levels(pixelData$ecoregionGroup)),
                                         ecoregion = gsub("_.*", "", levels(pixelData$ecoregionGroup)),
                                         ecoregionGroup = levels(pixelData$ecoregionGroup),
                                         stringsAsFactors = TRUE)

  sim$minRelativeB <- data.frame(ecoregionGroup = levels(pixelData$ecoregionGroup),
                                 X1 = 0.2, X2 = 0.4, X3 = 0.5,
                                 X4 = 0.7, X5 = 0.9)

  speciesEcoregion[, ecoregionGroup := factor(as.character(ecoregionGroup))]

  sim$speciesEcoregion <- speciesEcoregion

  sim$speciesLayers <- lapply(seq(numLayers(sim$speciesLayers)), function(x) {
    writeRaster(sim$speciesLayers[[x]],
                file.path(outputPath(sim), paste0(names(sim$speciesLayers)[x], ".tif")),
                datatype = "INT2U", overwrite = TRUE)
  }) %>% raster::stack()

  ##############################################################################
  ##  Collapse pixelCohortData to its cohortData : need pixelGroupMap
  sim$pixelGroupMap <- raster(sim$rasterToMatch)
  sim$pixelGroupMap[pixelData$pixelIndex] <- as.integer(pixelData$pixelGroup)

  sim$cohortData <- unique(pixelCohortData, by = c("pixelGroup", columnsForPixelGroups))
  sim$cohortData[ , `:=`(pixelIndex = NULL)]

  message(blue("Create pixelGroups based on: ", paste(columnsForPixelGroups, collapse = ", "),
               "\n  Resulted in", magenta(length(unique(sim$cohortData$pixelGroup))),
               "unique pixelGroup values"))
  LandR::assertERGs(sim$ecoregionMap, cohortData = sim$cohortData,
                    speciesEcoregion = speciesEcoregion,
                    minRelativeB = sim$minRelativeB)

  LandR::assertCohortData(sim$cohortData, sim$pixelGroupMap)

  LandR::assertUniqueCohortData(sim$cohortData, c("pixelGroup", "ecoregionGroup", "speciesCode"))

  message("Done Boreal_LBMRDataPrep: ", Sys.time())
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
    message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
    sim <- objectSynonyms(sim, list(c("studyAreaLarge", "studyArea")))
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
                            rasterToMatch = if (!needRTM) sim$rasterToMatch else NULL, ## TODO: biomass map needs rasterToMatch but it _is_ the rasterToMatch!!
                            maskWithRTM = TRUE,
                            useSAcrs = if (!needRTM) TRUE else FALSE,
                            method = "bilinear",
                            datatype = "INT2U",
                            filename2 = TRUE, overwrite = TRUE,
                            omitArgs = c("destinationPath", "targetFile", cacheTags, "stable"))
  }
  if (needRTM) {
    # if we need rasterToMatch, that means a) we don't have it, but b) we will have biomassMap
    # sim <- objectSynonyms(sim, list(c("rasterToMatch", "biomassMap")))
    sim$rasterToMatch <- sim$biomassMap
    studyArea <- sim$studyArea # temporary copy because it will be overwritten if it is suppliedElsewhere
    message("  Rasterizing the studyArea polygon map")
    if (!is(studyArea, "SpatialPolygonsDataFrame")) {
      dfData <- if (is.null(rownames(studyArea))) {
        polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
        data.frame("field" = as.character(seq_along(length(studyArea))), row.names = polyID)
      } else {
        polyID <- sapply(slot(studyArea, "polygons"), function(x) slot(x, "ID"))
        data.frame("field" = rownames(studyArea), row.names = polyID)
      }
      studyArea <- SpatialPolygonsDataFrame(studyArea, data = dfData)
    }
    if (!identical(crs(studyArea), crs(sim$rasterToMatch))) {
      studyArea <- spTransform(studyArea, crs(sim$rasterToMatch))
      studyArea <- fixErrors(studyArea)
    }
    #TODO: review whether this is necessary (or will break LandWeb if removed) see Git Issue #22
    # layers provided by David Andison sometimes have LTHRC, sometimes LTHFC ... chose whichever
    LTHxC <- grep("(LTH.+C)", names(studyArea), value = TRUE)
    fieldName <- if (length(LTHxC)) {
      LTHxC
    } else {
      if (length(names(studyArea)) > 1) {
        ## study region may be a simple polygon
        names(studyArea)[1]
      } else NULL
    }

    sim$rasterToMatch <- crop(fasterizeFromSp(studyArea, sim$rasterToMatch, fieldName),
                              studyArea)
    sim$rasterToMatch <- Cache(writeRaster, sim$rasterToMatch,
                               filename = file.path(dataPath(sim), "rasterToMatch.tif"),
                               datatype = "INT2U", overwrite = TRUE)
  }

  if (ncell(sim$rasterToMatch) < 1e4)
    stop("sim$rasterToMatch is too small, it should have more than 10,000 pixels")


  # rstLCC
  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- Cache(prepInputs,
                        targetFile = lcc2005Filename,
                        archive = asPath("LandCoverOfCanada2005_V1_4.zip"),
                        url = extractURL("rstLCC"),
                        destinationPath = dPath,
                        studyArea = sim$studyArea,
                        rasterToMatch = sim$rasterToMatch,
                        method = "bilinear",
                        datatype = "INT2U",
                        filename2 = TRUE, overwrite = TRUE,
                        userTags = c("prepInputsrstLCC_rtm", currentModule(sim)), # use at least 1 unique userTag
                        omitArgs = c("destinationPath", "targetFile"))

    projection(sim$rstLCC) <- projection(sim$rasterToMatch)
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
                             #filename2 = TRUE,
                             userTags = c("prepInputsEcoDistrict_SA", currentModule(sim), cacheTags), # use at least 1 unique userTag
                             omitArgs = c("destinationPath", "targetFile", "overwrite", "alsoExtract"))
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
                             filename2 = NULL,
                             overwrite = TRUE,
                             userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags), # use at least 1 unique userTag
                             omitArgs = c("destinationPath", "targetFile", "overwrite", "alsoExtract"))
    sim$standAgeMap[] <- asInteger(sim$standAgeMap[])
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      stop("If you provide sppColorVect, you MUST also provide sppEquiv")

    data("sppEquivalencies_CA", package = "LandR", envir = environment())
    sim$sppEquiv <- as.data.table(sppEquivalencies_CA)
    ## By default, Abies_las is renamed to Abies_sp
    sim$sppEquiv[KNN == "Abie_Las", LandR := "Abie_sp"]

    ## check spp column to use
    if (P(sim)$sppEquivCol == "Boreal") {
      message(paste("There is no 'sppEquiv' table supplied;",
                    "will attempt to use species listed under 'Boreal'",
                    "in the 'LandR::sppEquivalencies_CA' table"))
    } else {
      if (grepl(P(sim)$sppEquivCol, names(sim$sppEquiv))) {
        message(paste("There is no 'sppEquiv' table supplied,",
                      "will attempt to use species listed under", P(sim)$sppEquivCol,
                      "in the 'LandR::sppEquivalencies_CA' table"))
      } else {
        stop("You changed 'sppEquivCol' without providing 'sppEquiv',",
             "and the column name can't be found in the default table ('LandR::sppEquivalencies_CA').",
             "Please provide conforming 'sppEquivCol', 'sppEquiv' and 'sppColorVect'")
      }
    }

    ## remove empty lines/NAs
    sim$sppEquiv <- sim$sppEquiv[!"", on = P(sim)$sppEquivCol]
    sim$sppEquiv <- na.omit(sim$sppEquiv, P(sim)$sppEquivCol)

    ## add default colors for species used in model
    sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                  newVals = "Mixed", palette = "Accent")
  } else {
    if (is.null(sim$sppColorVect))
      stop("If you provide 'sppEquiv' you MUST also provide 'sppColorVect'")
  }

  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    sim$speciesLayers <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppEquivCol = P(sim)$sppEquivCol,
                               thresh = 5,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"))
  }

  # 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath = dPath, cacheTags = cacheTags)
    ## override longevity values - from
    sim$speciesTable[LandisCode == "PICE.GLA", Longevity := 400]
    sim$speciesTable[LandisCode == "PINU.CON.LAT", Longevity := 335]
    sim$speciesTable[LandisCode == "PICE.MAR", Longevity := 250]
    sim$speciesTable[LandisCode == "POPU.TRE", Longevity := 200]
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

  if (!suppliedElsewhere("columnsForPixelGroups", sim)) {
    sim$columnsForPixelGroups <- LandR::columnsForPixelGroups
  }

  return(invisible(sim))
}
