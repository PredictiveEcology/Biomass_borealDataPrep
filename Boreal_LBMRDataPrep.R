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
                    quote(biomass ~ logAge * speciesCode + (speciesCode | ecoregionGroup) + cover),
                    NA, NA,
                    "This formula is for estimating biomass from landcover, ecoregionGroup, speciesCode, age, and cover"),
    defineParameter("coverQuotedFormula", "name",
                    quote(cbind(coverPres, coverNum) ~ speciesCode + (1 | ecoregionGroup)),
                    NA, NA,
                    "This formula is for estimating biomass from landcover, ecoregion, and speciesCode"),
    defineParameter("pixelGroupAgeClass", "numeric", P(sim)$successionTimestep, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same"),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same"),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter("useCloudCacheForStats", "logical", TRUE, NA, NA,
                    "Some of the statistical models take long (at least 30 minutes, likely longer). If this is TRUE, then it will try to get previous cached runs from googledrive"),
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
  rstEcoregionMap <- fasterize::fasterize(sf::st_as_sf(ecoregionMap), raster = sim$rasterToMatch,
                                          field = "ECODISTRIC")
  ecoregionstatus <- data.table(active = "yes",
                                ecoregion = 1:1031)
  #  also rm 37, 38, 39 --> make them NA
  sim$LCC2005Adj <- sim$LCC2005

  # Rm rock and ice pixels
  pixelsToRm <- sim$LCC2005[] %in% 37:39 # these are lakes, rock and ice
  sim$rasterToMatch[pixelsToRm] <- NA
  sim$LCC2005Adj[pixelsToRm] <- NA
  rstEcoregionMap[pixelsToRm] <- NA

  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(rstEcoregionMap, sim$LCC2005Adj),
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
                   initialEcoregionCode = factor(factorValues2(ecoregionFiles$ecoregionMap,
                                                        ecoregionFiles$ecoregionMap[],
                                                        att = 5)),
                   totalBiomass = asInteger(sim$biomassMap[]) * 100, # change units
                   cover = sim$speciesLayers[],
                   pixelIndex = seq(ncell(sim$standAgeMap)),
                   lcc = factor(sim$LCC2005Adj[])
  )

  message(crayon::blue("This is the summary of the input data for age, ecoregionGroup, biomass, speciesLayers:"))
  print(summary(dt))
  dt <- dt[!is.na(lcc)]
  message(crayon::blue("Step 2: rm NAs, leaving", crayon::magenta(NROW(dt)), "pixels with data"))

  ### Create groupings
  coverColNames <- grep(colnames(dt), pattern = "cover", value = TRUE)
  newCoverColNames <- gsub("cover\\.", "", coverColNames)
  setnames(dt, old = coverColNames, new = newCoverColNames)
  message(crayon::blue("Step 3: Create initial cohortData object, with no pixelGroups yet"))

  cohortData <- data.table::melt(dt,
                                 value.name = "cover",
                                 measure.vars = newCoverColNames,
                                 variable.name = "speciesCode")
  cohortData[, coverOrig := cover]

  if (getOption("LandR.assertions"))
    #describeCohortData(cohortData)

  message(blue("Step 4a: assign biomass = 0 and age = 0 for pixels where cover = 0, ",
               "\n  because cover is most reliable dataset"))
  cohortData[cover == 0, `:=`(age = 0L, biomass = 0L)]
  message(blue("Step 4b: assign totalBiomass = 0 sum(cover) = 0 in a pixel, ",
               "\n  because cover is most reliable dataset"))
  cohortData <- cohortData[, sum(cover)==0, by = "pixelIndex"][V1 == TRUE][cohortData, on = "pixelIndex"][V1 == TRUE, totalBiomass := 0L]
  cohortData[, V1 := NULL]

  #####################
  ######################
  message(crayon::blue("Step 5: POSSIBLE ALERT -- assume deciduous cover is 1/2 the conversion to biomass as conifer"))
  cohortData[speciesCode == "Popu_sp", cover := asInteger(cover / 2)] # CRAZY TODO -- DIVIDE THE COVER BY 2 for DECIDUOUS -- will only affect mixed stands
  cohortData[ , cover := {
    sumCover <- sum(cover)
    if (sumCover > 100) {
      cover <- asInteger(cover/(sumCover + 0.0001) * 100L)
    }
    cover
  }, by = "pixelIndex"]


  # Biomass -- by cohort
  message(crayon::blue("Step 6: Divide total biomass of each pixel by the relative cover of the cohorts"))
  cohortData[ , biomass := asInteger(mean(totalBiomass) * cover / 100), by = "pixelIndex"] # /100 because cover is percent
  message(crayon::blue("Step 7: Round biomass to nearest P(sim)$pixelGroupBiomassClass"))
  cohortData[ , biomass := asInteger(ceiling(biomass / P(sim)$pixelGroupBiomassClass) *
                                       P(sim)$pixelGroupBiomassClass)] # /100 because cover is percent
  message(blue("Step 8a - Set biomass to 0 where cover > 0 and age = 0, because biomass is least quality dataset"))
  cohortData[ , totalBiomass := asInteger(totalBiomass)]

  ######################################################
  # Impute missing age
  ######################################################
  cohortDataMissingAge <- cohortData[, hasBadAge := all(age == 0 & cover > 0) | any(is.na(age)), by = "pixelIndex"][
    hasBadAge == TRUE]
  cohortDataMissingAgeUnique <- unique(cohortDataMissingAge,
                                         by = c("initialEcoregionCode", "speciesCode"))[
                                         , .(initialEcoregionCode, speciesCode)]
  cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[cohortData, on = c("initialEcoregionCode", "speciesCode"), nomatch = 0]
  ageQuotedFormula <- quote(age ~ biomass * speciesCode + (1 | initialEcoregionCode) + cover)
  cohortDataMissingAgeUnique <- cohortDataMissingAgeUnique[, .(biomass, age, speciesCode, initialEcoregionCode, cover)]
  system.time(outAge <- Cache(statsModel, form = ageQuotedFormula, .specialData = cohortDataMissingAgeUnique))

  print(outAge$rsq)

  cohortDataMissingAge[, imputedAge := pmax(0L, asInteger(predict(outAge$mod, newdata = cohortDataMissingAge)))]
  cohortData <- cohortDataMissingAge[, .(pixelIndex, imputedAge, speciesCode)][cohortData, on = c("pixelIndex", "speciesCode")]
  cohortData[!is.na(imputedAge), `:=`(age = imputedAge, logAge = log(imputedAge))]
  cohortData[, `:=`(imputedAge = NULL, hasBadAge = NULL)]

  #######################################################
  #
  #######################################################
  message(blue("Step 8b - Set recalculate totalBiomass as sum(biomass); many biomasses will have been set to 0 in previous steps"))
  cohortData[cover > 0 & age == 0, biomass := 0L]
  cohortData[, totalBiomass := sum(biomass), by = "pixelIndex"]

  # https://stats.stackexchange.com/questions/31300/dealing-with-0-1-values-in-a-beta-regression
  # cohortData[ , coverProp := (cover/100 * (NROW(cohortData) - 1) + 0.5) / NROW(cohortData)]

  ######### CHANGE 34 and 35 and 36 values -- burns and cities
  #     replace these with a neighbour class *that exists*
  browser()
  out3 <- convertUnwantedLCC(pixelClassesToReplace = 34:36,
                             rstLCC = sim$LCC2005Adj,
                             pixelCohortData = cohortData)

  cohortData34to36 <- cohortData[pixelIndex %in% out3$pixelIndex]
  cohortData34to36 <- out3[cohortData34to36, on = "pixelIndex"]

  if (getOption("LandR.assertions")) {
    allCodesAre34to36 <- all(grepl(".*34|.*35|.*36",
                                   as.character(cohortData34to36$initialEcoregionCode)))
    if (!allCodesAre34to36)
      stop("lcc classes were mismanaged; contact developers: code 234")
  }

  if (getOption("LandR.assertions")) {
    onlyExistingCodes <- all(unique(cohortData34to36$ecoregionGroup) %in% unique(cohortData$initialEcoregionCode))
    if (!onlyExistingCodes)
      stop("There are some ecoregionCodes created post replacement of 34 and 35")
  }

  cohortDataNo34to36 <- cohortData[!pixelIndex %in% out3$pixelIndex]
  cohortDataNo34to36[, ecoregionGroup := initialEcoregionCode]

  ##############################################################
  # Statistical estimation of SEP, maxBiomass and maxANPP
  ##############################################################
  cohortDataShort <- cohortDataNo34to36[, list(coverNum = .N,
                               coverPres = sum(cover > 0)),
                        by = c("ecoregionGroup", "speciesCode", "lcc")]

  message(crayon::blue("Step 9: Estimaing Species Establishment Probability using P(sim)$coverQuotedFormula, which is\n",
                       format(P(sim)$coverQuotedFormula)))

  modelCover <- cloudCache(statsModel, P(sim)$coverQuotedFormula, cohortDataShort, family = binomial,
                         # checksumsFileID = "1XznvuxsRixGxhYCicMr5mdoZlyYECY8C",
                         useCloud = P(sim)$useCloudCacheForStats,
                         cloudFolderID = "/folders/1wJXDyp5_XL2RubViWGAeTNDqGElfgkL8")
  message(crayon::blue("  The rsquared is: "))
  print(modelCover$rsq)


  # For biomass
  # For Cache -- doesn't need to cache all columns in the data.table -- only the ones in the model
  cohortDataNo34to36NoBiomass <- cohortDataNo34to36[biomass > 0, .(biomass, logAge, speciesCode, ecoregionGroup, lcc, cover)]
  message(crayon::blue("Step 10: Estimaing maxBiomass with P(sim)$biomassQuotedFormula, which is:\n",
          magenta(paste0(format(P(sim)$biomassQuotedFormula, appendLF = FALSE), collapse = ""))))
  modelBiomass <- cloudCache(statsModel, form = P(sim)$biomassQuotedFormula, .specialData = cohortDataNo34to36NoBiomass,
                         useCloud = P(sim)$useCloudCacheForStats,
                         cloudFolderID = "/folders/1wJXDyp5_XL2RubViWGAeTNDqGElfgkL8")
  message(crayon::blue("  The rsquared is: "))
  print(modelBiomass$rsq)

  # Create initial communities, i.e., pixelGroups
  # Rejoin back the pixels that were 34 and 35
  cohortData <- rbindlist(list(cohortData34to36, cohortDataNo34to36), use.names = TRUE, fill = TRUE)
  cohortData[, ecoregionGroup := factor(ecoregionGroup)] # refactor because the "_34" and "_35" ones are still levels
  columnsForPG <- c("ecoregionGroup", "speciesCode", "age", "biomass")

  cd <- cohortData[,c("pixelIndex", columnsForPG), with = FALSE]
  cohortData[, pixelGroup :=
               Cache(generatePixelGroups, cd, maxPixelGroup = 0,
                     columns = columnsForPG)]
  message(crayon::blue("Step 11: Create pixelGroups based on: ",
                       paste(columnsForPG, collapse = ", "),
                       "\n  Resulted in", magenta(length(unique(cohortData$pixelGroup))),
                       "unique pixelGroup values"))

  ########################################################################
  # create speciesEcoregionTable -- a single line for each combination of ecoregionGroup & speciesCode
  #   doesn't include combinations with biomass = 0 because those places can't have the species/ecoregion combo
  ########################################################################
  joinOn <- c("ecoregionGroup", "speciesCode")
  speciesEcoregionTable <- unique(cohortDataNo34to36NoBiomass, by = joinOn)
  speciesEcoregionTable[, c("biomass", "logAge", "cover") := NULL]
  speciesEcoregionTable[lcc %in% unique(cohortDataNo34to36NoBiomass$lcc)] # shouldn't do anything because already correct
  sim$species[, speciesCode := as.factor(species)]
  speciesEcoregionTable <- sim$species[, .(speciesCode, longevity)][speciesEcoregionTable, on = "speciesCode"]


  ########################################################################
  # Make predictions from statistical models for
  ########################################################################
  # SEP -- already is on the short dataset
  ########################################################################
  cohortDataShort[, SEP := modelCover$pred]

  ########################################################################
  # maxBiomass
  ########################################################################
  # Set age to the age of longevity and cover to 100%
  speciesEcoregionTable[, `:=`(logAge = log(longevity), cover = 100)]
  speciesEcoregionTable[ , maxBiomass := asInteger(predict(modelBiomass$mod,
                                                           newdata = speciesEcoregionTable,
                                                           type = "response"))]
  speciesEcoregionTable[maxBiomass < 0, maxBiomass := 0] # fix negative predictions
  message(crayon::blue("Step 12: Add maxANPP to speciesEcoregionTable -- currently --> maxBiomass/30"))

  ########################################################################
  # maxANPP
  ########################################################################
  speciesEcoregionTable[ , maxANPP := maxBiomass / 30]

  # Join cohortDataShort with SEP predictions to speciesEcoregionTable
  speciesEcoregionTable <- cohortDataShort[, .(ecoregionGroup, speciesCode, SEP)][speciesEcoregionTable, on = joinOn]

  ########################################################################
  # Clean up unneeded columns
  ########################################################################
  speciesEcoregionTable[ , `:=`(logAge = NULL, cover = NULL, longevity = NULL, #pixelIndex = NULL,
                                lcc = NULL)]

  #######################################
  if (!is.na(P(sim)$.plotInitialTime)) {
    uniqueSpeciesNames <- as.character(unique(speciesEcoregionTable$speciesCode))
    names(uniqueSpeciesNames) <- uniqueSpeciesNames
    speciesEcoregionTable2 <- copy(speciesEcoregionTable)
    speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionGroup)]
    maxBiomass <- stack(lapply(uniqueSpeciesNames, function(sp) {
      rasterizeReduced(speciesEcoregionTable2[speciesCode == sp], ecoregionFiles$ecoregionMap,
                       "maxBiomass", "ecoregionInt")
    }))
    Plot(maxBiomass, legendRange = c(0, max(maxValue(maxBiomass))))
  }

  if (ncell(sim$rasterToMatch) > 3e6) .gc()

  ## Assign cohortData as initialCommunities
  sim$initialCommunities <-  cohortData

  ## rebuild ecoregion and ecoregionMap objects -- some initial ecoregions disappeared (e.g., 34, 35, 36)
  sim$ecoregion <- data.table(active = "yes", ecoregionGroup = factor(levels(cohortData$ecoregionGroup)))

  pixelData <- unique(cohortData, by = "pixelIndex")
  sim$ecoregionMap <- raster(ecoregionFiles$ecoregionMap)
  sim$ecoregionMap[pixelData$pixelIndex] <- as.integer(pixelData$ecoregionGroup)
  levels(sim$ecoregionMap) <- data.frame(ID = seq(levels(pixelData$ecoregionGroup)),
                                         ecoregionGroup = levels(pixelData$ecoregionGroup),
                                         stringsAsFactors = TRUE)

  sim$minRelativeB <- data.frame(ecoregionGroup = unique(cohortData$ecoregionGroup),
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



statsModel <- function(form, .specialData, ...) {
  if ("family" %in% names(list(...))) {
    modelFn <- glmer
  } else {
    modelFn <- lmer
  }
  mod <- modelFn(
    formula = eval(form),
    data = .specialData,
    ...)

  list(mod = mod, pred = fitted(mod),
       rsq = MuMIn::r.squaredGLMM(mod))
}
