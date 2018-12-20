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
  version = numeric_version("1.3.3"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Boreal_LBMRDataPrep.Rmd"),
  reqdPkgs = list("data.table", "dplyr", "fasterize", "gdalUtils", "raster", "rgeos", "sp",
                  #"PredictiveEcology/LandR@development",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    #defineParameter("speciesPresence", "numeric", 50, NA, NA,
    #                "minimum percent cover required to classify a species as present"),
    defineParameter("minNumPixelsToEstMaxBiomass", "integer", 100, NA, NA,
                    "When estimating maximum biomass by species and ecoregion, this number indicates the minimum number of pixels with data required before a maximum is estimated."),
    defineParameter("quantileForMaxBiomass", "numeric", 0.99, NA, NA,
                    "When estimating maximum biomass by species and ecoregion, rather than take the absolute max(biomass), the quantile is taken. This gives the capacity to remove outliers."),
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
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
                 desc = "table of species equivalencies. See pemisc::sppEquivalencies_CA.", ## TODO: use LandR
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
    sim <- scheduleEvent(sim, P(sim)$.plotInitialTime, "Boreal_LBMRDataPrep", "plot")
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

  rasterToMatchBinary <- raster(sim$rasterToMatch)
  rasterToMatchBinary[] <- NA
  rasterToMatchBinary[!is.na(sim$rasterToMatch[])] <- 1

  message("2: ecoregionProducer: ", Sys.time())
  # Note: this ecoregionMap is NOT the Canadian EcoRegion -- it is for LBMR, which uses "ecoregion"
  ecoregionMap <- Cache(postProcess, sim$ecoDistrict, studyArea = sim$studyArea, filename2 = NULL)
  ecoregionstatus <- data.table(active = "yes",
                                ecoregion = 1:1031)
  ecoregionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(ecoregionMap, sim$rasterToMatch),
                          ecoregionName = "ECODISTRIC",
                          ecoregionActiveStatus = ecoregionstatus,
                          rasterToMatch = sim$rasterToMatch,
                          userTags = "stable")

  ecoRegionMap <- Cache(postProcess, sim$ecoRegion, studyArea = sim$studyArea, filename2 = NULL)
  ecoRegionFiles <- Cache(ecoregionProducer,
                          ecoregionMaps = list(ecoRegionMap, sim$rasterToMatch),
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
                            # This cuts up species pct into x groups, and age into x-year groups
                            percentileGrps = 10,
                            ecoregionMap = ecoregionFiles$ecoregionMap,
                            standAgeMap = sim$standAgeMap,
                            userTags = "stable")
  # remove unneeded columns of initialCommunities based on sim$species
  initialCommunities <- initialCommunities[, .(mapcode,
                                               description = NA,
                                               species,
                                               speciesPresence = as.integer(speciesPresence),
                                               age = as.integer(age1),
                                               pixelIndex = as.integer(pixelIndex))]

  ############################################################
  # Create initialCommunitiesMap
  ############################################################
  initialCommunitiesMap <- raster(speciesLayers[[1]])
  initialCommunitiesMap[initialCommunities$pixelIndex] <- as.integer(initialCommunities$mapcode) # integer is OK now that it is factor
  datatype <- if (length(levels(initialCommunities$mapCodeFac)) > 64e3)
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


  if (ncell(sim$rasterToMatch) > 3e6) .gc()

  ################################################################
  ## species traits inputs
  ################################################################
  message("6: Prepare 'species' table, i.e., species level traits", Sys.time())
  sim$species <- prepSpeciesTable(speciesTable = sim$speciesTable,
                                  speciesLayers = sim$speciesLayers,
                                  sppEquiv = sim$sppEquiv,
                                  sppEquivCol = P(sim)$sppEquivCol)


  ## filter communities to species that have traits
  sim$initialCommunities <- initialCommunities[initialCommunities$species %in% sim$species$species,]



  sim$minRelativeB <- data.frame(ecoregion = sim$ecoregion[active == "yes",]$ecoregion,
                                 X1 = 0.2, X2 = 0.4, X3 = 0.5,
                                 X4 = 0.7, X5 = 0.9)

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
                            rasterToMatch = sim$rasterToMatch,
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
    sim$standAgeMap <- Cache(prepInputs, #notOlderThan = Sys.time(),
                             targetFile = basename(standAgeMapFilename),
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
  }

  if (!suppliedElsewhere("sppEquiv", sim)) {
    data("sppEquivalencies_CA", package = "pemisc", envir = environment()) ## TODO: use LandR
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
                               sppEquivCol = "LandR",
                               # thresh = 10,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"))

    #options(opts)
    writeRaster(speciesLayersList$speciesLayers,
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
