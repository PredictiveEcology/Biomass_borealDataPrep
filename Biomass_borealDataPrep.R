defineModule(sim, list(
  name = "Biomass_borealDataPrep",
  description = "A data preparation module for parameterizing Biomass_core from open data sources, within the Boreal forest of Canada",
  keywords = c("LandWeb", "Biomass_core"),
  authors = c(
    person("Yong", "Luo", email = "yong.luo@canada.ca", role = c("aut")),
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@canada.ca", role = c("aut", "cre")),
    person(c("Ceres"), "Barros", email = "cbarros@mail.ubc.ca", role = c("ctb")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("ctb"))
  ),
  childModules = character(0),
  version = list(Biomass_borealDataPrep = numeric_version("1.4.0.9000"),
                 LandR = "0.0.3.9004", SpaDES.core = "1.0.0",
                 reproducible = "1.0.0.9006"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_borealDataPrep.Rmd"),
  reqdPkgs = list("crayon", "data.table", "dplyr", "fasterize", "plyr", "raster", "sp", "sf",
                  "SpaDES.tools", "reproducible",
                  "achubaty/amc@development",
                  "PredictiveEcology/LandR@development",
                  "PredictiveEcology/pemisc@development", "merTools"),
  parameters = rbind(
    defineParameter("biomassModel", "call",
                    quote(lme4::lmer(B ~ logAge * speciesCode + cover * speciesCode +
                                       (logAge + cover + speciesCode | ecoregionGroup))),
                    NA, NA,
                    paste("Model and formula for estimating biomass (B) from ecoregionGroup",
                          "(currently ecoDistrict * LandCoverClass), speciesCode,",
                          "logAge (gives a downward curving relationship), and cover.",
                          "Defaults to a LMEM, which can be slow if dealing with very large datasets",
                          "(e.g. 36 000 points take 20min).",
                          "For faster fitting try P(sim)$subsetDataBiomassModel == TRUE, or",
                          "quote(RcppArmadillo::fastLm(formula = B ~ logAge * speciesCode * ecoregionGroup",
                          "+ cover * speciesCode * ecoregionGroup)).",
                          "A custom model call can also be provided, as long as the 'data' argument",
                          "is NOT included.")),
    defineParameter("coverModel", "call",
                    quote(glm(cbind(coverPres, coverNum - coverPres) ~ speciesCode * ecoregionGroup,
                              family = binomial)),
                    NA, NA,
                    paste("Model and formula used for estimating cover from ecoregion and speciesCode",
                          "and potentially others. Defaults to a GLMEM if there are > 1 grouping levels.",
                          "A custom model call can also be provided, as long as the 'data' argument is NOT included")),

    # deciduous cover to biomass cover section ################################################
    defineParameter("coverPctToBiomassPctModel", "call",
                    quote(glm(I(log(B/100)) ~ logAge * I(log(totalBiomass/100)) * speciesCode * lcc)),
                    NA, NA,
                    paste("Model to estimate the relationship between % cover and % biomass, referred to as",
                          "deciduousCoverDiscount. It is a number between 0 and 1 that translates %cover,",
                          "as provided in several databases, to %biomass. It is assumed that all hardwoods",
                          "are equivalent and all softwoods are equivalent and that %cover of hardwoods will",
                          "be an overesimate of the %biomass of hardwoods. E.g., 30%cover of hardwoods",
                          "might translate to 20% biomass of hardwoods. The reason this discount exists is",
                          "because hardwoods in Canada have a much wider canopy than softwoods.")),
    defineParameter("fitDecidiousCoverDiscount", "logical",FALSE, NA, NA,
                    paste("If TRUE, this will re-estimate deciduousCoverDiscount. This may be unstable and",
                          "is not recommended currently. If FALSE, will use the current default")),
    defineParameter("deciduousCoverDiscount", "numeric",0.8418911, NA, NA,
                    paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                          "Will not be used if P(sim)$fitDeciduousCoverDiscount is TRUE")),
    ###########################################################################################

    defineParameter("ecoregionLayerField", "character", NULL, NA, NA,
                    paste("the name of the field used to distinguish ecoregions, if supplying a polygon.",
                          "The default is 'ECODISTRIC' where available (for legacy reasons), else the row numbers of",
                          "sim$ecoregionLayer. If this field is not numeric, it will be coerced to numeric")),
    defineParameter("forestedLCCClasses", "numeric", c(1:15, 20, 32, 34:35), 0, 39,
                    paste("The classes in the rstLCC layer that are 'treed' and will therefore be run in Biomass_core.",
                          "Defaults to forested classes in LCC2005 map.")),
    defineParameter("imputeBadAgeModel", "call",
                    quote(lme4::lmer(age ~ log(totalBiomass) * cover * speciesCode + (log(totalBiomass) * speciesCode | initialEcoregionCode))),
                    NA, NA,
                    paste("Model and formula used for imputing ages that are either missing or do not match well with",
                          "Biomass or Cover. Specifically, if Biomass or Cover is 0, but age is not, then age will be imputed.",
                          "Similarly, if Age is 0 and either Biomass or Cover is not, then age will be imputed")),
    defineParameter("LCCClassesToReplaceNN", "numeric", 34:35, NA, NA,
                    paste("This will replace these classes on the landscape with the closest forest class P(sim)$forestedLCCClasses.",
                          "If the user is using the default 2005 data product for rstLCC, then users may wish to",
                          "include 36 (cities -- if running a historic range of variation project), and 34:35 (burns)",
                          "Since this is about estimating parameters for growth, it doesn't make any sense to have",
                          "unique estimates for transient classes in most cases")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    "Cover that is equal to or below this number will be omitted from the dataset"),
    defineParameter("omitNonTreedPixels", "logical", TRUE, FALSE, TRUE,
                    "Should this module use only treed pixels, as identified by P(sim)$forestedLCCClasses?"),
    defineParameter("pixelGroupAgeClass", "numeric", params(sim)$Biomass_borealDataPrep$successionTimestep, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of ages that will be considered 'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same"),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    "When assigning pixelGroup membership, this defines the resolution of biomass that will be considered 'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same"),
    defineParameter("speciesUpdateFunction", "list",
                    list(quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable, sim$sppEquiv, P(sim)$sppEquivCol))),
                    NA, NA,
                    paste("Unnamed list of quoted functions that updates species table to customize values.",
                          "Default should always come first.")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in sim$specieEquivalency data.table to use as a naming convention"),
    defineParameter("subsetDataAgeModel", "numeric", 50, NA, NA,
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
                    "This describes the simulation time interval between save events"),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bind_rows(
    expectsInput("cloudFolderID", "character",
                 "The google drive location where cloudCache will store large statistical objects"),
    expectsInput("columnsForPixelGroups", "character",
                 "The names of the columns in cohortData that define unique pixelGroups. Default is c('ecoregionGroup', 'speciesCode', 'age', 'B') "),
    expectsInput("ecoregionLayer", "SpatialPolygonsDataFrame",
                 desc = paste("A SpatialPolygonsDataFrame that characterizes the unique ecological regions used to",
                              "parameterize the biomass, cover, and species establishment probability models.",
                              "It will be overlaid with landcover to generate classes for every ecoregion/LCC combination.",
                              "It must have same extent and crs as studyAreaLarge if suppplied by user.",
                              "It is superseded by sim$ecoregionRst if that object is supplied by the user"),
                 sourceURL = "http://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput('ecoregionRst', "RasterLayer",
                 desc = paste("A raster that characterizes the unique ecological regions used to",
                              "parameterize the biomass, cover, and species establishment probability models.",
                              "If this object is provided, it will supersede sim$ecoregionLayer.",
                              "It will be overlaid with landcover to generate classes for every ecoregion/LCC combination.",
                              "It must have same extent and crs as rasterToMatchLarge if suppplied by user - use reproducible::postProcess.",
                              "If it uses an attribute table, it must contain the field 'ecoregion' to represent raster values")),
    expectsInput("fireURL", "character",
                 desc = paste("A url to a fire database, such as the Canadian National Fire Database,",
                              "that is a zipped shapefile with fire polygons, an attribute (i.e., a column) named 'Year'.",
                              "If supplied (omitted with NULL or NA), this will be used to 'update' age pixels on standAgeMap",
                              "with 'time since fire' as derived from this fire polygons map"),
                 sourceURL = "https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2018_20191129.zip"),
    expectsInput("rstLCC", "RasterLayer",
                 desc = paste("A land classification map in study area. It must be 'corrected', in the sense that:\n",
                              "1) Every class must not conflict with any other map in this module\n",
                              "    (e.g., speciesLayers should not have data in LCC classes that are non-treed);\n",
                              "2) It can have treed and non-treed classes. The non-treed will be removed within this\n",
                              "    module if P(sim)$omitNonTreedPixels is TRUE;\n",
                              "3) It can have transient pixels, such as 'young fire'. These will be converted to a\n",
                              "    the nearest non-transient class, probabilistically if there is more than 1 nearest\n",
                              "    neighbour class, based on P(sim)$LCCClassesToReplaceNN.\n",
                              "The default layer used, if not supplied, is Canada national land classification in 2005.",
                              " The metadata (res, proj, ext, origin) need to match rasterToMatchLarge."),
                 sourceURL = "ftp://ftp.ccrs.nrcan.gc.ca/ad/NLCCLandCover/LandcoverCanada2005_250m/LandCoverOfCanada2005_V1_4.zip"),
    expectsInput("rasterToMatch", "RasterLayer",
                 desc = paste("A raster of the studyArea in the same resolution and projection as rawBiomassMap.",
                              "This is the scale used for all *outputs* for use in the simulation."),
                 sourceURL = NA),
    expectsInput("rasterToMatchLarge", "RasterLayer",
                 desc = paste("A raster of the studyAreaLarge in the same resolution and projection as rawBiomassMap.",
                              "This is the scale used for all *inputs* for use in the simulation."),
                 sourceURL = NA),
    expectsInput("rawBiomassMap", "RasterLayer",
                 desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001. If necessary, biomass values are rescaled to match changes in resolution.",
                              "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata."),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("speciesLayers", "RasterStack",
                 desc = paste("cover percentage raster layers by species in Canada species map.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2001 using a cover threshold of 10 -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/")),
    expectsInput("speciesTable", "data.table",
                 desc = "species attributes table, default is from Dominic Cyr and Yan Boulanger's project",
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species",
                 sourceURL = ""),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("sppNameVector", "character",
                 desc = "an optional vector of species names to be pulled from sppEquiv. If not provided, then species will be taken from the entire P(sim)$sppEquivCol in sppEquiv. See LandR::sppEquivalencies_CA.",
                 sourceURL = ""),
    expectsInput("standAgeMap", "RasterLayer",
                 desc =  paste("stand age map in study area.",
                               "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                               "kNN-derived biomass map from 2001 -",
                               "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("studyArea", "SpatialPolygonsDataFrame",
                 desc = paste("Polygon to use as the study area.",
                              "Defaults to  an area in Southwestern Alberta, Canada."),
                 sourceURL = ""),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (larger area than studyArea) used for parameter estimation,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."),
                 sourceURL = "")
  ),
  outputObjects = bind_rows(
    createsOutput("biomassMap", "RasterLayer",
                  desc = paste("total biomass raster layer in study area,",
                               "filtered for pixels covered by cohortData")),
    createsOutput("cohortData", "data.table",
                  desc = paste("initial community table, created from available biomass,",
                               "age and species cover data, as well as eco zonation information")),
    createsOutput("ecoregion", "data.table",
                  desc = "ecoregion look up table"),
    createsOutput("ecoregionMap", "RasterLayer",
                  desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table"),
    createsOutput("pixelGroupMap", "RasterLayer",
                  desc = "initial community map that has mapcodes match initial community table"),
    createsOutput("pixelFateDT", "data.table",
                  desc = paste("A small table that keeps track of the pixel removals and cause. This may help diagnose issues",
                  " related to understanding the creation of cohortData")),
    createsOutput("minRelativeB", "data.frame",
                  desc = "define the cut points to classify stand shadeness"),
    createsOutput("rawBiomassMap", "RasterLayer",
                  desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                               "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                               "from 2001. See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                               "for metadata")),
    createsOutput("species", "data.table",
                  desc = "a table that has species traits such as longevity..."),
    createsOutput("speciesEcoregion", "data.table",
                  desc = "define the maxANPP, maxB and establishprob change with both ecoregion and simulation time"),
    createsOutput("studyArea", "",
                  desc = paste("Polygon to use as the study area.",
                               "Defaults to  an area in Southwestern Alberta, Canada.")),
    createsOutput("sufficientLight", "data.frame",
                  desc = "define how the species with different shade tolerance respond to stand shadeness")
    # createsOutput("speciesEstablishmentProbMap", "RasterStack",
    #               paste("Species establishment probability as a map, ",
    #                     "by species. This is written to disk to save RAM space")),
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Biomass_borealDataPrep <- function(sim, eventTime, eventType, debug = FALSE) {
  if (eventType == "init") {
    sim <- createBiomass_coreInputs(sim)

    # schedule future event(s)
    sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_borealDataPrep", "save")
  } else if (eventType == "save") {
    sim <- Save(sim)
  } else {
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  }
  return(invisible(sim))
}

createBiomass_coreInputs <- function(sim) {
  # # ! ----- EDIT BELOW ----- ! #
  if (is.null(P(sim)$pixelGroupAgeClass))
    params(sim)[[currentModule(sim)]]$pixelGroupAgeClass <- P(sim)$successionTimestep

  cacheTags <- c(currentModule(sim), "init")

  message(blue("Starting to createBiomass_coreInputs in Biomass_borealDataPrep: ", Sys.time()))
  if (is.null(sim$speciesLayers))
    stop(red(paste("'speciesLayers' are missing in Biomass_borealDataPrep init event.\n",
                   "This is likely due to the module producing 'speciesLayers' being scheduled after Biomass_borealDataPrep.\n",
                   "Please check module order.")))

  ## check that input rasters all match
  compareRaster(sim$rasterToMatchLarge, sim$rawBiomassMap, sim$rstLCC,
                sim$speciesLayers, sim$standAgeMap, orig = TRUE)

  # sim$standAgeMap <- round(sim$standAgeMap / 20, 0) * 20 # use 20-year bins (#103)
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
  defaultQuote <- quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable,
                                                  sim$sppEquiv, P(sim)$sppEquivCol))
  if (P(sim)$speciesUpdateFunction[[1]] != defaultQuote) {
    stop("Make sure that the first entry in speciesUpdateFunction is the default expression")
  }

  for (fn in P(sim)$speciesUpdateFunction) {
    if (is(fn, "call")) {
      sim$species <- eval(fn)
    } else {
      stop("speciesUpdateFunction should be a list of functions.")
    }
  }

  if (getOption("LandR.verbose") > 0) {
    message("Adjusting species-level traits, part 2, for LandWeb")
    print(sim$species)
  }

  ## check that all species have trait values.
  missingTraits <- setdiff(names(sim$speciesLayers), sim$species$species)
  if (length(missingTraits) == length(names(sim$speciesLayers))) {
    stop("No trait values where found for ", paste(missingTraits, collapse = ", "), ".\n",
         "Please check the species list and traits table")
  } else if (length(missingTraits))
    warning("No trait values where found for ", paste(missingTraits, collapse = ", "), ".\n",
            "Please check the species list and traits table")

  ### make table of light shade tolerance  #######################
  sim$sufficientLight <- data.frame(speciesshadetolerance = 1:5,
                                    X0 = 1,
                                    X1 = c(0.5, rep(1, 4)),
                                    X2 = c(0, 0.5, rep(1, 3)),
                                    X3 = c(rep(0, 2), 0.5, rep(1, 2)),
                                    X4 = c(rep(0, 3), 0.5, 1),
                                    X5 = c(rep(0, 4), 1))

  ################################################################
  ## initialEcoregionMap
  ################################################################
  if (!compareCRS(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  if (!compareCRS(crs(sim$studyAreaLarge), crs(sim$rasterToMatchLarge))) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  rstLCCAdj <- sim$rstLCC

  ## Clean pixels for veg. succession model
  ## remove pixes with no spp data
  pixelsToRm <- is.na(sim$speciesLayers[[1]][])
  pixelFateDT <- pixelFate(fate = "Total number pixels", runningPixelTotal = ncell(sim$speciesLayers))
  pixelFateDT <- pixelFate(pixelFateDT, "NAs on sim$speciesLayers", sum(pixelsToRm))

  ## remove non-forested if asked by user
  if (P(sim)$omitNonTreedPixels) {
    if (is.null(P(sim)$forestedLCCClasses))
      stop("No P(sim)$forestedLCCClasses provided, but P(sim)$omitNonTreedPixels is TRUE.
           \nPlease provide a vector of forested classes in P(sim)$forestedLCCClasses")

    lccPixelsRemoveTF <- !(sim$rstLCC[] %in% P(sim)$forestedLCCClasses)
    pixelsToRm <- lccPixelsRemoveTF | pixelsToRm
    pixelFateDT <- pixelFate(pixelFateDT, "Non forested pixels (based on LCC classes)",
                             sum(lccPixelsRemoveTF) - tail(pixelFateDT$pixelsRemoved, 1))
  }

  rstLCCAdj[pixelsToRm] <- NA

  # The next function will remove the "zero" class on sim$ecoregionRst
  pixelFateDT <- pixelFate(pixelFateDT, "Removing 0 class in sim$ecoregionRst",
                           sum(sim$ecoregionRst[][!pixelsToRm] == 0, na.rm = TRUE))

  ecoregionFiles <- Cache(prepEcoregions,
                          ecoregionRst = sim$ecoregionRst,
                          ecoregionLayer = sim$ecoregionLayer,
                          ecoregionLayerField = P(sim)$ecoregionLayerField,
                          rasterToMatchLarge = sim$rasterToMatchLarge,
                          rstLCCAdj = rstLCCAdj,
                          pixelsToRm = pixelsToRm,
                          cacheTags = cacheTags)

  ################################################################
  ## put together pixelTable object
  ################################################################
  #  Round age to pixelGroupAgeClass
  # Internal data.table is changed; using memoise here causes the internal changes to
  #   come out to the pixelTable, which is not desired. Turn off memoising for one step
  opts <- options("reproducible.useMemoise" = FALSE)
  on.exit({
    options(opts)
  }, add = TRUE)
  pixelTable <- Cache(makePixelTable,
                      speciesLayers = sim$speciesLayers,
                      species = sim$species,
                      standAgeMap = sim$standAgeMap,
                      ecoregionFiles = ecoregionFiles,
                      biomassMap = sim$rawBiomassMap,
                      rasterToMatch = sim$rasterToMatchLarge,
                      rstLCC = rstLCCAdj,
                      # pixelGroupAgeClass = P(sim)$pixelGroupAgeClass,
                      userTags = c(cacheTags, "pixelTable"),
                      omitArgs = c("userTags"))
  options(opts)
  #######################################################
  # Make the initial pixelCohortData table
  #######################################################
  coverColNames <- paste0("cover.", sim$species$species)
  pixelCohortData <- Cache(makeAndCleanInitialCohortData, pixelTable,
                           sppColumns = coverColNames,
                           imputeBadAgeModel = P(sim)$imputeBadAgeModel,
                           #pixelGroupBiomassClass = P(sim)$pixelGroupBiomassClass,
                           #pixelGroupAgeClass = P(sim)$pixelGroupAgeClass,
                           minCoverThreshold = P(sim)$minCoverThreshold,
                           doSubset = P(sim)$subsetDataAgeModel,
                           userTags = c(cacheTags, "pixelCohortData"),
                           omitArgs = c("userTags"))

  pixelFateDT <- pixelFate(pixelFateDT, "makeAndCleanInitialCohortData rm cover < minThreshold",
                           tail(pixelFateDT$runningPixelTotal,1) - NROW(unique(pixelCohortData$pixelIndex)))
  #######################################################
  # Partition totalBiomass into individual species B, via estimating how %cover and %biomass
  #   are related
  #######################################################
  message(blue("Partitioning totalBiomass per pixel into cohort B as:"))
  if (isTRUE(P(sim)$fitDecidiousCoverDiscount)) {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))

    pixelCohortData[, lcc := as.factor(lcc)]

    plot.it <- FALSE
    sam <- subsetDT(pixelCohortData, by = c("speciesCode", "lcc"),
                    doSubset = P(sim)$subsetDataAgeModel,
                    indices = TRUE)
    pi <- unique(pixelCohortData[sam]$pixelIndex)
    sam <- which(pixelCohortData$pixelIndex %in% pi)

    system.time(out <- optimize(interval = c(0.1, 1), f = coverOptimFn, bm = P(sim)$coverPctToBiomassPctModel,
                                pixelCohortData = pixelCohortData,
                                subset = sam, maximum = FALSE))
    deciduousCoverDiscount <- out$minimum
    if (plot.it) {
      coverModel <- coverOptimFn(out$minimum, pixelCohortData, P(sim)$subsetDataAgeModel,
                                 P(sim)$coverPctToBiomassPctModel,
                                 returnRsq = FALSE)
      sam1 <- sample(NROW(pixelCohortData), 1e5)
      dev()
      par(mfrow = c(1,2))
      plot(predict(coverModel$modelBiomass1$mod, newdata = coverModel$pixelCohortData[sam1]),
           log(coverModel$pixelCohortData$B/100)[sam1], pch = ".")
      abline(a = 0, b = 1)

      coverModel1 <- coverOptimFn(1, pixelCohortData, P(sim)$subsetDataAgeModel,
                                  P(sim)$coverPctToBiomassPctModel,
                                  returnRsq = FALSE)
      dev()
      plot(predict(coverModel1$modelBiomass1$mod, newdata = coverModel1$pixelCohortData[sam1]),
           log(coverModel1$pixelCohortData$B/100)[sam1], pch = ".")
      abline(a = 0, b = 1)

      pcd <- pixelCohortData
      bb <- pcd[sample(sam)]
      cc <- bb[,cover3:=cover*c(1,out$minimum)[decid+1]][
        ,actualX:=cover3/sum(cover3)/(cover/100), by = "pixelIndex"]
      setkey(cc, pixelIndex)
      mean(cc[speciesCode == "Popu_Tre"]$actualX)
    }

  } else {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))
    deciduousCoverDiscount <- P(sim)$deciduousCoverDiscount
    message(blue("using previously estimated decidiousCoverDiscount:", round(deciduousCoverDiscount,3)))
  }
  pixelCohortData <- partitionBiomass(x = deciduousCoverDiscount, pixelCohortData)
  set(pixelCohortData, NULL, "B", asInteger(pixelCohortData$B/P(sim)$pixelGroupBiomassClass)*
        P(sim)$pixelGroupBiomassClass)
  set(pixelCohortData, NULL, c("decid", "cover2"), NULL)
  set(pixelCohortData, NULL, "cover", asInteger(pixelCohortData$cover))

  ####################################################### replace 34 and 35 and
  #36 values -- burns and cities -- to a neighbour class *that exists*.
  # 1. We need
  #to have a spatial estimate of maxBiomass everywhere there is forest; we can't
  #have gaps The pixels that are 34, 35 or 36 are places for which we don't want
  #maxBiomass associated with their LCC ... i.e., we don't want a maximum
  #biomass associated with 34 and 35 because those classes are transient. They
  #will transition to another class before they arrive at a tree maximum
  #biomass. So, 34 and 35 should not have estimates of maxBiomass 36 is urban.
  #So, we SHOULD remove these pixels from our studies, except if we are doing
  #NRV studies (e.g., LandWeb wanted to replace 36 with some forest class) We
  #decided that we should not use 34 or 35 in our models of Biomass because the
  #only objective of these models is to estimate maxBiomass, so don't use 34 or
  #35 To associate the pixels that were 34 or 35 with a maxBiomass , we need to
  #give them a "forest class" that they might "become" after they grow out of
  #being 34 or 35. The pixels where there were 34 and 35 nevertheless have
  #Biomass estimates in them from KNN and other sources. We leave those as is.
  #######################################################
  uwc <- P(sim)$LCCClassesToReplaceNN

  message("Replace ", paste(uwc, collapse = ", "),
          " values -- ", "burns"[any(uwc %in% 34:35)], " and cities"[any(uwc %in% 36)],
          " -- to a neighbour class *that exists*")

  rmZeroBiomassQuote <- quote(totalBiomass > 0)
  ## version 1: from before March 2019 - Ceres noticed it created issues with fitting modelCover
  ## March 2020: seems to be the preferred behaviour?
  availableCombinations <- unique(pixelCohortData[eval(rmZeroBiomassQuote),
                                                  .(speciesCode, initialEcoregionCode, pixelIndex)])
  ## version 2: Ceres's fix from March 2019 to solve issues with modelCover fitting (?)
  # availableCombinations <- unique(pixelCohortData[,
  #                                                 .(speciesCode, initialEcoregionCode, pixelIndex)])
  ## version 3: Feb 2020 Eliot's fix that is WRONG - this behaviour is being achieved in convertUnwantedLCC and creates empty tables if done here
  # availableCombinations <- unique(pixelCohortData[!(lcc %in% uwc),
  #                                                 .(speciesCode, initialEcoregionCode, pixelIndex)])
  newLCCClasses <- Cache(convertUnwantedLCC,
                         classesToReplace = P(sim)$LCCClassesToReplaceNN,
                         rstLCC = rstLCCAdj,
                         availableERC_by_Sp = availableCombinations,
                         userTags = c(cacheTags, "newLCCClasses", "stable"),
                         omitArgs = c("userTags"))

  ## split pixelCohortData into 2 parts -- one with the former 34:36 pixels, one without
  #    The one without 34:36 can be used for statistical estimation, but not the one with
  cohortData34to36 <- pixelCohortData[pixelIndex %in% newLCCClasses$pixelIndex]
  cohortData34to36 <- merge(newLCCClasses, cohortData34to36, all.x = TRUE,
                            all.y = FALSE, by = "pixelIndex")
  cohortDataNo34to36 <- pixelCohortData[!pixelIndex %in% newLCCClasses$pixelIndex]
  setnames(cohortDataNo34to36, "initialEcoregionCode", "ecoregionGroup")
  cohortDataNo34to36NoBiomass <- cohortDataNo34to36[eval(rmZeroBiomassQuote),
                                                    .(B, logAge, speciesCode, ecoregionGroup, lcc, cover)]
  cohortDataNo34to36NoBiomass <- unique(cohortDataNo34to36NoBiomass)

  ## make sure ecoregionGroups match
  assert1(cohortData34to36, pixelCohortData, rmZeroBiomassQuote)

  ##############################################################
  # Statistical estimation of establishprob, maxB and maxANPP
  ##############################################################
  cohortDataShort <- cohortDataNo34to36[, list(coverPres = sum(cover > 0)),
                                        by = c("ecoregionGroup", "speciesCode")]
  # find coverNum for each known class
  aa <- table(pixelTable$initialEcoregionCode)
  dt1 <- data.table(ecoregionGroup = factor(names(aa)), coverNum = as.integer(unname(aa)))
  allCombos <- expand.grid(ecoregionGroup = dt1$ecoregionGroup, speciesCode = unique(cohortDataShort$speciesCode))
  setDT(allCombos)
  dt1 <- dt1[allCombos, on = "ecoregionGroup", nomatch = 0]
  cohortDataShortNoCover <- cohortDataShort[dt1, on = c("ecoregionGroup", "speciesCode"), nomatch = NA]

  #cohortDataShortNoCover <- cohortDataShort[coverPres == 0] #
  cohortDataShort <- cohortDataShortNoCover[coverPres > 0] # remove places where there is 0 cover
  cohortDataShortNoCover <- cohortDataShortNoCover[is.na(coverPres)][, coverPres := 0]
  # will be added back as establishprob = 0
  message(blue("Estimating Species Establishment Probability using P(sim)$coverModel, which is"))
  message(magenta(paste0(format(P(sim)$coverModel, appendLF = FALSE), collapse = "")))

  # for backwards compatibility -- change from parameter to object
  if (is.null(sim$cloudFolderID))
    if (!is.null(P(sim)$cloudFolderID))
      sim$cloudFolderID <- P(sim)$cloudFolderID

  useCloud <- if (!is.null(sim$cloudFolderID)) {
    (getOption("reproducible.useCache", FALSE) && P(sim)$useCloudCacheForStats)
  } else {
    FALSE
  }

  # Remove all cases where there is 100% presence in an ecoregionGroup -- causes failures in binomial models
  cdsWh <- cohortDataShort$coverPres == cohortDataShort$coverNum
  cds <- Copy(cohortDataShort)
  cds <- cds[!cdsWh]
  modelCover <- Cache(statsModel,
                      modelFn = P(sim)$coverModel,
                      # modelFn = cm,
                      uniqueEcoregionGroup = .sortDotsUnderscoreFirst(as.character(unique(cohortDataShort$ecoregionGroup))),
                      sumResponse = sum(cohortDataShort$coverPres, cohortDataShort$coverNum, na.rm = TRUE),
                      .specialData = cds,
                      useCloud = useCloud,
                      cloudFolderID = sim$cloudFolderID,
                      # useCache = "overwrite",
                      showSimilar = getOption("reproducible.showSimilar", FALSE),
                      userTags = c(cacheTags, "modelCover"),
                      omitArgs = c("showSimilar", "useCache", ".specialData", "useCloud", "cloudFolderID"))
  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelCover$rsq, 4))), function(x) message(blue(x)))
  if (isTRUE(any(cdsWh))) {
    cds[, pred := fitted(modelCover$mod, response = "response")]
    cohortDataShort <- cds[, -c("coverPres", "coverNum")][cohortDataShort,
                                                          on = c('ecoregionGroup', 'speciesCode'), nomatch = NA]
    cohortDataShort[is.na(pred), pred := 1]
    modelCover <- cohortDataShort$pred
  }

  ## For biomass
  ### Subsample cases where there are more than 50 points in an ecoregionGroup * speciesCode
  totalBiomass <- sum(cohortDataNo34to36NoBiomass$B, na.rm = TRUE)
  cohortDataNo34to36NoBiomass <- subsetDT(cohortDataNo34to36NoBiomass,
                                          by = c("ecoregionGroup", "speciesCode"),
                                          doSubset = P(sim)$subsetDataBiomassModel)

  ### For Cache -- doesn't need to cache all columns in the data.table -- only the ones in the model
  ### force parameter values to avoid more checks
  # If using mixed effect model, see here for good discussion of
  #  shrinkage https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
  message(blue("Estimating biomass using P(sim)$biomassModel as:\n"),
          magenta(paste0(format(P(sim)$biomassModel, appendLF = FALSE), collapse = "")))
  modelBiomass <- Cache(
    statsModel,
    modelFn = P(sim)$biomassModel,
    uniqueEcoregionGroup = .sortDotsUnderscoreFirst(as.character(unique(cohortDataNo34to36NoBiomass$ecoregionGroup))),
    sumResponse = totalBiomass,
    .specialData = cohortDataNo34to36NoBiomass,
    useCloud = useCloud,
    # useCache = "overwrite",
    cloudFolderID = sim$cloudFolderID,
    showSimilar = getOption("reproducible.showSimilar", FALSE),
    userTags = c(cacheTags, "modelBiomass", paste0("subsetSize:", P(sim)$subsetDataBiomassModel)),
    omitArgs = c("showSimilar", ".specialData", "useCloud", "cloudFolderID", "useCache")
  )

  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelBiomass$rsq, 4))), function(x) message(blue(x)))

  ########################################################################
  # create speciesEcoregion -- a single line for each combination of ecoregionGroup & speciesCode
  #   doesn't include combinations with B = 0 because those places can't have the species/ecoregion combo
  ########################################################################
  message(blue("Create speciesEcoregion using modelCover and modelBiomass to estimate species traits"))
  speciesEcoregion <- makeSpeciesEcoregion(cohortDataNoBiomass = cohortDataNo34to36NoBiomass,
                                           cohortDataShort = cohortDataShort,
                                           cohortDataShortNoCover = cohortDataShortNoCover,
                                           species = sim$species,
                                           modelCover = modelCover,
                                           modelBiomass = modelBiomass,
                                           successionTimestep = P(sim)$successionTimestep,
                                           currentYear = time(sim))

  #######################################
  if (!is.na(P(sim)$.plotInitialTime)) {
    uniqueSpeciesNames <- as.character(unique(speciesEcoregion$speciesCode))
    names(uniqueSpeciesNames) <- uniqueSpeciesNames
    speciesEcoregionTable2 <- copy(speciesEcoregion)
    speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionGroup)]
    maxB <- raster::stack(lapply(uniqueSpeciesNames, function(sp) {
      rasterizeReduced(speciesEcoregionTable2[speciesCode == sp], ecoregionFiles$ecoregionMap,
                       "maxB", "ecoregionInt")
    }))
    curDev <- dev.cur()
    quickPlot::dev(6, width = 18, height = 10)
    Plot(maxB, legendRange = c(0, max(maxValue(maxB))))
    quickPlot::dev(curDev)
  }

  if (ncell(sim$rasterToMatchLarge) > 3e6) .gc()

  ########################################################################
  # Create initial communities, i.e., pixelGroups
  ########################################################################
  # Rejoin back the pixels that were 34 and 35
  set(cohortData34to36, NULL, "initialEcoregionCode", NULL)
  pixelCohortData <- rbindlist(list(cohortData34to36, cohortDataNo34to36),
                               use.names = TRUE, fill = TRUE)

  ########################################################################
  # "Downsize" to studyArea after estimating parameters on studyAreaLarge
  ########################################################################
  ## 1. Subset pixels (IDs) on rasterToMatchLarge, using rasterToMatch
  ## 2. Subset data.tables using the pixel IDs / ecoregion/species combinations
  ##    that are common across the two rasters
  ## 3. Re-do pixel ID numbering so that it matches the final rasterToMatch
  ## Note: if SA and SALarge are the same, no subsetting will take place.

  if (!identical(extent(sim$rasterToMatch), extent(sim$rasterToMatchLarge))) {
    message(blue("Subsetting to studyArea"))
    rasterToMatchLarge <- sim$rasterToMatchLarge
    rasterToMatchLarge <- setValues(rasterToMatchLarge, seq(ncell(rasterToMatchLarge)))
    rasterToMatchLarge <- Cache(postProcess,
                                x = rasterToMatchLarge,
                                rasterToMatch = sim$rasterToMatch,
                                maskWithRTM = TRUE,
                                filename2 = NULL,
                                userTags = c(cacheTags, "rasterToMatchLarge"),
                                omitArgs = c("userTags"))

    if (!compareRaster(rasterToMatchLarge, sim$rasterToMatch, orig = TRUE, stopiffalse = FALSE))
      stop("Downsizing to rasterToMatch after estimating parameters didn't work.
           Please debug Biomass_borealDataPrep::createBiomass_coreInputs()")

    ## subset pixels that are in studyArea/rasterToMatch only
    pixToKeep <- na.omit(getValues(rasterToMatchLarge))
    pixelCohortData <- pixelCohortData[pixelIndex %in% pixToKeep]

    # re-do pixelIndex (it now needs to match rasterToMatch)
    newPixelIndexDT <- data.table(pixelIndex = getValues(rasterToMatchLarge),
                                  newPixelIndex = as.integer(1:ncell(rasterToMatchLarge)))

    pixelCohortData <- newPixelIndexDT[pixelCohortData, on = "pixelIndex"]
    pixelCohortData[, pixelIndex := NULL]
    setnames(pixelCohortData, old = "newPixelIndex", new = "pixelIndex")
    rm(rasterToMatchLarge)
  }

  if (ncell(sim$rasterToMatch) > 3e6) .gc()

  ## subset ecoregionFiles$ecoregionMap to smaller area.
  ecoregionFiles$ecoregionMap <- Cache(postProcess,
                                       x = ecoregionFiles$ecoregionMap,
                                       rasterToMatch = sim$rasterToMatch,
                                       maskWithRTM = TRUE,
                                       filename2 = NULL,
                                       userTags = c(cacheTags, "ecoregionMap"),
                                       omitArgs = c("userTags"))

  maxAgeHighQualityData <- -1
  if (length(extractURL("fireURL"))) {
    firstFireYear = as.numeric(gsub("^.+nbac_(.*)_to.*$", "\\1", extractURL("fireURL")))
    if (!is.na(firstFireYear)) {
      maxAgeHighQualityData <- start(sim) - firstFireYear
      youngRows <- pixelCohortData$age <= maxAgeHighQualityData
      young <- pixelCohortData[youngRows == TRUE]
      # whYoungBEqZero <- which(young$B == 0)
      whYoungAgeEqZero <- which(young$age == 0)
      if (length(whYoungAgeEqZero)) {
        youngWAgeEqZero <- young[whYoungAgeEqZero]
        youngNoAgeEqZero <- young[-whYoungAgeEqZero]
      }
      young <- Cache(updateYoungBiomasses, youngNoAgeEqZero,
                                    biomassModel = modelBiomass$mod)
      set(young, NULL, setdiff(colnames(young), colnames(pixelCohortData)), NULL)

      # put the B = 0
      if (length(whYoungAgeEqZero)) {
        young <- rbindlist(list(young, youngWAgeEqZero), use.names = TRUE)
      }
      pixelCohortData <- rbindlist(list(pixelCohortData[youngRows == FALSE],
                                        young), use.names = TRUE)
    }
  }

  ## make cohortDataFiles: pixelCohortData (rm unnecessary cols, subset pixels with B>0,
  ## generate pixelGroups, add ecoregionGroup and totalBiomass) and cohortData
  cohortDataFiles <- makeCohortDataFiles(pixelCohortData, columnsForPixelGroups, speciesEcoregion,
                                         pixelGroupBiomassClass = P(sim)$pixelGroupBiomassClass,
                                         pixelGroupAgeClass = P(sim)$pixelGroupAgeClass,
                                         minAgeForGrouping = maxAgeHighQualityData,
                                         pixelFateDT = pixelFateDT
  )

  sim$cohortData <- cohortDataFiles$cohortData
  pixelCohortData <- cohortDataFiles$pixelCohortData
  pixelFateDT <- cohortDataFiles$pixelFateDT

  rm(cohortDataFiles)

  ## make a table of available active and inactive (no biomass) ecoregions
  sim$ecoregion <- makeEcoregionDT(pixelCohortData, speciesEcoregion)
  #sim$ecoregion <- ecoregionFiles$ecoregion ## TODO: don't use this one yet (ever?)

  ## make biomassMap, ecoregionMap, minRelativeB, pixelGroupMap (at the scale of rasterToMatch)
  sim$biomassMap <- makeBiomassMap(pixelCohortData, sim$rasterToMatch)
  sim$ecoregionMap <- makeEcoregionMap(ecoregionFiles, pixelCohortData)
  sim$minRelativeB <- makeMinRelativeB(pixelCohortData)
  sim$pixelGroupMap <- makePixelGroupMap(pixelCohortData, sim$rasterToMatch)

  ## make sure speciesLayers match RTM (since that's what is used downstream in simulations)
  ## TODO: use postProcess?

  message(blue("Writing sim$speciesLayers to disk as they are likely no longer needed in RAM"))
  sim$speciesLayers <- Cache(postProcess, sim$speciesLayers,
                             rasterToMatch = sim$rasterToMatch,
                             maskWithRTM = TRUE, filename2 = paste0(names(sim$speciesLayers), ".tif"))
  # if (!compareRaster(sim$speciesLayers, sim$rasterToMatch, stopiffalse = FALSE))
  #   sim$speciesLayers <- cropInputs(sim$speciesLayers, sim$rasterToMatch)
  # sim$speciesLayers <- Cache(maskInputs, sim$speciesLayers,
  #                                 rasterToMatch = sim$rasterToMatch,
  #                                 maskWithRTM = TRUE)

  ## double check these rasters all match RTM
  compareRaster(sim$biomassMap, sim$ecoregionMap, sim$pixelGroupMap, sim$rasterToMatch, sim$speciesLayers)

  ## rm ecoregions that may not be present in rasterToMatch
  ## make ecoregionGroup a factor and export speciesEcoregion to sim
  speciesEcoregion <- speciesEcoregion[ecoregionGroup %in% pixelCohortData$ecoregionGroup]
  speciesEcoregion[, ecoregionGroup := factor(as.character(ecoregionGroup))]
  sim$speciesEcoregion <- speciesEcoregion

  ## write species layers to disk
  # sim$speciesLayers <- lapply(seq(numLayers(sim$speciesLayers)), function(x) {
  #   writeRaster(sim$speciesLayers[[x]],
  #               file.path(outputPath(sim), paste0(names(sim$speciesLayers)[x], ".tif")),
  #               datatype = "INT2U", overwrite = TRUE)
  # }) %>% raster::stack()

  ## do assertions
  message(blue("Create pixelGroups based on: ", paste(columnsForPixelGroups, collapse = ", "),
               "\n  Resulted in", magenta(length(unique(sim$cohortData$pixelGroup))),
               "unique pixelGroup values"))
  LandR::assertERGs(sim$ecoregionMap, cohortData = sim$cohortData,
                    speciesEcoregion = sim$speciesEcoregion,
                    minRelativeB = sim$minRelativeB)

  LandR::assertCohortData(sim$cohortData, sim$pixelGroupMap)

  message("Done Biomass_borealDataPrep: ", Sys.time())
  sim$pixelFateDT <- pixelFateDT
  out <- messageDF(pixelFateDT, 3, "blue")
  #out <- lapply(capture.output(sim$pixelFateDT), function(x) message(blue(x)))

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
  whThisMod <- which(unlist(lapply(a@dependencies, function(x) x@name)) == "Biomass_borealDataPrep")
  objNames <- a@dependencies[[whThisMod]]@inputObjects$objectName
  objExists <- !unlist(lapply(objNames, function(x) is.null(sim[[x]])))
  names(objExists) <- objNames

  # Filenames
  lcc2005Filename <- file.path(dPath, "LCC2005_V1_4a.tif")

  ## Study area(s) ------------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
    sim <- objectSynonyms(sim, list(c("studyAreaLarge", "studyArea")))
  }

  if (!identical(crs(sim$studyArea), crs(sim$studyAreaLarge))) {
    warning("studyArea and studyAreaLarge have different projections.\n
            studyAreaLarge will be projected to match crs(studyArea)")
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$studyArea))
  }

  ## check whether SA is within SALarge
  ## convert to temp sf objects
  studyArea <- st_as_sf(sim$studyArea)
  studyAreaLarge <- st_as_sf(sim$studyAreaLarge)

  #this is necessary if studyArea and studyAreaLarge are multipolygon objects
  if (nrow(studyArea) > 1) {
    studyArea <- st_union(studyArea) %>%
      st_as_sf(.)
  }

  if (nrow(studyAreaLarge) > 1) {
    studyAreaLarge <- st_union(studyArea) %>%
      st_as_sf(.)
  }

  if (length(st_within(studyArea, studyAreaLarge))[[1]] == 0)
    stop("studyArea is not fully within studyAreaLarge.
         Please check the aligment, projection and shapes of these polygons")
  rm(studyArea, studyAreaLarge)

  ## Raster(s) to match ------------------------------------------------
  needRTM <- FALSE
  if (is.null(sim$rasterToMatch) || is.null(sim$rasterToMatchLarge)) {
    if (!suppliedElsewhere("rasterToMatch", sim) ||
        !suppliedElsewhere("rasterToMatchLarge", sim)) {      ## if one is not provided, re do both (safer?)
      needRTM <- TRUE
      message("There is no rasterToMatch/rasterToMatchLarge supplied; will attempt to use rawBiomassMap")
    } else {
      stop("rasterToMatch/rasterToMatchLarge is going to be supplied, but ", currentModule(sim), " requires it ",
           "as part of its .inputObjects. Please make it accessible to ", currentModule(sim),
           " in the .inputObjects by passing it in as an object in simInit(objects = list(rasterToMatch = aRaster)",
           " or in a module that gets loaded prior to ", currentModule(sim))
    }
  }

  if (!suppliedElsewhere("rawBiomassMap", sim) || needRTM) {
    sim$rawBiomassMap <- Cache(prepInputs,
                               url = extractURL("rawBiomassMap"),
                               destinationPath = dPath,
                               studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                               rasterToMatch = if (!needRTM) sim$rasterToMatchLarge else NULL,
                               maskWithRTM = if (!needRTM) TRUE else FALSE,
                               useSAcrs = FALSE,     ## never use SA CRS
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = TRUE, overwrite = TRUE,
                               userTags = c(cacheTags, "rawBiomassMap"),
                               omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
  }

  if (needRTM) {
    ## if we need rasterToMatch/rasterToMatchLarge, that means a) we don't have it, but b) we will have rawBiomassMap
    ## even if one of the rasterToMatch is present re-do both.

    if (is.null(sim$rasterToMatch) != is.null(sim$rasterToMatchLarge))
      warning(paste0("One of rasterToMatch/rasterToMatchLarge is missing. Both will be created \n",
                     "from rawBiomassMap and studyArea/studyAreaLarge.\n
                     If this is wrong, provide both rasters"))

    sim$rasterToMatchLarge <- sim$rawBiomassMap
    RTMvals <- getValues(sim$rasterToMatchLarge)
    sim$rasterToMatchLarge[!is.na(RTMvals)] <- 1

    sim$rasterToMatchLarge <- Cache(writeOutputs, sim$rasterToMatchLarge,
                                    filename2 = file.path(cachePath(sim), "rasters", "rasterToMatchLarge.tif"),
                                    datatype = "INT2U", overwrite = TRUE,
                                    userTags = c(cacheTags, "rasterToMatchLarge"),
                                    omitArgs = c("userTags"))

    sim$rasterToMatch <- Cache(postProcess,
                               x = sim$rawBiomassMap,
                               studyArea = sim$studyArea,
                               rasterToMatch = sim$rasterToMatchLarge,
                               useSAcrs = FALSE,
                               maskWithRTM = FALSE,   ## mask with SA
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = file.path(cachePath(sim), "rasterToMatch.tif"),
                               overwrite = TRUE,
                               userTags = c(cacheTags, "rasterToMatch"),
                               omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))

    ## covert to 'mask'
    RTMvals <- getValues(sim$rasterToMatch)
    sim$rasterToMatch[!is.na(RTMvals)] <- 1
  }

  ## if using custom raster resolution, need to allocate biomass proportionally to each pixel
  ## if no rawBiomassMap/RTM/RTMLarge were suppliedElsewhere, the "original" pixel size respects
  ## whatever resolution comes with the rawBiomassMap data
  simPixelSize <- unique(asInteger(res(sim$rasterToMatchLarge)))
  origPixelSize <- 250L # unique(res(sim$rawBiomassMap)) ## TODO: figure out a good way to not hardcode this

  if (simPixelSize != origPixelSize) { ## make sure we are comparing integers, else else %!=%
    rescaleFactor <- (origPixelSize / simPixelSize)^2
    sim$rawBiomassMap <- sim$rawBiomassMap / rescaleFactor
  }

  # if (ncell(sim$rasterToMatch) < 1e4)
  # stop("sim$rasterToMatch is too small, it should have more than 10,000 pixels")

  ## TODO: KEEP THIS HERE OR ONLY INIT?
  if (!identical(crs(sim$studyArea), crs(sim$rasterToMatch))) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  if (!identical(crs(sim$studyAreaLarge), crs(sim$rasterToMatchLarge))) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  ## Land cover raster ------------------------------------------------
  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- Cache(prepInputs,
                        targetFile = lcc2005Filename,
                        archive = asPath("LandCoverOfCanada2005_V1_4.zip"),
                        url = extractURL("rstLCC"),
                        destinationPath = dPath,
                        studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                        # studyArea = sim$studyArea,
                        rasterToMatch = sim$rasterToMatchLarge,
                        # rasterToMatch = sim$rasterToMatch,
                        maskWithRTM = TRUE,
                        method = "bilinear",
                        datatype = "INT2U",
                        filename2 = TRUE, overwrite = TRUE,
                        userTags = c("prepInputsrstLCC_rtm", currentModule(sim)), # use at least 1 unique userTag
                        omitArgs = c("destinationPath", "targetFile", "userTags"))

    if (!identical(projection(sim$rstLCC),
                   projection(sim$rasterToMatchLarge)))
      projection(sim$rstLCC) <- projection(sim$rasterToMatchLarge) ## Ceres: this shouldn't be necessary anymore
  }

  ## Ecodistrict ------------------------------------------------
  if (!suppliedElsewhere("ecoregionLayer", sim)) {
    sim$ecoregionLayer <- Cache(prepInputs,
                                targetFile = "ecodistricts.shp",
                                archive = asPath("ecodistrict_shp.zip"),
                                url = extractURL("ecoregionLayer", sim),
                                alsoExtract = "similar",
                                destinationPath = dPath,
                                studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                                overwrite = TRUE,
                                useSAcrs = TRUE, # this is required to make ecoZone be in CRS of studyArea
                                fun = "raster::shapefile",
                                userTags = c("prepInputsEcoDistrict_SA", currentModule(sim), cacheTags))
  }

  ## Stand age map ------------------------------------------------
  if (!suppliedElsewhere("standAgeMap", sim)) {
    sim$standAgeMap <- Cache(prepInputsStandAgeMap,
                             destinationPath = dPath,
                             ageURL = extractURL("standAgeMap"),
                             ageFun = "raster::raster",
                             studyArea = raster::aggregate(sim$studyAreaLarge),
                             #studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                             rasterToMatch = sim$rasterToMatchLarge,
                             # rasterToMatch = sim$rasterToMatch,
                             maskWithRTM = TRUE,
                             method = "bilinear",
                             datatype = "INT2U",
                             filename2 = NULL,
                             overwrite = TRUE,
                             fireURL = extractURL("fireURL"),
                             fireFun = "sf::st_read",
                             fireField = "YEAR",
                             startTime = start(sim),
                             userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
                             omitArgs = c("destinationPath", "targetFile", "overwrite",
                                          "alsoExtract", "userTags"))
  }

  ## Species equivalencies table -------------------------------------------

  if (!suppliedElsewhere("sppEquiv", sim)) {
    if (!is.null(sim$sppColorVect))
      message("No 'sppColorVect' provided; using default colour palette: Accent")

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
    if (is.null(sim$sppColorVect)) {
      ## add default colors for species used in model
      sim$sppColorVect <- sppColors(sim$sppEquiv, P(sim)$sppEquivCol,
                                    newVals = "Mixed", palette = "Accent")
      message("No 'sppColorVect' provided; using default colour palette: Accent")
    }
  }

  ## Species raster layers -------------------------------------------
  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opts <- options(reproducible.useCache = "overwrite")
    sim$speciesLayers <- Cache(loadkNNSpeciesLayers,
                               dPath = dPath,
                               rasterToMatch = sim$rasterToMatchLarge,
                               # rasterToMatch = sim$rasterToMatch,
                               studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                               sppEquiv = sim$sppEquiv,
                               knnNamesCol = "KNN",
                               sppNameVector = sim$sppNameVector,
                               sppEquivCol = P(sim)$sppEquivCol,
                               thresh = 10,
                               url = extractURL("speciesLayers"),
                               userTags = c(cacheTags, "speciesLayers"),
                               omitArgs = c("userTags"))
  }

  # 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath = dPath,
                                        cacheTags = c(cacheTags, "speciesTable"))
  }

  if (!suppliedElsewhere("columnsForPixelGroups", sim)) {
    sim$columnsForPixelGroups <- LandR::columnsForPixelGroups
  }

  return(invisible(sim))
}


#' @importFrom sf st_cast st_transform
#' @importFrom fasterize fasterize
#' @importFrom raster crs
prepInputsFireYear <- function(..., rasterToMatch, field) {
  a <- Cache(prepInputs, ...)
  gg <- st_cast(a, "MULTIPOLYGON") # collapse them into a single multipolygon
  d <- st_transform(gg, crs(rasterToMatch))
  fasterize(d, raster = rasterToMatch, field = field)

}

prepInputsStandAgeMap <- function(..., ageURL, ageFun, maskWithRTM,
                                  method, datatype, filename2,
                                  fireURL, fireFun,
                                  rasterToMatch, fireField, startTime) {
  standAgeMap <- Cache(prepInputs, ...,
                       maskWithRTM = maskWithRTM, method = method,
                       datatype = datatype, filename2 = filename2,
                       url = ageURL, fun = ageFun, rasterToMatch = rasterToMatch)
  standAgeMap[] <- asInteger(standAgeMap[])
  if (!(missing(fireURL) || is.null(fireURL) || is.na(fireURL))) {
    fireYear <- Cache(prepInputsFireYear, ...,
                      url = fireURL,
                      fun = fireFun,
                      rasterToMatch = rasterToMatch,
                      field = fireField
    )
    toChange <- !is.na(fireYear[]) & fireYear[] <= asInteger(startTime)
    standAgeMap[] <- asInteger(standAgeMap[])
    standAgeMap[toChange] <- asInteger(startTime) - asInteger(fireYear[][toChange])
  }
  standAgeMap

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

updateYoungBiomasses <- function(young, biomassModel) {
  if (is(biomassModel, "merMod")) {
    columns <- c("ecoregionGroup", "logAge", "speciesCode", "cover")
    young2 <- unique(young, by = columns)
    message(green("  -- Calculating bootstrap estimates around B; will replace B in young data if it is beyond 95% CI"))
    message(green("     This will take a bit"))
    reproducible::Require("merTools")
    PI.time <- system.time(
      PI <- merTools::predictInterval(merMod = biomassModel, newdata = young2,
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
