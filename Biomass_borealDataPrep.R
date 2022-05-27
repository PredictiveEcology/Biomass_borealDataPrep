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
  version = list(Biomass_borealDataPrep = "1.5.0.9000"),
  spatialExtent = raster::extent(rep(NA_real_, 4)),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_borealDataPrep.Rmd"),
  reqdPkgs = list("assertthat", "crayon", "data.table", "dplyr", "fasterize", "plyr", "raster",
                  "sp", "sf", "merTools", "SpaDES.tools",
                  "PredictiveEcology/reproducible@development (>=1.1.1.9004)",
                  "achubaty/amc@development (>=0.1.6.9000)",
                  "PredictiveEcology/LandR@development (>=0.0.11.9008)",
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    defineParameter("biomassModel", "call",
                    quote(lme4::lmer(B ~ logAge * speciesCode + cover * speciesCode +
                                       (logAge + cover | ecoregionGroup))),
                    NA, NA,
                    paste("Model and formula for estimating biomass (B) from ecoregionGroup",
                          "(currently ecoregionLayer * LandCoverClass), speciesCode,",
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
                          "be an overesimate of the %biomass of hardwoods. E.g., 30% cover of hardwoods",
                          "might translate to 20% biomass of hardwoods. The reason this discount exists is",
                          "because hardwoods in Canada have a much wider canopy than softwoods.")),
    defineParameter("deciduousCoverDiscount", "numeric", 0.8418911, NA, NA,
                    paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                          "Will not be used if P(sim)$fitDeciduousCoverDiscount is TRUE")),
    defineParameter("fitDeciduousCoverDiscount", "logical", FALSE, NA, NA,
                    paste("If TRUE, this will re-estimate deciduousCoverDiscount. This may be unstable and",
                          "is not recommended currently. If FALSE, will use the current default")),
    ###########################################################################################

    defineParameter("ecoregionLayerField", "character", NULL, NA, NA,
                    paste("the name of the field used to distinguish ecoregions, if supplying a polygon.",
                          "Defaults to NULL and tries to use  'ECODISTRIC' where available (for legacy reasons), or the row numbers of",
                          "sim$ecoregionLayer. If this field is not numeric, it will be coerced to numeric")),
    defineParameter("exportModels", "character", "none", NA, NA,
                    paste("Controls whether models used to estimate maximum B/ANPP ('biomassModel') and species establishment",
                          "('coverModel') probabilities are exported for posterior analyses or not. This may be important",
                          "when models fail to converge or hit singularity (but can still be used to make predictions) and",
                          "the user wants to investigate them further. Can be set to 'none' (no models are exported), 'all'",
                          "(both are exported), 'biomassModel' or 'coverModel'.")),
    defineParameter("forestedLCCClasses", "numeric", c(1:15, 20, 32, 34:35), 0, 39,
                    paste("The classes in the rstLCC layer that are 'treed' and will therefore be run in Biomass_core.",
                          "Defaults to forested classes in LCC2005 map.")),
    defineParameter("imputeBadAgeModel", "call",
                    quote(lme4::lmer(age ~ log(totalBiomass) * cover * speciesCode + (log(totalBiomass) | initialEcoregionCode))),
                    NA, NA,
                    paste("Model and formula used for imputing ages that are either missing or do not match well with",
                          "Biomass or Cover. Specifically, if Biomass or Cover is 0, but age is not, then age will be imputed.",
                          "Similarly, if Age is 0 and either Biomass or Cover is not, then age will be imputed")),
    defineParameter("LCCClassesToReplaceNN", "numeric", 34:35, NA, NA,
                    paste("This will replace these classes on the landscape with the closest forest class P(sim)$forestedLCCClasses.",
                          "If the user is using the default 2005 data product for rstLCC, then users may wish to",
                          "include 36 (cities -- if running a historic range of variation project), and 34:35 (burns)",
                          "Since this is about estimating parameters for growth, it doesn't make any sense to have",
                          "unique estimates for transient classes in most cases. If no classes are to be replaced, pass",
                          "'LCCClassesToReplaceNN' = numeric(0) when supplying parameters.")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    "Cover that is equal to or below this number will be omitted from the dataset"),
    defineParameter("minRelativeBFunction", "call", quote(LandR::makeMinRelativeB(pixelCohortData)),
                    NA, NA,
                    paste("A quoted function that makes the table of min. relative B determining",
                          "a stand shade level for each ecoregionGroup. Using the internal object",
                          "`pixelCohortData` is advisable to access/use the list of ecoregionGroups",
                          "per pixel. The function must output a data.frame with 6 columns, named 'ecoregionGroup'",
                          "and 'X1' to 'X5', with one line per ecoregionGroup code ('ecoregionGroup'), and",
                          "the min. relative biomass for each stand shade level X1-5. The default function uses",
                          "values from LANDIS-II available at:",
                          paste0("https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/LandisInputs/BSW/",
                                 "biomass-succession-main-inputs_BSW_Baseline.txt%7E."))),
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
    defineParameter("speciesTableAreas", "character", c("BSW", "BP", "MC"), NA, NA,
                    paste("One or more of the Ecoprovince short forms that are in the `speciesTable` file,",
                          "e.g., BSW, MC etc. Default is good for Alberta and maybe other places.")),
    defineParameter("subsetDataAgeModel", "numeric", 50, NA, NA,
                    "the number of samples to use when subsampling the biomass data model; if TRUE, uses 50"),
    defineParameter("subsetDataBiomassModel", "numeric", NULL, NA, NA,
                    "the number of samples to use when subsampling the biomass data model; if TRUE, uses 50"),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter("useCloudCacheForStats", "logical", TRUE, NA, NA,
                    paste("Some of the statistical models take long (at least 30 minutes, likely longer).",
                          "If this is TRUE, then it will try to get previous cached runs from googledrive.")),
    defineParameter(".plotInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first plot event should occur"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If NA, a hash of studyArea will be used."),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
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
    expectsInput("ecoregionRst", "RasterLayer",
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
                 sourceURL = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"),
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
                              "This is the scale used for all *outputs* for use in the simulation.")),
    expectsInput("rasterToMatchLarge", "RasterLayer",
                 desc = paste("A raster of the studyAreaLarge in the same resolution and projection as rawBiomassMap.",
                              "This is the scale used for all *inputs* for use in the simulation.")),
    expectsInput("rawBiomassMap", "RasterLayer",
                 desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001 (in tonnes/ha). If necessary, biomass values are rescaled to match changes in resolution.",
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
                 desc = "named character vector of hex colour codes corresponding to each species"),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See LandR::sppEquivalencies_CA."),
    expectsInput("sppNameVector", "character",
                 desc = paste("an optional vector of species names to be pulled from sppEquiv. If not provided,",
                              "then species will be taken from the entire P(sim)$sppEquivCol in sppEquiv.",
                              "See LandR::sppEquivalencies_CA.")),
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
                              "Defaults to  an area in Southwestern Alberta, Canada.")),
    expectsInput("studyAreaLarge", "SpatialPolygonsDataFrame",
                 desc = paste("multipolygon (larger area than studyArea) used for parameter estimation,",
                              "with attribute LTHFC describing the fire return interval.",
                              "Defaults to a square shapefile in Southwestern Alberta, Canada."))
  ),
  outputObjects = bindrows(
    createsOutput("biomassMap", "RasterLayer",
                  desc = paste("total biomass raster layer in study area,",
                               "filtered for pixels covered by cohortData. Units in g/m2")),
    createsOutput("cohortData", "data.table",
                  desc = paste("initial community table, created from available biomass (g/m2),",
                               "age and species cover data, as well as eco zonation information")),
    createsOutput("ecoregion", "data.table",
                  desc = "ecoregion look up table"),
    createsOutput("ecoregionMap", "RasterLayer",
                  desc = "ecoregion map that has mapcodes match ecoregion table and speciesEcoregion table"),
    createsOutput("pixelGroupMap", "RasterLayer",
                  desc = "initial community map that has mapcodes match initial community table"),
    createsOutput("pixelFateDT", "data.table",
                  desc = paste("A small table that keeps track of the pixel removals and cause. This may help diagnose issues",
                               "related to understanding the creation of cohortData")),
    createsOutput("minRelativeB", "data.frame",
                  desc = "define the cut points to classify stand shadeness"),
    createsOutput("rawBiomassMap", "RasterLayer",
                  desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                               "Service, National Forest Inventory, kNN-derived total aboveground biomass map (in tonnes/ha)",
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
                  desc = paste("define how the species with different shade tolerance respond to stand shadeness.",
                               "Table values follow LANDIS-II test traits available at: ",
                               paste0("https://raw.githubusercontent.com/LANDIS-II-Foundation/",
                                      "Extensions-Succession/master/biomass-succession-archive/",
                                      "trunk/tests/v6.0-2.0/biomass-succession_test.txt")))
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
  origDTthreads <- data.table::getDTthreads()
  data.table::setDTthreads(min(origDTthreads, 2)) # seems to only improve up to 2 threads
  on.exit(setDTthreads(origDTthreads))

  # # ! ----- EDIT BELOW ----- ! #
  if (is.null(P(sim)$pixelGroupAgeClass))
    params(sim)[[currentModule(sim)]]$pixelGroupAgeClass <- P(sim)$successionTimestep

  cacheTags <- c(currentModule(sim), "init")

  message(blue("Starting to createBiomass_coreInputs in Biomass_borealDataPrep: ", Sys.time()))
  if (is.null(sim$speciesLayers))
    stop(red(paste("'speciesLayers' are missing in Biomass_borealDataPrep init event.\n",
                   "This is likely due to the module producing 'speciesLayers' being scheduled after Biomass_borealDataPrep.\n",
                   "Please check module order.")))

  if (!all(P(sim)$LCCClassesToReplaceNN %in% P(sim)$forestedLCCClasses)) {
    stop("All 'LCCClassesToReplaceNN' should be included in 'forestedLCCClasses'.")
  }

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
                                  sppEquiv = sim$sppEquiv,
                                  areas = P(sim)$speciesTableAreas,
                                  sppEquivCol = P(sim)$sppEquivCol)

  # sim$species <- prepSpeciesTable(speciesTable = sim$speciesTable,
  #                                 # speciesLayers = sim$speciesLayers,
  #                                 sppEquiv = sim$sppEquiv[get(P(sim)$sppEquivCol) %in%
  #                                                           names(sim$speciesLayers)],
  #                                 sppEquivCol = P(sim)$sppEquivCol)

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
      stop("speciesUpdateFunction should be a list of quoted function expressions e.g.:
           list(quote(LandR::speciesTableUpdate(...)), quote(speciesTableUpdateCustom(...)))")
    }
  }

  if (getOption("LandR.verbose") > 0) {
    message("Adjusting species-level traits, part 2, for LandWeb")
    print(sim$species)
  }

  ## check that all species have trait values.
  missingTraits <- setdiff(names(sim$speciesLayers), sim$species$species)
  if (length(missingTraits) == length(names(sim$speciesLayers))) {
    stop("No trait values were found for ", paste(missingTraits, collapse = ", "), ".\n",
         "Please check the species list and traits table")
  } else if (length(missingTraits))
    stop("No trait values were found for ", paste(missingTraits, collapse = ", "), ".\n",
            "Missing traits will result in species removal from simulation.\n
            Please check the species list and traits table")

  ### make table of light shade tolerance  #######################
  ## D. Cyr's version: seems to exacerbate no. of cohorts in our simulations
  ## https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/LandisInputs/BSW/biomass-succession-main-inputs_BSW_Baseline.txt%7E
  ## a prob of 0.5 over 10yrs virtually always results in the successful establishment of a cohort: 1 - (1 - 0.5)^10 = 0.9990234
  # sim$sufficientLight <- data.frame(speciesshadetolerance = 1:5,
  #                                   X0 = 1,
  #                                   X1 = c(0.5, rep(1, 4)),
  #                                   X2 = c(0, 0.5, rep(1, 3)),
  #                                   X3 = c(rep(0, 2), 0.5, rep(1, 2)),
  #                                   X4 = c(rep(0, 3), 0.5, 1),
  #                                   X5 = c(rep(0, 4), 1))

  ## LANDIS-test table (see source in metadata desc.)
  sim$sufficientLight <- data.frame(speciesshadetolerance = 1:5,
                                    X0 = c(rep(1, 4), 0),
                                    X1 = c(0, rep(1, 3), 0),
                                    X2 = c(0, 0, rep(1, 3)),
                                    X3 = c(rep(0, 3), rep(1, 2)),
                                    X4 = c(rep(0, 4), 1),
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
  ## remove pixels with no species data
  # pixelsToRm <- rowSums(!is.na(sim$speciesLayers[])) == 0 # keep
  pixelsToRm <- is.na(sim$speciesLayers[[1]][]) # seems to be OK because seem to be NA on each layer for a given pixel
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
                           tail(pixelFateDT$runningPixelTotal, 1) -
                             NROW(unique(pixelCohortData$pixelIndex)))

  #######################################################
  # Partition totalBiomass into individual species B, via estimating how %cover and %biomass
  #   are related
  #######################################################
  message(blue("Partitioning totalBiomass per pixel into cohort B as:"))
  if (isTRUE(P(sim)$fitDeciduousCoverDiscount)) {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))

    pixelCohortData[, lcc := as.factor(lcc)]

    plot.it <- FALSE
    sam <- subsetDT(pixelCohortData, by = c("speciesCode", "lcc"),
                    doSubset = P(sim)$subsetDataAgeModel,
                    indices = TRUE)
    pi <- unique(pixelCohortData[sam]$pixelIndex)
    sam <- which(pixelCohortData$pixelIndex %in% pi)

    system.time({
      out <- optimize(interval = c(0.1, 1), f = coverOptimFn, bm = P(sim)$coverPctToBiomassPctModel,
                      pixelCohortData = pixelCohortData, subset = sam, maximum = FALSE)
    })
    params(sim)$Biomass_borealDataPrep$deciduousCoverDiscount <- out$minimum
    if (plot.it) {
      cover2BiomassModel <- coverOptimFn(out$minimum, pixelCohortData, P(sim)$subsetDataAgeModel,
                                         P(sim)$coverPctToBiomassPctModel, returnAIC = FALSE)
      sam1 <- sample(NROW(pixelCohortData), 1e5)
      dev()
      par(mfrow = c(1,2))
      plot(predict(cover2BiomassModel$modelBiomass1$mod,
                   newdata = cover2BiomassModel$pixelCohortData[sam1]),
           log(cover2BiomassModel$pixelCohortData$B / 100)[sam1], pch = ".")
      abline(a = 0, b = 1)

      cover2BiomassModel1 <- coverOptimFn(1, pixelCohortData, P(sim)$subsetDataAgeModel,
                                          P(sim)$coverPctToBiomassPctModel,
                                          returnAIC = FALSE)
      dev()
      plot(predict(cover2BiomassModel1$modelBiomass1$mod,
                   newdata = cover2BiomassModel1$pixelCohortData[sam1]),
           log(cover2BiomassModel1$pixelCohortData$B / 100)[sam1], pch = ".")
      abline(a = 0, b = 1)

      pcd <- pixelCohortData
      bb <- pcd[sample(sam)]
      cc <- bb[, cover3 := cover * c(1, out$minimum)[decid + 1]][
        , actualX := cover3 / sum(cover3) / (cover / 100), by = "pixelIndex"]
      setkey(cc, pixelIndex)
      mean(cc[speciesCode == "Popu_Tre"]$actualX)
    }
  } else {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))
    message(blue("using previously estimated deciduousCoverDiscount:",
                 round(P(sim)$deciduousCoverDiscount, 3)))
  }

  pixelCohortData <- partitionBiomass(x = P(sim)$deciduousCoverDiscount, pixelCohortData)
  set(pixelCohortData, NULL, "B", asInteger(pixelCohortData$B/P(sim)$pixelGroupBiomassClass) *
        P(sim)$pixelGroupBiomassClass)
  set(pixelCohortData, NULL, c("decid", "cover2"), NULL)
  set(pixelCohortData, NULL, "cover", asInteger(pixelCohortData$cover))

  #######################################################
  #replace 34 and 35 and 36 values -- burns and cities -- to a neighbour class *that exists*.
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
  if (length(P(sim)$LCCClassesToReplaceNN)) {
    uwc <- P(sim)$LCCClassesToReplaceNN

    message("Replace ", paste(uwc, collapse = ", "),
            " values -- ", "burns"[any(uwc %in% 34:35)], " and cities"[any(uwc %in% 36)],
            " -- to a neighbour class *that exists*")

    ## version 1: from before March 2019 - Ceres noticed it created issues with fitting modelCover
    ## March 2020: seems to be the preferred behaviour?
    ## June 2020: this leads to ignoring pixels with classes to be converted that have cover > 0
    # availableCombinations <- unique(pixelCohortData[eval(rmZeroBiomassQuote),
    # .(speciesCode, initialEcoregionCode, pixelIndex)])
    ## version 2: Ceres's fix from March 2019 to solve issues with modelCover fitting (?)
    ## June 2020: Ceres re-activated this so that pixels with B == 0 and cover > 0 could be converted if need be
    availableCombinations <- unique(pixelCohortData[, .(speciesCode, initialEcoregionCode, pixelIndex)])
    ## version 3: Feb 2020 Eliot's fix that is WRONG - this behaviour is being achieved in convertUnwantedLCC and creates empty tables if done here
    # availableCombinations <- unique(pixelCohortData[!(lcc %in% uwc),
    #                                                 .(speciesCode, initialEcoregionCode, pixelIndex)])

    newLCCClasses <- Cache(convertUnwantedLCC,
                           classesToReplace = P(sim)$LCCClassesToReplaceNN,
                           rstLCC = rstLCCAdj,
                           availableERC_by_Sp = availableCombinations,
                           userTags = c(cacheTags, "newLCCClasses", "stable"),
                           omitArgs = c("userTags"))
  } else {
    newLCCClasses <- data.table(pixelIndex = numeric(), ecoregionGroup = numeric())
  }

  ## split pixelCohortData into 2 parts -- one with the former 34:36 pixels, one without
  #    The one without 34:36 can be used for statistical estimation, but not the one with
  cohortData34to36 <- pixelCohortData[pixelIndex %in% newLCCClasses$pixelIndex]
  cohortData34to36 <- merge(newLCCClasses, cohortData34to36, all.x = TRUE,
                            all.y = FALSE, by = "pixelIndex")
  cohortDataNo34to36 <- pixelCohortData[!pixelIndex %in% newLCCClasses$pixelIndex]
  if (!length(P(sim)$LCCClassesToReplaceNN)) {
    if (!identical(cohortDataNo34to36, pixelCohortData))
      stop("No LCC classes were listed for replacement, but some pixels may have been lost")
  }
  setnames(cohortDataNo34to36, "initialEcoregionCode", "ecoregionGroup")
  rmZeroBiomassQuote <- quote(totalBiomass > 0)
  cohortDataNo34to36Biomass <- cohortDataNo34to36[eval(rmZeroBiomassQuote),
                                                  .(B, logAge, speciesCode, ecoregionGroup, lcc, cover)]
  cohortDataNo34to36Biomass <- unique(cohortDataNo34to36Biomass)

  ## make sure ecoregionGroups match
  ## remember to match rmZeroBiomassQuote the rule used to filter `availableCombinations` (NULL if none)
  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert1(cohortData34to36, pixelCohortData, rmZeroBiomassQuote = NULL,
            classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(cohortDataNo34to36, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }

  ##############################################################
  # Statistical estimation of establishprob, maxB and maxANPP
  ##############################################################
  cohortDataShort <- cohortDataNo34to36[, list(coverPres = sum(cover > 0)),
                                        by = c("ecoregionGroup", "speciesCode")]
  ## find coverNum for each known class
  ## TODO: Ceres: I feel we should be using the converted classes here...
  ## otherwise  pixelTable is reintroducing the converted classes in cohortDataShortNoCover which is confusing
  ## even if they end up being removed later by `makeSpeciesEcoregion`
  ## because they end up with NAs that are converted to 0s
  # aa <- table(pixelTable$initialEcoregionCode)

  ## add new ecoregions to pixelTable, before calc. table
  tempDT <- rbind(cohortData34to36[, .(pixelIndex, ecoregionGroup)],
                  cohortDataNo34to36[, .(pixelIndex, ecoregionGroup)])
  pixelTable <- tempDT[pixelTable, on = .(pixelIndex)]
  aa <- table(as.character(pixelTable$ecoregionGroup))   ## as.character avoids counting levels that don't exist anymore

  dt1 <- data.table(ecoregionGroup = factor(names(aa)), coverNum = as.integer(unname(aa)))
  allCombos <- expand.grid(ecoregionGroup = dt1$ecoregionGroup, speciesCode = unique(cohortDataShort$speciesCode))
  setDT(allCombos)
  dt1 <- dt1[allCombos, on = "ecoregionGroup", nomatch = 0]
  cohortDataShortNoCover <- cohortDataShort[dt1, on = c("ecoregionGroup", "speciesCode"), nomatch = NA]

  #cohortDataShortNoCover <- cohortDataShort[coverPres == 0] #
  cohortDataShort <- cohortDataShortNoCover[coverPres > 0] # remove places where there is 0 cover
  cohortDataShortNoCover <- cohortDataShortNoCover[is.na(coverPres)][, coverPres := 0]
  # will be added back as establishprob = 0

  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert2(cohortDataShort, classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(cohortDataShortNoCover, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }

  message(blue("Estimating Species Establishment Probability using P(sim)$coverModel, which is"))
  message(magenta(paste0(format(P(sim)$coverModel, appendLF = FALSE), collapse = "")))

  # for backwards compatibility -- change from parameter to object
  if (is.null(sim$cloudFolderID))
    if (!is.null(P(sim)$cloudFolderID))
      sim$cloudFolderID <- P(sim)$cloudFolderID

  useCloud <- if (!is.null(sim$cloudFolderID)) {
    (isTRUE(getOption("reproducible.useCache", FALSE)) && P(sim)$useCloudCacheForStats)
  } else {
    FALSE
  }

  # Remove all cases where there is 100% presence in an ecoregionGroup -- causes failures in binomial models
  cdsWh <- cohortDataShort$coverPres == cohortDataShort$coverNum
  cds <- Copy(cohortDataShort)
  cds <- cds[!cdsWh]

  modelCover <- Cache(
    statsModel,
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
    omitArgs = c("showSimilar", "useCache", ".specialData", "useCloud", "cloudFolderID")
  )
  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelCover$rsq, 4))), function(x) message(blue(x)))

  ## export model before overriding happens
  if (any(P(sim)$exportModels %in% c("all", "coverModel")))
    sim$modelCover <- modelCover

  if (isTRUE(any(cdsWh))) {
    cds[, pred := fitted(modelCover$mod, response = "response")]
    cohortDataShort <- cds[, -c("coverPres", "coverNum")][cohortDataShort,
                                                          on = c('ecoregionGroup', 'speciesCode'), nomatch = NA]
    cohortDataShort[is.na(pred), pred := 1]
    modelCover <- cohortDataShort$pred
  }

  ## For biomass
  ### Subsample cases where there are more than 50 points in an ecoregionGroup * speciesCode
  totalBiomass <- sum(cohortDataNo34to36Biomass$B, na.rm = TRUE)
  cohortDataNo34to36Biomass <- subsetDT(cohortDataNo34to36Biomass,
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
    uniqueEcoregionGroups = .sortDotsUnderscoreFirst(as.character(unique(cohortDataNo34to36Biomass$ecoregionGroup))),
    sumResponse = totalBiomass,
    .specialData = cohortDataNo34to36Biomass,
    useCloud = useCloud,
    # useCache = "overwrite",
    cloudFolderID = sim$cloudFolderID,
    showSimilar = getOption("reproducible.showSimilar", FALSE),
    userTags = c(cacheTags, "modelBiomass", paste0("subsetSize:", P(sim)$subsetDataBiomassModel)),
    omitArgs = c("showSimilar", ".specialData", "useCloud", "cloudFolderID", "useCache")
  )

  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelBiomass$rsq, 4))), function(x) message(blue(x)))

  if (any(P(sim)$exportModels %in% c("all", "biomassModel")))
    sim$modelBiomass <- modelBiomass

  ########################################################################
  # create speciesEcoregion -- a single line for each combination of ecoregionGroup & speciesCode
  #   doesn't include combinations with B = 0 because those places can't have the species/ecoregion combo
  ########################################################################
  ## cohortDataNo34to36Biomass ends up determining which ecoregion combinations end up in
  ## species ecoregion, thus removing converted/masked classes present cohortDataShortNoCover
  message(blue("Create speciesEcoregion using modelCover and modelBiomass to estimate species traits"))
  speciesEcoregion <- makeSpeciesEcoregion(cohortDataBiomass = cohortDataNo34to36Biomass,
                                           cohortDataShort = cohortDataShort,
                                           cohortDataShortNoCover = cohortDataShortNoCover,
                                           species = sim$species,
                                           modelCover = modelCover,
                                           modelBiomass = modelBiomass,
                                           successionTimestep = P(sim)$successionTimestep,
                                           currentYear = time(sim))
  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert2(speciesEcoregion, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }

  ## check that all species have maxB/maxANPP
  assertSppMaxBMaxANPP(speciesEcoregion)

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
    newDev <- if (!is.null(dev.list())) max(dev.list()) + 1 else 1
    quickPlot::dev(newDev, width = 18, height = 10)
    Plot(maxB, legendRange = c(0, max(maxValue(maxB), na.rm = TRUE)))
    quickPlot::dev(curDev)
  }

  if (ncell(sim$rasterToMatchLarge) > 3e6) .gc()

  ########################################################################
  # Create initial communities, i.e., pixelGroups
  ########################################################################
  # Rejoin back the pixels that were 34:36
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

  if (sum(is.na(getValues(sim$rasterToMatch))) != sum(is.na(getValues(sim$rasterToMatchLarge)))) {
    message(blue("Subsetting to studyArea"))
    rasterToMatchLarge <- sim$rasterToMatchLarge
    rasterToMatchLarge <- setValues(rasterToMatchLarge, seq(ncell(rasterToMatchLarge)))
    rasterToMatchLarge <- Cache(postProcess,
                                x = rasterToMatchLarge,
                                rasterToMatch = sim$rasterToMatch,
                                maskWithRTM = TRUE,
                                filename2 = NULL,
                                #useCache = "overwrite",
                                userTags = c(cacheTags, "rasterToMatchLarge"),
                                omitArgs = c("userTags"))

    assertthat::assert_that(sum(is.na(getValues(rasterToMatchLarge))) < ncell(rasterToMatchLarge)) ## i.e., not all NA

    if (!compareRaster(rasterToMatchLarge, sim$rasterToMatch, orig = TRUE, stopiffalse = FALSE))
      stop("Downsizing to rasterToMatch after estimating parameters didn't work.",
           "Please debug Biomass_borealDataPrep::createBiomass_coreInputs().")

    ## subset pixels that are in studyArea/rasterToMatch only
    pixToKeep <- na.omit(getValues(rasterToMatchLarge)) # these are the old indices of RTML
    pixelCohortData <- pixelCohortData[pixelIndex %in% pixToKeep]

    ## re-do pixelIndex (it now needs to match rasterToMatch)
    newPixelIndexDT <- data.table(pixelIndex = getValues(rasterToMatchLarge),
                                  newPixelIndex = as.integer(1:ncell(rasterToMatchLarge))) %>%
      na.omit(.)

    pixelCohortData <- newPixelIndexDT[pixelCohortData, on = "pixelIndex"]
    pixelCohortData[, pixelIndex := NULL]
    setnames(pixelCohortData, old = "newPixelIndex", new = "pixelIndex")
    rm(pixToKeep, rasterToMatchLarge)

    assertthat::assert_that(NROW(pixelCohortData) > 0)

    if (ncell(sim$rasterToMatch) > 3e6) .gc()
  }
  ## subset ecoregionFiles$ecoregionMap to smaller area.
  ecoregionFiles$ecoregionMap <- Cache(postProcess,
                                       x = ecoregionFiles$ecoregionMap,
                                       rasterToMatch = sim$rasterToMatch,
                                       maskWithRTM = TRUE,
                                       filename2 = NULL,
                                       userTags = c(cacheTags, "ecoregionMap"),
                                       omitArgs = c("userTags"))

  maxAgeHighQualityData <- -1

  # If this module used a fire database to extract better young ages, then we
  #   can use those high quality younger ages to help with our biomass estimates
  if (length(extractURL("fireURL"))) {
    # fireURL <- "https://cwfis.cfs.nrcan.gc.ca/downloads/nbac/nbac_1986_to_2019_20200921.zip"
    # This was using the nbac filename to figure out what the earliest year in the
    #   fire dataset was. Since that is not actually used here, it doesn't really
    #   matter what the fire dataset was. Basically, this section is updating
    #   young ages that are way outside of their biomass. Can set this to 1986 to just
    #   give a cutoff
    firstFireYear <- 1986 # as.numeric(gsub("^.+nbac_(.*)_to.*$", "\\1", fireURL))
    maxAgeHighQualityData <- start(sim) - firstFireYear
    ## if maxAgeHighQualityData is lower than 0, it means it's prior to the first fire Year
    ## or not following calendar year
    if (!is.na(maxAgeHighQualityData) & maxAgeHighQualityData >= 0) {
      youngRows <- pixelCohortData$age <= maxAgeHighQualityData
      young <- pixelCohortData[youngRows == TRUE]

      # whYoungBEqZero <- which(young$B == 0)
      whYoungZeroToMaxHighQuality <- which(young$age > 0)
      if (length(whYoungZeroToMaxHighQuality) > 0) {
        youngWAgeEqZero <- young[-whYoungZeroToMaxHighQuality]
        youngNoAgeEqZero <- young[whYoungZeroToMaxHighQuality]

        young <- Cache(updateYoungBiomasses,
                       young = youngNoAgeEqZero,
                       biomassModel = modelBiomass$mod,
                       userTags = c(cacheTags, "updateYoungBiomasses"),
                       omitArgs = c("userTags"))
        set(young, NULL, setdiff(colnames(young), colnames(pixelCohortData)), NULL)

        young <- rbindlist(list(young, youngWAgeEqZero), use.names = TRUE)

      }
      pixelCohortData <- rbindlist(list(pixelCohortData[youngRows == FALSE], young), use.names = TRUE)
    } else {
      ## return maxAgeHighQualityData to -1
      maxAgeHighQualityData <- -1
    }
  }

  # Fill in any remaining B values that are still NA -- the previous chunk filled in B for young cohorts only
  if (anyNA(pixelCohortData$B)) {
    theNAsBiomass <- is.na(pixelCohortData$B)
    message(blue(" -- ", sum(theNAsBiomass),"cohort(s) has NA for Biomass: being replaced with model-derived estimates"))
    set(pixelCohortData, which(theNAsBiomass), "B",
        asInteger(predict(modelBiomass$mod, newdata = pixelCohortData[theNAsBiomass])))
  }

  ## make cohortDataFiles: pixelCohortData (rm unnecessary cols, subset pixels with B>0,
  ## generate pixelGroups, add ecoregionGroup and totalBiomass) and cohortData
  cohortDataFiles <- makeCohortDataFiles(pixelCohortData, columnsForPixelGroups, speciesEcoregion,
                                         pixelGroupBiomassClass = P(sim)$pixelGroupBiomassClass,
                                         pixelGroupAgeClass = P(sim)$pixelGroupAgeClass,
                                         minAgeForGrouping = maxAgeHighQualityData,
                                         pixelFateDT = pixelFateDT)

  sim$cohortData <- cohortDataFiles$cohortData
  pixelCohortData <- cohortDataFiles$pixelCohortData
  pixelFateDT <- cohortDataFiles$pixelFateDT

  rm(cohortDataFiles)
  assertthat::assert_that(NROW(pixelCohortData) > 0)
  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert2(pixelCohortData, classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(sim$cohortData, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }
  ## make a table of available active and inactive (no biomass) ecoregions
  sim$ecoregion <- makeEcoregionDT(pixelCohortData, speciesEcoregion)

  ## make biomassMap, ecoregionMap, minRelativeB, pixelGroupMap (at the scale of rasterToMatch)
  sim$biomassMap <- makeBiomassMap(pixelCohortData, sim$rasterToMatch)
  sim$ecoregionMap <- makeEcoregionMap(ecoregionFiles, pixelCohortData)

  sim$pixelGroupMap <- makePixelGroupMap(pixelCohortData, sim$rasterToMatch)

  if (is(P(sim)$minRelativeBFunction, "call")) {
    sim$minRelativeB <- eval(P(sim)$minRelativeBFunction)
  } else {
    stop("minRelativeBFunction should be a quoted function expression, using `pixelCohortData`, e.g.:\n",
         "    quote(LandR::makeMinRelativeB(pixelCohortData))")
  }
  ## make sure speciesLayers match RTM (since that's what is used downstream in simulations)
  message(blue("Writing sim$speciesLayers to disk as they are likely no longer needed in RAM"))

  sim$speciesLayers <- Cache(postProcess, sim$speciesLayers,
                             rasterToMatch = sim$rasterToMatch,
                             maskWithRTM = TRUE,
                             filename2 = .suffix(file.path(outputPath(sim), 'speciesLayers.grd'),
                                                 paste0("_", P(sim)$.studyAreaName)),
                             overwrite = TRUE,
                             userTags = c(cacheTags, "speciesLayersRTM"),
                             omitArgs = c("userTags"))

  ## double check these rasters all match RTM
  compareRaster(sim$biomassMap, sim$ecoregionMap, sim$pixelGroupMap, sim$rasterToMatch, sim$speciesLayers)

  ## rm ecoregions that may not be present in rasterToMatch
  ## make ecoregionGroup a factor and export speciesEcoregion to sim
  onMatch <- c("ecoregionGroup", "speciesCode")
  toRm <- speciesEcoregion[!sim$cohortData, on = onMatch]
  speciesEcoregion <- speciesEcoregion[!toRm, on = onMatch]
  sim$speciesEcoregion <- speciesEcoregion
  sim$speciesEcoregion$ecoregionGroup <- factor(as.character(sim$speciesEcoregion$ecoregionGroup))

  ## do assertions
  message(blue("Create pixelGroups based on: ", paste(columnsForPixelGroups, collapse = ", "),
               "\n  Resulted in", magenta(length(unique(sim$cohortData$pixelGroup))),
               "unique pixelGroup values"))
  assertSpeciesEcoregionCohortDataMatch(sim$cohortData, sim$speciesEcoregion, doAssertion = TRUE)

  # LandR::assertERGs(sim$ecoregionMap, cohortData = sim$cohortData,
  #                  speciesEcoregion = sim$speciesEcoregion,
  #                  minRelativeB = sim$minRelativeB, doAssertion = TRUE)

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

  if (!compareCRS(sim$studyArea, sim$studyAreaLarge)) {
    warning("studyArea and studyAreaLarge have different projections.\n
            studyAreaLarge will be projected to match crs(studyArea)")
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$studyArea))
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyAreaLarge)
    message("The .studyAreaName is not supplied; derived name from sim$studyAreaLarge: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
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
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), { ## TODO: re-enable verify
      #necessary for KNN
      sim$rawBiomassMap <- Cache(prepInputs,
                                 url = extractURL("rawBiomassMap"),
                                 destinationPath = dPath,
                                 studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                                 rasterToMatch = if (!needRTM) sim$rasterToMatchLarge else NULL,
                                 maskWithRTM = if (!needRTM) TRUE else FALSE,
                                 useSAcrs = FALSE,     ## never use SA CRS
                                 method = "bilinear",
                                 datatype = "INT2U",
                                 filename2 = .suffix("rawBiomassMap.tif", paste0("_", P(sim)$.studyAreaName)),
                                 overwrite = TRUE,
                                 userTags = c(cacheTags, "rawBiomassMap"),
                                 omitArgs = c("destinationPath", "targetFile", "userTags", "stable"))
    })
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

    sim$rasterToMatchLarge <- Cache(
      writeOutputs,
      sim$rasterToMatchLarge,
      filename2 = .suffix(file.path(dPath, "rasterToMatchLarge.tif"),
                          paste0("_", P(sim)$.studyAreaName)),
      datatype = "INT2U",
      overwrite = TRUE,
      userTags = c(cacheTags, "rasterToMatchLarge"),
      omitArgs = c("userTags")
    )

    sim$rasterToMatch <- Cache(postProcess,
                               x = sim$rawBiomassMap,
                               studyArea = sim$studyArea,
                               # rasterToMatch = sim$rasterToMatchLarge,   ## Ceres: this messes up the extent. if we are doing this it means BOTH RTMs come from biomassMap, so no need for RTMLarge here.
                               useSAcrs = FALSE,
                               # maskWithRTM = FALSE,   ## mask with SA
                               method = "bilinear",
                               datatype = "INT2U",
                               filename2 = .suffix(file.path(dPath, "rasterToMatch.tif"),
                                                   paste0("_", P(sim)$.studyAreaName)),
                               overwrite = TRUE,
                               useCache = "overwrite",
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
  if (!compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- spTransform(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  if (!compareCRS(sim$studyAreaLarge, sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- spTransform(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  ## Land cover raster ------------------------------------------------
  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- prepInputsLCC(
      destinationPath = dPath,
      studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
      rasterToMatch = sim$rasterToMatchLarge,
      filename2 = .suffix("rstLCC.tif", paste0("_", P(sim)$.studyAreaName)),
      overwrite = TRUE,
      userTags = c("rstLCC", currentModule(sim), P(sim)$.studyAreaName))
  }

  if (!compareRaster(sim$rstLCC, sim$rasterToMatchLarge)) {
    sim$rstLCC <- projectRaster(sim$rstLCC, to = sim$rasterToMatchLarge)
  }

  ## Ecodistrict ------------------------------------------------
  if (!suppliedElsewhere("ecoregionLayer", sim)) {
    sim$ecoregionLayer <- Cache(prepInputs,
                                targetFile = "ecodistricts.shp",
                                archive = asPath("ecodistrict_shp.zip"),
                                url = extractURL("ecoregionLayer", sim),
                                alsoExtract = "similar",
                                destinationPath = dPath,
                                filename2 = NULL,
                                studyArea = sim$studyAreaLarge,   ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
                                overwrite = TRUE,
                                useSAcrs = TRUE, # this is required to make ecoZone be in CRS of studyArea
                                fun = "raster::shapefile",
                                userTags = c("prepInputsEcoDistrict_SA", currentModule(sim), cacheTags))
  }

  ## Stand age map ------------------------------------------------
  if (!suppliedElsewhere("standAgeMap", sim)) {
    httr::with_config(config = httr::config(ssl_verifypeer = 0L), {
      sim$standAgeMap <- Cache(LandR::prepInputsStandAgeMap,
                               destinationPath = dPath,
                               ageURL = extractURL("standAgeMap"),
                               studyArea = raster::aggregate(sim$studyAreaLarge),
                               rasterToMatch = sim$rasterToMatchLarge,
                               filename2 = .suffix("standAgeMap.tif", paste0("_", P(sim)$.studyAreaName)),
                               overwrite = TRUE,
                               fireURL = extractURL("fireURL"),
                               fireField = "YEAR",
                               startTime = start(sim),
                               userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
                               omitArgs = c("destinationPath", "targetFile", "overwrite",
                                            "alsoExtract", "userTags"))
    })
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

    ## make sure empty pixels inside study area have 0 cover, instead of NAs.
    ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
    tempRas <- sim$rasterToMatchLarge
    tempRas[!is.na(tempRas[])] <- 0
    sim$speciesLayers <- cover(sim$speciesLayers, tempRas)
    rm(tempRas)
  }

  # 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath = dPath, cacheTags = c(cacheTags, "speciesTable"))
  }

  if (!suppliedElsewhere("columnsForPixelGroups", sim)) {
    sim$columnsForPixelGroups <- LandR::columnsForPixelGroups
  }

  return(invisible(sim))
}
