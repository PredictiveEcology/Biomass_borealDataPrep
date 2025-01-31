defineModule(sim, list(
  name = "Biomass_borealDataPrep",
  description = paste("A data preparation module for parameterizing `Biomass_core` from open data sources,",
                      "within the Boreal forest of Canada."),
  keywords = c("LandWeb", "Biomass_core"),
  authors = c(
    person("Yong", "Luo", email = "Yong.Luo@gov.bc.ca", role = c("aut")),
    person(c("Eliot", "J", "B"), "McIntire", email = "eliot.mcintire@nrcan-rncan.gc.ca", role = c("aut", "cre")),
    person(c("Ceres"), "Barros", email = "ceres.barros@ubc.ca", role = c("aut")),
    person(c("Alex", "M."), "Chubaty", email = "achubaty@for-cast.ca", role = c("aut"))
  ),
  childModules = character(0),
  version = list(Biomass_borealDataPrep = "1.5.7.9004"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_borealDataPrep.Rmd"),
  loadOrder = list(after = c("Biomass_speciesData"),
                   before = c("Biomass_core")),
  reqdPkgs = list("assertthat", "crayon", "data.table", "dplyr", "fasterize",  "ggplot2",
                  "merTools", "plyr", "rasterVis", "sf", "terra",
                  "reproducible (>= 2.1.0)",
                  "SpaDES.core (>= 2.1.0)", "SpaDES.tools (>= 2.0.0)",
                  "PredictiveEcology/LandR@development (>= 1.1.5.9025)",
                  "PredictiveEcology/SpaDES.project@development (>= 0.0.8.9026)", ## TODO: update this once merged
                  "PredictiveEcology/pemisc@development"),
  parameters = rbind(
    ## maxB, maxANPP, SEP estimation section ------------------------------------------------
    defineParameter("biomassModel", "call",
                    quote(lme4::lmer(B ~ logAge * speciesCode + cover * speciesCode +
                                       (logAge + cover | ecoregionGroup))),
                    NA, NA,
                    paste("Model and formula for estimating biomass (B) from `ecoregionGroup`",
                          "(currently `ecoregionLayer` * `LandCoverClass`), `speciesCode`,",
                          "`logAge` (gives a downward curving relationship), and `cover`.",
                          "Defaults to a LMEM, which can be slow if dealing with very large datasets",
                          "(e.g., 36,000 points takes 20 minutes).",
                          "For faster fitting try `P(sim)$subsetDataBiomassModel == TRUE`, or",
                          "`quote(RcppArmadillo::fastLm(formula = B ~ logAge * speciesCode * ecoregionGroup",
                          "+ cover * speciesCode * ecoregionGroup))`.",
                          "A custom model call can also be provided, as long as the 'data' argument",
                          "is NOT included.")),
    defineParameter("coverModel", "call",
                    quote(glm(cbind(coverPres, coverNum - coverPres) ~ speciesCode * ecoregionGroup,
                              family = binomial)),
                    NA, NA,
                    paste("Model and formula used for estimating cover from `ecoregionGroup` and `speciesCode`",
                          "and potentially others. Defaults to a GLMEM if there are > 1 grouping levels.",
                          "A custom model call can also be provided, as long as the 'data' argument is NOT included")),
    defineParameter("fixModelBiomass", "logical", FALSE, NA, NA,
                    paste("should `biomassModel` be fixed in the case of non-convergence?",
                          "Only scaling of variables and attempting to fit with a new optimizer (bobyqa, see `?lme4`)",
                          "are implemented at this time.")),
    defineParameter("subsetDataAttempts", "integer", 3L, 1L, 10L,
                    paste("How many times should `biomassModel` be attempted to fit with a new data subset in case of",
                          "non-convergence? Each time, the data is resampled (if `subsetDataBiomassModel = TRUE`)",
                          "and the model re-fit with the original data, scaled variables and/or a different optimizer",
                          "if `fixModelBiomass = TRUE`. Model refiting with original data, rescaled variables and/or a new",
                          "optimizer occurs up to three times for each data subset, regardless of this parameter's value.")),
    defineParameter("subsetDataBiomassModel", "integer", 50L, NA_integer_, NA_integer_,
                    paste("the number of samples to use when subsampling the biomass data model (`biomassModel`);",
                          "Can be `TRUE`/`FALSE`/`NULL` or numeric; if `TRUE`, uses 50, the default.",
                          "If `FALSE`/`NULL` no subsetting is done.")),
    ## deciduous cover to biomass cover section ------------------------------------------------
    defineParameter("coverPctToBiomassPctModel", "call",
                    quote(glm(I(log(B/100)) ~ logAge * I(log(totalBiomass/100)) * speciesCode * lcc)),
                    NA, NA,
                    paste("Model to estimate the relationship between % cover and % biomass, referred to as",
                          "`P(sim)$fitDeciduousCoverDiscount` It is a number between 0 and 1 that translates % cover,",
                          "as provided in several databases, to % biomass. It is assumed that all hardwoods",
                          "are equivalent and all softwoods are equivalent and that % cover of hardwoods will",
                          "be an overesimate of the % biomass of hardwoods. E.g., 30% cover of hardwoods",
                          "might translate to 20% biomass of hardwoods. The reason this discount exists is",
                          "because hardwoods in Canada have a much wider canopy than softwoods.")),
    defineParameter("deciduousCoverDiscount", "numeric", 0.8418911, NA, NA,
                    paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                          "Will not be used if `P(sim)$fitDeciduousCoverDiscount == TRUE`")),
    defineParameter("fitDeciduousCoverDiscount", "logical", FALSE, NA, NA,
                    paste("If TRUE, this will re-estimate `P(sim)$fitDeciduousCoverDiscount` This may be unstable and",
                          "is not recommended currently. If `FALSE`, will use the current default")),
    ## -------------------------------------------------------------------------------------------
    defineParameter("dataYear", "numeric", 2001, NA, NA,
                    paste("Used to override the default 'sourceURL' of KNN datasets (species cover, stand biomass",
                          "and stand age), which point to 2001 data, to fetch KNN data for another year. Currently,",
                          "the only other possible year is 2011. Will also select NTEMS landcover from appropriate year.")),
    defineParameter("ecoregionLayerField", "character", NULL, NA, NA,
                    paste("the name of the field used to distinguish ecoregions, if supplying a polygon.",
                          "Defaults to `NULL` and tries to use  'ECODISTRIC' where available (for legacy reasons), or the row numbers of",
                          "`sim$ecoregionLayer`. If this field is not numeric, it will be coerced to numeric.")),
    defineParameter("exportModels", "character", "none", NA, NA,
                    paste("Controls whether models used to estimate maximum B/ANPP (`biomassModel`) and species establishment",
                          "(`coverModel`) probabilities are exported for posterior analyses or not. This may be important",
                          "when models fail to converge or hit singularity (but can still be used to make predictions) and",
                          "the user wants to investigate them further. Can be set to 'none' (no models are exported), 'all'",
                          "(both are exported), 'biomassModel' or 'coverModel'. BEWARE: because this is intended for posterior",
                          "model inspection, the models will be exported with data, which may mean very large simList(s)!")),
    defineParameter("forestedLCCClasses", "numeric", c(81, 210, 220, 230, 240), 0, NA,
                    paste("The classes in the `rstLCC` layer that are 'treed' and will therefore be run in Biomass_core.",
                          "Defaults to forested classes in NTEMS map (210 = conif, 220 deciduous, 230 mixed) plus",
                          "LandR-generated 240 class, which is recently disturbed forest.")),
    defineParameter("imputeBadAgeModel", "call",
                    quote(lme4::lmer(age ~ log(totalBiomass) * cover * speciesCode + (log(totalBiomass) | initialEcoregionCode))),
                    NA, NA,
                    paste("Model and formula used for imputing ages that are either missing or do not match well with",
                          "biomass or cover. Specifically, if biomass or cover is 0, but age is not, or if age is missing (`NA`),",
                          "then age will be imputed. Note that this is independent from replacing ages inside fire perimeters",
                          "(see `P(sim)$overrideAgeInFires`)")),
    defineParameter("LCCClassesToReplaceNN", "numeric", 240, NA, NA,
                    paste("This will replace these classes on the landscape with the closest forest class `P(sim)$forestedLCCClasses`.",
                          "If the user is using the LCC 2005 land-cover data product for `rstLCC`, then they may wish to",
                          "include 36 (cities -- if running a historic range of variation project), and 34:35 (burns)",
                          "Since this is about estimating parameters for growth, it doesn't make any sense to have",
                          "unique estimates for transient classes in most cases. If no classes are to be replaced, pass",
                          "`'LCCClassesToReplaceNN' = numeric(0)` when supplying parameters.")),
    defineParameter("minCoverThreshold", "numeric", 5, 0, 100,
                    "Pixels with total cover that is equal to or below this number will be omitted from the dataset"),
    defineParameter("minRelativeBFunction", "call", quote(LandR::makeMinRelativeB(pixelCohortData)),
                    NA, NA,
                    paste(
                      "A quoted function that makes the table of min. relative B determining",
                      "a stand shade level for each `ecoregionGroup`. Using the internal object",
                      "`pixelCohortData` is advisable to access/use the list of `ecoregionGroup`s per pixel.",
                      "The function must output a `data.frame` with 6 columns, named `ecoregionGroup`",
                      "and 'X1' to 'X5', with one line per `ecoregionGroup` code, and",
                      "the min. relative biomass for each stand shade level X1-5.",
                      "The default function uses values from LANDIS-II available at:",
                      paste0("https://github.com/dcyr/LANDIS-II_IA_generalUseFiles/blob/master/",
                             "LandisInputs/BSW/biomass-succession-main-inputs_BSW_Baseline.txt"),
                      "and applies them to all ecolocations (`ecoregionGroup` codes)"
                    )),
    defineParameter("omitNonTreedPixels", "logical", TRUE, FALSE, TRUE,
                    "Should this module use only treed pixels, as identified by `P(sim)$forestedLCCClasses`?"),
    defineParameter("overrideAgeInFires", "logical", TRUE, NA, NA,
                    paste("should stand age values inside fire perimeters be replaced with number of years since last fire?")),
    defineParameter("overrideBiomassInFires", "logical", TRUE, NA, NA,
                    paste("should B values be re-estimated using *Biomass_core* for pixels within the fire perimeters",
                          "for which age was replaced with time since last fire? Ignored if `P(sim)$overrideAgeInFires = FALSE`. ",
                          "See `firePerimeters` input object and `P(sim)$overrideAgeInFires` for further detail.")),
    defineParameter("pixelGroupAgeClass", "numeric", params(sim)$Biomass_borealDataPrep$successionTimestep, NA, NA,
                    paste("When assigning `pixelGroup` membership, this defines the resolution of ages that will be considered",
                          "'the same pixelGroup', e.g., if it is 10, then 6 and 14 will be the same")),
    defineParameter("pixelGroupBiomassClass", "numeric", 100, NA, NA,
                    paste("When assigning pixelGroup membership, this defines the resolution of biomass that will be considered",
                          "'the same pixelGroup', e.g., if it is 100, then 5160 and 5240 will be the same")),
    defineParameter("rmImputedPix", "logical", FALSE, NA, NA,
                    "Should `sim$imputedPixID` be removed from the simulation?"),
    defineParameter("speciesUpdateFunction", "list",
                    list(quote(LandR::speciesTableUpdate(sim$species, sim$speciesTable, sim$sppEquiv, P(sim)$sppEquivCol))),
                    NA, NA,
                    paste("Unnamed list of (one or more) quoted functions that updates species table to customize values.",
                          "By default, `LandR::speciesTableUpdate` is used to change longevity and shade tolerance values,",
                          "using values appropriate to Boreal Shield West (BSW), Boreal Plains (BP) and Montane Cordillera (MC)",
                          "ecoprovinces (see `?LandR::speciesTableUpdate` for details). Set to `NULL` if default trait values from",
                          "`speciesTable` are to be kept instead. The user can supply other or additional functions to change",
                          "trait values (see `LandR::updateSpeciesTable`)")),
    defineParameter("sppEquivCol", "character", "Boreal", NA, NA,
                    "The column in `sim$speciesEquivalency` data.table to use as a naming convention."),
    defineParameter("speciesTableAreas", "character", c("BSW", "BP", "MC"), NA, NA,
                    paste("One or more of the Ecoprovince short forms that are in the `speciesTable` file,",
                          "e.g., BSW, MC etc. Default is good for Alberta and other places in the western Canadian boreal forests.")),
    defineParameter("subsetDataAgeModel", "numeric", 50, NA, NA,
                    paste("the number of samples to use when subsampling the age data model and when fitting `coverPctToBiomassPctModel`;",
                          "Can be `TRUE`/`FALSE`/`NULL` or numeric; if `TRUE`, uses 50, the default.",
                          "If `FALSE`/`NULL` no subsetting is done.")),
    defineParameter("successionTimestep", "numeric", 10, NA, NA, "defines the simulation time step, default is 10 years"),
    defineParameter("useCloudCacheForStats", "logical", TRUE, NA, NA,
                    paste("Some of the statistical models take long (at least 30 minutes, likely longer).",
                          "If this is `TRUE`, then it will try to get previous cached runs from googledrive.")),
    defineParameter("vegLeadingProportion", "numeric", 0.8, 0, 1,
                    desc = "a number that defines whether a species is leading for a given pixel"),
    defineParameter(".plotInitialTime", "numeric", start(sim), NA, NA,
                    "This is here for backwards compatibility. Please use `.plots`"),
    defineParameter(".plots", "character", NA, NA, NA,
                    "This describes the type of 'plotting' to do. See `?Plots` for possible types. To omit, set to NA"),
    defineParameter(".plotInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between plot events"),
    defineParameter(".saveInitialTime", "numeric", NA, NA, NA,
                    "This describes the simulation time at which the first save event should occur"),
    defineParameter(".saveInterval", "numeric", NA, NA, NA,
                    "This describes the simulation time interval between save events"),
    defineParameter(".seed", "list", NULL, NA, NA,
                    paste("Named list of seeds to use for each event (names). E.g., `list('init' = 123)` will `set.seed(123)`",
                          "at the start of the init event and unset it at the end. Defaults to `NULL`, meaning that",
                          "no seeds will be set")),
    defineParameter(".sslVerify", "integer", as.integer(unname(curl::curl_options("^ssl_verifypeer$"))), NA_integer_, NA_integer_,
                    paste("Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading KNN",
                          "(NFI) datasets. Set to 0L if necessary to bypass checking the SSL certificate (this",
                          "may be necessary when NFI's website SSL certificate is not correctly configured).")),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If `NA`, a hash of studyArea will be used."),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    desc = "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
  ),
  inputObjects = bindrows(
    expectsInput("cloudFolderID", "character",
                 "The google drive location where cloudCache will store large statistical objects"),
    expectsInput("columnsForPixelGroups", "character",
                 paste("The names of the columns in `cohortData` that define unique `pixelGroup`s.",
                       "Default is `c('ecoregionGroup', 'speciesCode', 'age')`;",
                       "see `?LandR::columnsForPixelGroups`).")),
    expectsInput("ecoregionLayer", "sf",
                 desc = paste("A `sf` polygon object that characterizes the unique ecological regions (`ecoregionGroup`) used to",
                              "parameterize the biomass, cover, and species establishment probability models.",
                              "It will be overlaid with landcover to generate classes for every ecoregion/LCC combination.",
                              "It must have same extent and crs as `studyAreaLarge`.",
                              "It is superseded by `sim$ecoregionRst` if that object is supplied by the user"),
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("ecoregionRst", "SpatRaster",
                 desc = paste("A raster that characterizes the unique ecological regions used to",
                              "parameterize the biomass, cover, and species establishment probability models.",
                              "If this object is provided, it will supercede `sim$ecoregionLayer`.",
                              "It will be overlaid with landcover to generate classes for every ecoregion/LCC combination.",
                              "It must have same extent and crs as `rasterToMatchLarge` if supplied by user - use `reproducible::postProcess`.",
                              "If it uses an attribute table, it must contain the field 'ecoregion' to represent raster values")),
    expectsInput("firePerimeters", "SpatRaster",
                 desc = paste("Fire perimeters raster, with fire year information used to 'update' stand",
                              "age using time since last fire as the imputed value. Only used if",
                              "`P(sim)$overrideAgeInFires = TRUE`. Biomass will also be updated in these pixels",
                              "if `P(sim)$overrideBiomassInFires = TRUE` and the last fire was later than 1985.",
                              "Defaults to using fire perimeters in the Canadian National Fire Database, downloaded",
                              "as a zipped shapefile with fire polygons, an attribute (i.e., a column) named 'YEAR',",
                              "which is used to rasterize to the study area."),
                 sourceURL = "https://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip"),
    expectsInput("imputedPixID", "integer",
                 desc = paste("A vector of pixel IDs - matching rasterMatch IDs - that suffered data imputation.",
                              "Data imputation may be in age (to match last fire event post 1950s, or 0 cover),",
                              "biomass (to match fire-related imputed ages; correct for missing values or for 0 age/cover),",
                              "land cover (to convert non-forested classes into to nearest forested class).",
                              "If `standAgeMap` had imputed data, then this is expected to be created at that time.",
                              " It will be added as an attribute to `sim$standAgeMap`"),
                 sourceURL = NA),
    expectsInput("rstLCC", "SpatRaster",
                 desc = paste("A land classification map in study area. It must be 'corrected', in the sense that:\n",
                              "1) Every class must not conflict with any other map in this module\n",
                              "    (e.g., `speciesLayers` should not have data in LCC classes that are non-treed);\n",
                              "2) It can have treed and non-treed classes. The non-treed will be removed within this\n",
                              "    module if `P(sim)$omitNonTreedPixels` is `TRUE`;\n",
                              "3) It can have transient pixels, such as 'young fire'. These will be converted to a\n",
                              "    the nearest non-transient class, probabilistically if there is more than 1 nearest\n",
                              "    neighbour class, based on `P(sim)$LCCClassesToReplaceNN`.\n",
                              "The default layer used, if not supplied, is Canada national land classification in 2010.",
                              " The metadata (res, proj, ext, origin) need to match `rasterToMatchLarge`."),
                 sourceURL = NA), ## uses P(sim)$rstLCCYear and LandR::prepInputsLCC() defaults
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = paste("A raster of the `studyArea` in the same resolution and projection as `rawBiomassMap`.",
                              "This is the scale used for all *outputs* for use in the simulation.",
                              "If not supplied will be forced to match the *default* `rawBiomassMap`.")),
    expectsInput("rasterToMatchLarge", "SpatRaster",
                 desc = paste("A raster of the `studyAreaLarge` in the same resolution and projection as `rawBiomassMap`.",
                              "This is the scale used for all *inputs* for use in the simulation.",
                              "If not supplied will be forced to match the *default* `rawBiomassMap`.")),
    expectsInput("rawBiomassMap", "SpatRaster",
                 desc = paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
                              "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
                              "from 2001 (in tonnes/ha), unless 'dataYear' != 2001. See",
                              "https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990",
                              "for metadata."),
                 sourceURL = paste0("https://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Biomass_TotalLiveAboveGround_v1.tif")),
    expectsInput("speciesLayers", "SpatRaster",
                 desc = paste("cover percentage raster layers by species in Canada species map.",
                              "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                              "kNN-derived species cover maps from 2001 using a cover threshold of 10 -",
                              "see https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/2001-attributes_attributs-2001/")),
    expectsInput("speciesTable", "data.table",
                 desc = paste("a table of invariant species traits with the following trait colums:",
                              "'species', 'Area', 'longevity', 'sexualmature', 'shadetolerance',",
                              "'firetolerance', 'seeddistance_eff', 'seeddistance_max', 'resproutprob',",
                              "'resproutage_min', 'resproutage_max', 'postfireregen', 'leaflongevity',",
                              "'wooddecayrate', 'mortalityshape', 'growthcurve', 'leafLignin',",
                              "'hardsoft'. Names can differ, but not the column order.",
                              "Default is from Dominic Cyr and Yan Boulanger's project."),
                 sourceURL = "https://raw.githubusercontent.com/dcyr/LANDIS-II_IA_generalUseFiles/master/speciesTraits.csv"),
    expectsInput("sppColorVect", "character",
                 desc = "named character vector of hex colour codes corresponding to each species"),
    expectsInput("sppEquiv", "data.table",
                 desc = "table of species equivalencies. See `?LandR::sppEquivalencies_CA`."),
    expectsInput("sppNameVector", "character",
                 desc = paste("an optional vector of species names to be pulled from `sppEquiv`.",
                              "Species names must match `P(sim)$sppEquivCol` column in `sppEquiv`.",
                              "If not provided, then species will be taken from",
                              "the entire `P(sim)$sppEquivCol` column in `sppEquiv`.",
                              "See `LandR::sppEquivalencies_CA`.")),
    expectsInput("standAgeMap", "SpatRaster",
                 desc =  paste("stand age map in study area. Must have a 'imputedPixID' attribute (a  vector of pixel IDs)",
                               "indicating which pixels suffered age imputation. If no pixel ages were imputed, please set",
                               "this attribute to `integer(0)`.",
                               "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                               "kNN-derived biomass map from 2001, unless 'dataYear' != 2001.",
                               "See https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990 for metadata"),
                 sourceURL = paste0("http://ftp.maps.canada.ca/pub/nrcan_rncan/Forests_Foret/",
                                    "canada-forests-attributes_attributs-forests-canada/",
                                    "2001-attributes_attributs-2001/",
                                    "NFI_MODIS250m_2001_kNN_Structure_Stand_Age_v1.tif")),
    expectsInput("studyArea", "sfc",
                 desc = paste("Polygon to use as the study area. Must be supplied by the user. Can also be a SpatVector.")),
    expectsInput("studyAreaLarge", "sfc",
                 desc = paste("multipolygon (potentially larger than `studyArea`) used for parameter estimation,",
                              "Must be supplied by the user. If larger than `studyArea`, it must fully contain it.",
                              "Can also be a SpatVector."))
  ),
  outputObjects = bindrows(
    createsOutput("biomassMap", "SpatRaster",
                  paste("total biomass raster layer in study area,",
                        "filtered for pixels covered by cohortData. Units in $g/m^2$")),
    createsOutput("cohortData", "data.table",
                  paste("initial community table, containing corrected biomass ($g/m^2$), age and",
                        "species cover data, as well as ecolocation and `pixelGroup` information. This table defines",
                        "the initial community composition and structure used by `Biomass_core`")),
    createsOutput("ecoregion", "data.table",
                  paste("`ecoregionGroup` look up table")),
    createsOutput("ecoregionMap", "SpatRaster",
                  paste("`ecoregionGroup` map that has mapcodes match `ecoregion` table and `speciesEcoregion` table")),
    createsOutput("firePerimeters", "SpatRaster",
                  paste("As the input object `firePerimeters`, but potentially cropped/masked/projected to match `rasterToMatchLarge`")),
    createsOutput("imputedPixID", "integer",
                  paste("A vector of pixel IDs - matching rasterMatch IDs - that suffered data imputation.",
                        "Data imputation may be in age (to match last fire event post 1950s, or 0 cover),",
                        "biomass (to match fire-related imputed ages, correct for missing values or for 0 age/cover),",
                        "land cover (to convert non-forested classes into to nearest forested class)")),
    createsOutput("pixelGroupMap", "SpatRaster",
                  "initial community map that has mapcodes (`pixelGroup` IDs) match `cohortData`"),
    createsOutput("pixelFateDT", "data.table",
                  paste("A small table that keeps track of the pixel removals and cause.",
                        "This may help diagnose issues related to understanding the creation of `cohortData`.")),
    createsOutput("minRelativeB", "data.frame",
                  paste("minimum relative biomass thresholds that determine a shade level in each",
                        "pixel. `X0-5` represent site shade classes from no-shade (0) to maximum shade (5).")),
    createsOutput("modelCover", "data.frame",
                  paste("If `P(sim)$exportModels` is 'all', or 'cover',",
                        "fitted cover model, as defined by `P(sim)$coverModel`.")),
    createsOutput("modelBiomass", "data.frame",
                  paste("If `P(sim)$exportModels` is 'all', or 'biomass',",
                        "fitted biomass model, as defined by `P(sim)$biomassModel`")),
    # createsOutput("rawBiomassMap", "SpatRaster",
    #               paste("total biomass raster layer in study area. Defaults to the Canadian Forestry",
    #                     "Service, National Forest Inventory, kNN-derived total aboveground biomass map",
    #                     "(in tonnes/ha) from 2001, unless `dataYear != 2001`.",
    #                     "See <https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990>",
    #                     "for metadata")),
    createsOutput("rstLCC", "SpatRaster",
                  paste("As the input object `rstLCC`, but potentially cropped/projected/masked",
                        "to match `rasterToMatchLarge`")),
    createsOutput("species", "data.table",
                  paste("Table that of invariant species traits.",
                        "Will have the same traits as the input `speciesTable`,",
                        "with values adjusted where necessary.")),
    createsOutput("speciesLayers", "SpatRaster",
                  paste("cover percentage raster layers by species in Canada species map.",
                        "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                        "kNN-derived species cover maps from 2001 using a cover threshold of 10 -",
                        "see <https://open.canada.ca/data/en/dataset/ec9e2659-1c29-4ddb-87a2-6aced147a990>",
                        "for metadata.")),
    createsOutput("speciesEcoregion", "data.table",
                  paste("table of spatially-varying species traits (`maxB`, `maxANPP`, `establishprob`),",
                        "defined by species and `ecoregionGroup` (i.e. ecolocation)")),
    createsOutput("standAgeMap", "SpatRaster",
                  paste("As the input object `standAgeMap`, but potentially cropped, projected,",
                        "masked to match `rasterToMatchLarge`.")),
    createsOutput("studyArea", "sfc",
                  paste("As the input object `studyArea`, but potentially projected to match `rasterToMatch` CRS.")),
    createsOutput("studyAreaLarge", "sfc",
                  paste("As the input object `studyAreaLarge`, but potentially projected to match `studyArea`",
                        "and `rasterToMatch` CRS.")),
    createsOutput("sufficientLight", "data.frame",
                  paste("Probability of germination for species shade tolerance (in `species`)",
                        "and shade level`(defined by `minRelativeB`) combinations.",
                        "Table values follow LANDIS-II test traits available at:",
                        paste0("<https://raw.githubusercontent.com/LANDIS-II-Foundation/",
                               "Extensions-Succession/master/biomass-succession-archive/",
                               "trunk/tests/v6.0-2.0/biomass-succession_test.txt>")))
  )
))

## event types
#   - type `init` is required for initialiazation

doEvent.Biomass_borealDataPrep <- function(sim, eventTime, eventType, debug = FALSE) {

  ## open a plotting device so that Biomass_core doesn't plot on top of it if it's too small.
  ## needs to be outside of init, in case init event is cached.
  if (anyPlotting(P(sim)$.plots) && any("screen" %in% P(sim)$.plots)) {
    dev()
    clearPlot()
    mod$plotWindow <- dev.cur()
  }

  switch(
    eventType,
    init = {
      sim <- createBiomass_coreInputs(sim)

      # schedule future event(s)
      sim <- scheduleEvent(sim, P(sim)$.saveInitialTime, "Biomass_borealDataPrep", "save")

      if (anyPlotting(P(sim)$.plots)) {
        plottingFn(sim)
      }
    },
    save = {
      sim <- Save(sim)
    },
    warning(paste("Undefined event type: '", current(sim)[1, "eventType", with = FALSE],
                  "' in module '", current(sim)[1, "moduleName", with = FALSE], "'", sep = ""))
  )
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
  if (is.null(sim$speciesLayers)) {
    stop(red(paste(
      "'speciesLayers' are missing in Biomass_borealDataPrep init event.\n",
      "This is likely due to the module producing 'speciesLayers' being scheduled after Biomass_borealDataPrep.\n",
      "Please check module order."
    )))
  }

  if (!all(P(sim)$LCCClassesToReplaceNN %in% P(sim)$forestedLCCClasses)) {
    stop("All 'LCCClassesToReplaceNN' should be included in 'forestedLCCClasses'.")
  }

  ## check that input rasters all match
  # Too many times this was failing with non-Terra # Eliot March 8, 2022
  # Now it fails with terra: Ceres Jul 08 2022
  # opt <- options("reproducible.useTerra" = FALSE)
  # on.exit(options(opt), add = TRUE)
  if (!.compareRas(sim$standAgeMap, sim$rasterToMatchLarge, res = TRUE)) {
    ## note that extents may never align if the resolution and projection do not allow for it
    ## this is not working, need to use projectRaster
    sim$standAgeMap <- Cache(postProcess,
                             sim$standAgeMap,
                             to = sim$rasterToMatchLarge,
                             overwrite = TRUE)
    attr(sim$standAgeMap, "imputedPixID") <- sim$imputedPixID
  }

  if (!.compareRas(sim$rstLCC, sim$rasterToMatchLarge, res = TRUE)) {
    sim$rstLCC <- Cache(postProcess,
                        sim$rstLCC,
                        to = sim$rasterToMatchLarge,
                        overwrite = TRUE)
  }

  if (P(sim)$overrideAgeInFires) {
    if (!.compareRas(sim$firePerimeters, sim$rasterToMatchLarge, res = TRUE, stopOnError = FALSE)) {
      sim$firePerimeters <- Cache(postProcess,
                                  sim$firePerimeters,
                                  to = sim$rasterToMatchLarge,
                                  overwrite = TRUE)
    }
  }
  # options(opt)
  if (!.compareRas(sim$speciesLayers, sim$rasterToMatchLarge, res = TRUE)) {
    sim$speciesLayers <- Cache(postProcessTerra,
                               sim$speciesLayers,
                               to = sim$rasterToMatchLarge,
                               overwrite = TRUE)
  }

  if (!.compareRas(sim$rasterToMatchLarge, sim$rawBiomassMap, sim$rstLCC,
                   sim$speciesLayers, sim$standAgeMap, res = TRUE)) {
    stop("sim$rasterToMatchLarge, sim$rawBiomassMap, sim$rstLCC,
                   sim$speciesLayers, sim$standAgeMap properties do not match")
  }

  ## species traits inputs ---------------------------------------
  message(blue("Prepare 'species' table, i.e., species level traits", Sys.time()))

  sim$species <- Cache(prepSpeciesTable(speciesTable = sim$speciesTable,
                                        sppEquiv = sim$sppEquiv,
                                        areas = P(sim)$speciesTableAreas,
                                        sppEquivCol = P(sim)$sppEquivCol))

  ## override species table values -------------------------------
  if (!is.null(P(sim)$speciesUpdateFunction)) {
    for (fn in P(sim)$speciesUpdateFunction) {
      if (is(fn, "call")) {
        sim$species <- eval(fn)
      } else {
        stop("speciesUpdateFunction should be a list of one or more quoted function expressions e.g.:\n",
             "list(quote(LandR::speciesTableUpdate(...)), quote(speciesTableUpdateCustom(...)))")
      }
    }
  }

  if (getOption("LandR.verbose") > 0) {
    message("Adjusting species-level traits, part 2")
    print(sim$species)
  }

  ## check that all species have trait values.
  missingTraits <- setdiff(names(sim$speciesLayers), sim$species$species)
  if (length(missingTraits) == length(names(sim$speciesLayers))) {
    stop("No trait values were found for ", paste(missingTraits, collapse = ", "), ".\n",
         "Please check the species list and traits table")
  } else if (length(missingTraits)) {
    spps <- grep("_Spp", missingTraits, ignore.case = TRUE)
    if (length(spps)) {
      toRm <- grep("_Spp", names(sim$speciesLayers))
      sim$speciesLayers <- sim$speciesLayers[[-toRm]] # works on Raster or SpatRaster
      message("No trait values were found for ", paste(missingTraits, collapse = ", "), ".\n",
              " Since this is a 'genus-level' designation (_Spp), omitting it. ",
              " Please ensure that is the correct behaviour")
    }
    missingTraits <- setdiff(names(sim$speciesLayers), sim$species$species)
    if (length(missingTraits))
      stop("No trait values were found for ", paste(missingTraits, collapse = ", "), ".\n",
           "Missing traits will result in species removal from simulation.\n
            Please check the species list and traits table")
  }

  ## filter table in case sppEquiv has more species tahn those being modelled
  sim$species <- sim$species[species %in% names(sim$speciesLayers)]

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

  ## initialEcoregionMap -----------------------------------------
  if (!.compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- projectTo(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }

  if (!.compareCRS(sim$studyAreaLarge, sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- projectTo(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  ## Clean pixels for veg. succession model
  ## remove pixels with no species data or non-forested LCC
  pixelsToRm <- nonForestedPixels(sim$speciesLayers, P(sim)$omitNonTreedPixels,
                                  P(sim)$forestedLCCClasses, sim$rstLCC)
  rstLCCAdj <- sim$rstLCC
  rstLCCAdj[pixelsToRm] <- NA

  pixelFateDT <- pixelFate(fate = "Total number pixels", runningPixelTotal = ncell(sim$speciesLayers))
  pixelFateDT <- pixelFate(pixelFateDT, "NAs on sim$speciesLayers", sum(pixelsToRm))
  if (P(sim)$omitNonTreedPixels) {
    pixelFateDT <- pixelFate(pixelFateDT, "Non forested pixels (based on LCC classes)",
                             sum(!(as.vector(sim$rstLCC[]) %in% P(sim)$forestedLCCClasses)) -
                               tail(pixelFateDT$pixelsRemoved, 1))
  }

  ## The next function will remove the "zero" class on sim$ecoregionRst
  pixelFateDT <- pixelFate(pixelFateDT, "Removing 0 class in sim$ecoregionRst",
                           sum(as.vector(sim$ecoregionRst[])[!pixelsToRm] == 0, na.rm = TRUE))
  ecoregionFiles <- prepEcoregions(
    ecoregionRst = sim$ecoregionRst,
    ecoregionLayer = sim$ecoregionLayer,
    ecoregionLayerField = P(sim)$ecoregionLayerField,
    rasterToMatchLarge = sim$rasterToMatchLarge,
    rstLCCAdj = rstLCCAdj,
    pixelsToRm = pixelsToRm,
    cacheTags = c(cacheTags, "prepEcoregionFiles")
  ) |>
    Cache()

  ## create pixelTable object ------------------------------------
  #  Round age to pixelGroupAgeClass
  # Internal data.table is changed; using memoise here causes the internal changes to
  #   come out to the pixelTable, which is not desired. Turn off memoising for one step
  opt <- options("reproducible.useMemoise" = FALSE)
  on.exit(try(options(opt), silent = TRUE), add = TRUE)

  pixelTable <- makePixelTable(
    speciesLayers = sim$speciesLayers,
    standAgeMap = sim$standAgeMap,
    ecoregionFiles = ecoregionFiles,
    biomassMap = sim$rawBiomassMap,
    rasterToMatch = sim$rasterToMatchLarge,
    rstLCC = rstLCCAdj
  ) |>
    Cache(userTags = c(cacheTags, "pixelTable"), omitArgs = c("userTags"))
  options(opt)
  pixelTable[, rasterToMatch := NULL]

  ## create initial pixelCohortData table ---------------
  coverColNames <- paste0("cover.", sim$species$species)
  pixelCohortData <- makeAndCleanInitialCohortData(
    inputDataTable = pixelTable,
    sppColumns = coverColNames,
    imputeBadAgeModel = P(sim)$imputeBadAgeModel,
    minCoverThreshold = P(sim)$minCoverThreshold,
    doSubset = P(sim)$subsetDataAgeModel
  ) |>
    Cache(userTags = c(cacheTags, "pixelCohortData"), omitArgs = c("userTags"))
  assertCohortDataAttr(pixelCohortData)

  sim$imputedPixID <- unique(c(sim$imputedPixID, attr(pixelCohortData, "imputedPixID")))
  pixelFateDT <- pixelFate(pixelFateDT, "makeAndCleanInitialCohortData rm cover < minThreshold",
                           tail(pixelFateDT$runningPixelTotal, 1) -
                             NROW(unique(pixelCohortData$pixelIndex)))

  ## partition totalBiomass into individual species B -----------------------------------------
  ## via estimating how %cover and %biomass are related
  message(blue("Partitioning totalBiomass per pixel into cohort B as:"))
  if (isTRUE(P(sim)$fitDeciduousCoverDiscount)) {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))

    # pixelCohortData[, lcc := as.factor(lcc)]
    #
    # plot.it <- FALSE
    # sam <- subsetDT(pixelCohortData, by = c("speciesCode", "lcc"),
    #                 doSubset = P(sim)$subsetDataAgeModel,
    #                 indices = TRUE)
    # pi <- unique(pixelCohortData[sam]$pixelIndex)
    # sam <- which(pixelCohortData$pixelIndex %in% pi)
    #
    # system.time({
    #   out <- optimize(interval = c(0.1, 1), f = coverOptimFn, bm = P(sim)$coverPctToBiomassPctModel,
    #                   pixelCohortData = pixelCohortData, subset = sam, maximum = FALSE)
    # })
    # params(sim)$Biomass_borealDataPrep$deciduousCoverDiscount <- out$minimum
    #
    # if (plot.it) {
    #   cover2BiomassModel <- coverOptimFn(out$minimum, pixelCohortData, P(sim)$subsetDataAgeModel,
    #                                      P(sim)$coverPctToBiomassPctModel, returnAIC = FALSE)
    #   sam1 <- sample(NROW(pixelCohortData), 1e5)
    #   dev()
    #   par(mfrow = c(1,2))
    #   plot(predict(cover2BiomassModel$modelBiomass1$mod,
    #                newdata = cover2BiomassModel$pixelCohortData[sam1]),
    #        log(cover2BiomassModel$pixelCohortData$B / 100)[sam1], pch = ".")
    #   abline(a = 0, b = 1)
    #
    #   cover2BiomassModel1 <- coverOptimFn(1, pixelCohortData, P(sim)$subsetDataAgeModel,
    #                                       P(sim)$coverPctToBiomassPctModel,
    #                                       returnAIC = FALSE)
    #   dev()
    #   plot(predict(cover2BiomassModel1$modelBiomass1$mod,
    #                newdata = cover2BiomassModel1$pixelCohortData[sam1]),
    #        log(cover2BiomassModel1$pixelCohortData$B / 100)[sam1], pch = ".")
    #   abline(a = 0, b = 1)
    #
    #   pcd <- pixelCohortData
    #   bb <- pcd[sample(sam)]
    #   cc <- bb[, cover3 := cover * c(1, out$minimum)[decid + 1]][
    #     , actualX := cover3 / sum(cover3) / (cover / 100), by = "pixelIndex"]
    #   setkey(cc, pixelIndex)
    #   mean(cc[speciesCode == "Popu_Tre"]$actualX)
    # }

    params(sim)$Biomass_borealDataPrep$deciduousCoverDiscount <- Cache(deciduousCoverDiscountFun,
                                                                       pixelCohortData = pixelCohortData,
                                                                       coverPctToBiomassPctModel = P(sim)$coverPctToBiomassPctModel,
                                                                       subsetDataAgeModel = P(sim)$subsetDataAgeModel,
                                                                       userTags = c(cacheTags, "decidCoverDisc"),
                                                                       omitArgs = c("userTags"))

  } else {
    message(magenta(paste0(format(P(sim)$coverPctToBiomassPctModel, appendLF = FALSE))))
    message(blue("using previously estimated deciduousCoverDiscount:",
                 round(P(sim)$deciduousCoverDiscount, 3)))
  }

  pixelCohortData <- Cache(partitionBiomass(x = P(sim)$deciduousCoverDiscount, pixelCohortData))
  set(pixelCohortData, NULL, "B", asInteger(pixelCohortData$B/P(sim)$pixelGroupBiomassClass) *
        P(sim)$pixelGroupBiomassClass)
  set(pixelCohortData, NULL, "cover", asInteger(pixelCohortData$cover))

  ## replace unwanted LCC classes ----------------------------------------------------------------
  ## replace 34 and 35 and 36 values -- burns and cities -- to a neighbour class *that exists*.
  ## 1. We need to have a spatial estimate of maxBiomass everywhere there is forest;
  ## we can't have gaps. The pixels that are 34, 35 or 36 are places for which we don't want
  ## maxBiomass associated with their LCC ... i.e., we don't want a maximum
  ## biomass associated with 34 and 35 because those classes are transient.
  ## They will transition to another class before they arrive at a tree maximum biomass.
  ## So, 34 and 35 should not have estimates of maxBiomass 36 is urban.
  ## So, we SHOULD remove these pixels from our studies, except if we are doing
  ## NRV studies (e.g., LandWeb wanted to replace 36 with some forest class) We
  ## decided that we should not use 34 or 35 in our models of Biomass because the
  ## only objective of these models is to estimate maxBiomass, so don't use 34 or
  ## 35 To associate the pixels that were 34 or 35 with a maxBiomass , we need to
  ## give them a "forest class" that they might "become" after they grow out of
  ## being 34 or 35. The pixels where there were 34 and 35 nevertheless have
  ## Biomass estimates in them from KNN and other sources. We leave those as is.
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

  sim$imputedPixID <- unique(c(sim$imputedPixID, newLCCClasses$pixelIndex))
  ## split pixelCohortData into 2 parts -- one with the former 34:36 pixels, one without
  #    The one without 34:36 can be used for statistical estimation, but not the one with
  cohortData34to36 <- pixelCohortData[pixelIndex %in% newLCCClasses$pixelIndex]
  cohortData34to36 <- merge(newLCCClasses, cohortData34to36, all.x = TRUE,
                            all.y = FALSE, by = "pixelIndex")
  cohortDataNo34to36 <- pixelCohortData[!pixelIndex %in% newLCCClasses$pixelIndex]
  if (!length(P(sim)$LCCClassesToReplaceNN)) {
    if (!all.equal(cohortDataNo34to36, pixelCohortData, check.attributes = FALSE))
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

  ## Statistical estimation of establishprob, maxB and maxANPP ----------------------
  cohortDataShort <- cohortDataNo34to36[, list(coverPres = sum(cover > 0)),
                                        by = c("ecoregionGroup", "speciesCode")]
  ## find coverNum for each known class
  ## add new ecoregions to pixelTable, before calc. table
  cohortDataShortNoCover <-
    (function(x) {
      tempDT <- rbind(cohortData34to36[, .(pixelIndex, ecoregionGroup)],
                      cohortDataNo34to36[, .(pixelIndex, ecoregionGroup)])
      pixelTable <- tempDT[pixelTable, on = .(pixelIndex)]

      aa <- table(as.character(pixelTable$ecoregionGroup))   ## as.character avoids counting levels that don't exist anymore

      dt1 <- data.table(ecoregionGroup = factor(names(aa)), coverNum = as.integer(unname(aa)))
      allCombos <- expand.grid(ecoregionGroup = dt1$ecoregionGroup, speciesCode = unique(cohortDataShort$speciesCode))
      setDT(allCombos)
      dt1 <- dt1[allCombos, on = "ecoregionGroup", nomatch = 0]
      cohortDataShortNoCover <- cohortDataShort[dt1, on = c("ecoregionGroup", "speciesCode"), nomatch = NA]
    })() |>
    Cache(.functionName = "cohortDataShortNoCover",
          .cacheExtra = list(cohortData34to36[, .(pixelIndex, ecoregionGroup)],
                             cohortDataNo34to36[, .(pixelIndex, ecoregionGroup)],
                             pixelTable,
                             cohortDataShort))

  # cohortDataShortNoCover <- cohortDataShort[coverPres == 0]
  cohortDataShort <- cohortDataShortNoCover[coverPres > 0] # remove places where there is 0 cover
  cohortDataShortNoCover <- cohortDataShortNoCover[is.na(coverPres)][, coverPres := 0]
  ##  will be added back as establishprob = 0

  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert2(cohortDataShort, classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(cohortDataShortNoCover, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }

  message(blue("Estimating Species Establishment Probability using P(sim)$coverModel, which is"))
  message(magenta(paste0(format(P(sim)$coverModel, appendLF = FALSE), collapse = "")))

  useCloud <- if (!is.null(sim$cloudFolderID)) {
    (isTRUE(getOption("reproducible.useCache", FALSE)) && P(sim)$useCloudCacheForStats)
  } else {
    FALSE
  }

  ## Remove all cases where there is 100% presence in an ecoregionGroup -- causes failures in binomial models
  cdsWh <- cohortDataShort$coverPres == cohortDataShort$coverNum
  cds <- Copy(cohortDataShort)
  cds <- cds[!cdsWh]

  modelCover <- Cache(
    statsModel,
    modelFn = P(sim)$coverModel,
    # modelFn = cm,
    uniqueEcoregionGroups = .sortDotsUnderscoreFirst(as.character(unique(cohortDataShort$ecoregionGroup))),
    sumResponse = sum(cohortDataShort$coverPres, cohortDataShort$coverNum, na.rm = TRUE),
    .specialData = cds,
    .cacheExtra = levels(cohortDataShort$speciesCode), # in case sppEquivCol changes
    useCloud = useCloud,
    cloudFolderID = sim$cloudFolderID,
    # useCache = "overwrite",
    showSimilar = getOption("reproducible.showSimilar", FALSE),
    userTags = c(cacheTags, "modelCover"),
    omitArgs = c("showSimilar", "useCache", ".specialData", "useCloud", "cloudFolderID")
  )
  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelCover$rsq, 4))), function(x) {
    message(blue(x))
  })

  ## export model before overriding happens
  if (any(P(sim)$exportModels %in% c("all", "coverModel"))) {
    sim$modelCover <- modelCover
  }

  if (isTRUE(any(cdsWh))) {
    cds[, pred := fitted(modelCover$mod, response = "response")]
    cohortDataShort <- cds[, -c("coverPres", "coverNum")][cohortDataShort,
                                                          on = c("ecoregionGroup", "speciesCode"), nomatch = NA]
    cohortDataShort[is.na(pred), pred := 1]
    modelCover <- cohortDataShort$pred
  }

  ## For biomass
  ### Subsample cases where there are more than 50 points in an ecoregionGroup * speciesCode
  totalBiomass <- sum(cohortDataNo34to36Biomass$B, na.rm = TRUE)

  ## There are several reasons why the modelBiomass can fail;
  ##   1) inappropriate sub-sample
  ##   2) fit algorithm
  ## Run two nested loops to do both of these things
  modelBiomassTags <- c(cacheTags, "modelBiomass",
                        paste0("subsetSize:", P(sim)$subsetDataBiomassModel))
  maxDataSubsetTries <- ifelse(isTRUE(P(sim)$subsetDataBiomassModel > 0),
                               P(sim)$subsetDataAttempts, 1)
  for (tryBiomassDataSubset in 1:maxDataSubsetTries) {
    cohortDataNo34to36BiomassSubset <- subsetDT(cohortDataNo34to36Biomass,
                                                by = c("ecoregionGroup", "speciesCode"),
                                                doSubset = P(sim)$subsetDataBiomassModel)

    ### For Cache -- doesn't need to cache all columns in the data.table -- only the ones in the model
    ### force parameter values to avoid more checks
    # If using mixed effect model, see here for good discussion of
    #  shrinkage https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
    message(blue("Estimating biomass using P(sim)$biomassModel as:"), "\n",
            magenta(paste0(format(P(sim)$biomassModel, appendLF = FALSE), collapse = "")))

    ## NOTE: we are NOT using logB because the relationship between B~age should be hump-shaped
    ## (or at least capped at high age values). Ideally, we would want a non-linear model

    # Default values of args to modelBiomass -- prior to any attempts to fix
    ueg <- .sortDotsUnderscoreFirst(as.character(unique(cohortDataNo34to36BiomassSubset$ecoregionGroup)))
    specDat <- cohortDataNo34to36BiomassSubset
    modelFn <- P(sim)$biomassModel
    sumResponse <- c(totalBiomass)

    fixModelBiomass <- P(sim)$fixModelBiomass
    timePriorToFit <- Sys.time()
    cohortDataNo34to36BiomassSubset2 <- copy(cohortDataNo34to36BiomassSubset)

    tryControl <- FALSE
    needRescaleModelB <- FALSE
    scaledVarsModelB <- NULL
    for (tryBiomassModel in 1:3) { # try thrice -- default, then once to rescale, once to refit
      modelBiomass <- Cache(
        statsModel,
        modelFn = modelFn,
        uniqueEcoregionGroups = ueg,
        .cacheExtra = sumResponse, # only digest on this (formerly used sumResponse arg; now .cacheExtra is in Cache)
        .specialData = specDat,
        useCloud = useCloud,
        # useCache = "overwrite",
        cloudFolderID = sim$cloudFolderID,
        # showSimilar = getOption("reproducible.showSimilar", FALSE),
        userTags = c(modelBiomassTags,
                     paste0("subsetSize:", P(sim)$subsetDataBiomassModel)),
        omitArgs = c("showSimilar", ".specialData", "useCloud", "cloudFolderID", "useCache")
      )

      modMessages <- modelBiomass$mod@optinfo$conv$lme4$messages
      needRedo <- (length(modMessages) > 0 & fixModelBiomass)
      if (needRedo && (!tryControl || !needRescaleModelB)) {
        modCallChar <- paste(deparse(P(sim)$biomassModel), collapse = "")
        if (any(grepl("Rescale", modMessages)) & !needRescaleModelB) {
          message(blue("Trying to rescale variables to refit P(sim)$biomassModel"))
          ## save this in separate objects for later
          logAge_sc <- scale(cohortDataNo34to36BiomassSubset$logAge)
          cover_sc <- scale(cohortDataNo34to36BiomassSubset$cover)

          scaledVarsModelB <- list(logAge = logAge_sc, cover = cover_sc)
          ## remove attributes with as.numeric
          ## don't change the original data
          cohortDataNo34to36BiomassSubset2[, `:=`(logAge = as.numeric(logAge_sc),
                                                  cover = as.numeric(cover_sc))]
          needRescaleModelB <- TRUE
          ueg <- .sortDotsUnderscoreFirst(as.character(unique(cohortDataNo34to36BiomassSubset2$ecoregionGroup)))
        } else {
          message(blue("Trying to refit P(sim)$biomassModel with 'bobyqa' optimizer"))
          ## redo model call with new optimizer
          modCallChar <- paste(deparse(P(sim)$biomassModel), collapse = "")
          if (grepl("lme4::lmer", modCallChar)) {
            modCallChar <-  sub(")$", ", control = lme4::lmerControl(optimizer = 'bobyqa'))", modCallChar)
          } else if (grepl("lme4::glmer", modCallChar)) {
            modCallChar <-  sub(")$", ", = lme4::glmerControl(optimizer = 'bobyqa'))", modCallChar)
          } else {
            message(blue("P(sim)$biomassModel does not call 'lme4::lmer' or 'lme4::glmer' explicitly",
                         "preventing an attempt to use a different optimizer."))
          }
          tryControl <- TRUE
        }
        userTagsToClear <- c("statsModel", modelBiomassTags[1:3])
        suppressMessages(clearCache(userTags = userTagsToClear, #after = timePriorToFit,
                                    ask = FALSE))
        specDat <- cohortDataNo34to36BiomassSubset2
        modelBiomassTags <- c("refit", "modelBiomass",
                              paste(c(if (needRescaleModelB) "rescaled",
                                      if (tryControl) "control"), collapse = "_"))
        modelFn <- str2lang(modCallChar)
        sumResponse <- c(totalBiomass, tryControl, needRescaleModelB)

        ## break out of while, even after trying to rescale and fit with bobyqa
        if (needRescaleModelB & tryControl)
          fixModelBiomass <- FALSE
      } else {
        if (tryControl && needRescaleModelB)
          warning("Biomass model did not converge and automated attempts to fix also failed.",
                  " This will need more attention.")
        break
      }
    } ## End of tryBiomassModel
    if (!needRedo)
      break
  } ## End of tryBiomassData

  if (!is.null(scaledVarsModelB)) {
    modelBiomass$scaledVarsModelB <- scaledVarsModelB
  }

  if (isTRUE(tryBiomassDataSubset == maxDataSubsetTries) && isTRUE(needRedo)) {
    warning("The biomass model did not converge with ", tryBiomassDataSubset,
            " attempts of data subsetting and changing lme algorithm.")
  }

  message(blue("  The rsquared is: "))
  out <- lapply(capture.output(as.data.frame(round(modelBiomass$rsq, 4))), function(x) {
    message(blue(x))
  })

  if (any(P(sim)$exportModels %in% c("all", "biomassModel")))
    sim$modelBiomass <- modelBiomass

  ## remove logB
  # cohortDataNo34to36BiomassSubset[, logB := NULL]

  ## create speciesEcoregion ---------------------------------------------
  ## a single line for each combination of ecoregionGroup & speciesCode;
  ## doesn't include combinations with B = 0 because those places can't have the species/ecoregion combo
  ## cohortDataNo34to36BiomassSubset ends up determining which ecoregion combinations end up in
  ## species ecoregion, thus removing converted/masked classes present cohortDataShortNoCover
  message(blue("Create speciesEcoregion using modelCover and modelBiomass to estimate species traits"))
  speciesEcoregion <- makeSpeciesEcoregion(cohortDataBiomass = cohortDataNo34to36BiomassSubset,
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

  # if (!is.na(P(sim)$.plotInitialTime)) {
  #   uniqueSpeciesNames <- as.character(unique(speciesEcoregion$speciesCode))
  #   names(uniqueSpeciesNames) <- uniqueSpeciesNames
  #   speciesEcoregionTable2 <- copy(speciesEcoregion)
  #   speciesEcoregionTable2[, ecoregionInt := as.integer(ecoregionGroup)]
  #   maxB <- raster::stack(lapply(uniqueSpeciesNames, function(sp) {
  #     rasterizeReduced(speciesEcoregionTable2[speciesCode == sp], ecoregionFiles$ecoregionMap,
  #                      "maxB", "ecoregionInt")
  #   }))
  #   curDev <- dev.cur()
  #   newDev <- if (!is.null(dev.list())) max(dev.list()) + 1 else 1
  #   quickPlot::dev(newDev, width = 18, height = 10)
  #   Plot(maxB, legendRange = c(0, max(maxValue(maxB), na.rm = TRUE)))
  #   quickPlot::dev(curDev)
  # }

  if (ncell(sim$rasterToMatchLarge) > 3e7) replicate(3, gc())

  ## Create initial communities, i.e., pixelGroups -----------------------
  ## Rejoin back the pixels that were 34:36
  set(cohortData34to36, NULL, "initialEcoregionCode", NULL)
  pixelCohortData <- rbindlist(list(cohortData34to36, cohortDataNo34to36),
                               use.names = TRUE, fill = TRUE)

  ## "Downsize" to studyArea after estimating parameters on studyAreaLarge --------------
  ## 1. Subset pixels (IDs) on rasterToMatchLarge, using rasterToMatch
  ## 2. Subset data.tables using the pixel IDs / ecoregion/species combinations
  ##    that are common across the two rasters
  ## 3. Re-do pixel ID numbering so that it matches the final rasterToMatch
  ## Note: if SA and SALarge are the same, no subsetting will take place.
  if (sum(is.na(as.vector(values(sim$rasterToMatch)))) != sum(is.na(as.vector(values(sim$rasterToMatchLarge))))) {
    message(blue("Subsetting to studyArea"))
    rasterToMatchLarge <- sim$rasterToMatchLarge
    rasterToMatchLarge <- setValues(rasterToMatchLarge, seq(ncell(rasterToMatchLarge)))

    rasterToMatchLargeCropped <- Cache(postProcess,
                                       x = rasterToMatchLarge,
                                       to = sim$rasterToMatch,
                                       datatype = assessDataType(rasterToMatchLarge),
                                       method = "near",
                                       userTags = c(cacheTags, "rasterToMatchLargeCropped"),
                                       omitArgs = c("userTags"))

    rtmlc_int <- LandR::asInt(rasterToMatchLargeCropped)
    assertthat::assert_that(all(na.omit(as.vector(rasterToMatchLargeCropped - rtmlc_int)) == 0))
    rm(rtmlc_int)
    assertthat::assert_that(sum(is.na(as.vector(rasterToMatchLargeCropped))) < ncell(rasterToMatchLargeCropped))
    ## i.e., not all NA

    if (!.compareRas(rasterToMatchLargeCropped, sim$rasterToMatch)) {
      stop("Downsizing to rasterToMatch after estimating parameters didn't work.",
           "Please debug Biomass_borealDataPrep::createBiomass_coreInputs().")
    }

    ## subset pixels that are in studyArea/rasterToMatch only
    pixToKeep <- na.omit(as.vector(values(rasterToMatchLargeCropped))) # these are the old indices of RTML
    pixelCohortData <- pixelCohortData[pixelIndex %in% pixToKeep]

    ## re-do pixelIndex (it now needs to match rasterToMatch)
    newPixelIndexDT <- data.table(pixelIndex = as.vector(values(rasterToMatchLargeCropped)),
                                  newPixelIndex = as.integer(1:ncell(rasterToMatchLargeCropped))) |>
      na.omit()

    pixelCohortData <- newPixelIndexDT[pixelCohortData, on = "pixelIndex"]
    pixelCohortData[, pixelIndex := NULL]
    setnames(pixelCohortData, old = "newPixelIndex", new = "pixelIndex")

    assertthat::assert_that(NROW(pixelCohortData) > 0)

    ## now convert imputedPixID to RTM
    sim$imputedPixID <- newPixelIndexDT[pixelIndex %in% sim$imputedPixID, newPixelIndex]

    rm(pixToKeep, rasterToMatchLargeCropped, newPixelIndexDT)
    if (ncell(sim$rasterToMatch) > 3e7) replicate(3, gc())
  }
  ## subset ecoregionFiles$ecoregionMap to smaller area.

  ecoregionFiles$ecoregionMap <- Cache(postProcess,
                                       x = ecoregionFiles$ecoregionMap,
                                       to = sim$rasterToMatch,
                                       writeTo = NULL,
                                       userTags = c(cacheTags, "ecoregionMap"),
                                       omitArgs = c("userTags"))

  if (is(P(sim)$minRelativeBFunction, "call")) {
    sim$minRelativeB <- eval(P(sim)$minRelativeBFunction)
  } else {
    stop("minRelativeBFunction should be a quoted function expression, using `pixelCohortData`, e.g.:\n",
         "    quote(LandR::makeMinRelativeB(pixelCohortData))")
  }

  maxAgeHighQualityData <- -1

  maxRawB <- max(values(sim$rawBiomassMap), na.rm = TRUE) * 100 ## match units in cohortData (t/ha ==> g/m^2)
  # maxRawB <- maxValue(sim$rawBiomassMap) * 100 ## match units in cohortData (t/ha ==> g/m^2)

  ## If this module used a fire database to extract better young ages, then we
  ##   can use those high quality younger ages to help with our biomass estimates

  if (isTRUE(P(sim)$overrideBiomassInFires)) {
    if (isFALSE(P(sim)$overrideAgeInFires)) {
      message(blue("'P(sim)$overrideBiomassInFires' is TRUE but 'P(sim)$overrideAgeInFires' if FALSE."))
      message(blue("B values will NOT be re-estimated inside fire perimeters."))
    } else {
      message(blue("Overriding B values (originally from 'rawBiomassMap') within the fire perimeters",
                   "defined in 'firePerimeters'."))
      message(blue("To skip this step, set 'P(sim)$overrideBiomassInFires' to FALSE."))

      firstFireYear <- min(as.vector(sim$firePerimeters[]), na.rm = TRUE) # 1986 # as.numeric(gsub("^.+nbac_(.*)_to.*$", "\\1", fireURL))
      ## this is not necessary when using min(),
      ## but will be kept in case we use something else in the future
      whichFiresTooOld <- which(as.vector(sim$firePerimeters[]) < firstFireYear)

      if (length(whichFiresTooOld)) {
        message("There were fires in the database older than ", firstFireYear, ";",
                " The data from these will not be used because firstFireYear = ", firstFireYear)
        sim$firePerimeters[whichFiresTooOld] <- NA
      }

      maxAgeHighQualityData <- start(sim) - firstFireYear
      ## if maxAgeHighQualityData is lower than 0, it means it's prior to the first fire Year
      ## or not following calendar year

      if (isTRUE(maxAgeHighQualityData >= 0)) {
        ## identify young in the pixelCohortData
        youngRows <- pixelCohortData$age <= maxAgeHighQualityData
        young <- pixelCohortData[youngRows == TRUE]

        youngRows2 <- !is.na(as.vector(sim$firePerimeters[young$pixelIndex]))
        young <- young[youngRows2]

        # whYoungBEqZero <- which(young$B == 0)
        whYoungZeroToMaxHighQuality <- which(young$age > 0)

        if (length(whYoungZeroToMaxHighQuality) > 0) {
          youngWAgeEqZero <- young[-whYoungZeroToMaxHighQuality]
          youngNoAgeEqZero <- young[whYoungZeroToMaxHighQuality]

          message("Running 'spinup' on pixels that are within fire polygons and whose age < ",
                  maxAgeHighQualityData)
          young <- Cache(spinUpPartial,
                         pixelCohortData = youngNoAgeEqZero,
                         speciesEcoregion = speciesEcoregion,
                         maxAge = maxAgeHighQualityData,
                         minRelativeB = sim$minRelativeB,
                         species = sim$species,
                         sppEquiv = sim$sppEquiv,
                         sppEquivCol = P(sim)$sppEquivCol,
                         sppColorVect = sim$sppColorVect,
                         paths = paths(sim),
                         currentModule = currentModule(sim),
                         modules = modules(sim), ## will also check modules in paths$moduelPath
                         userTags = c(cacheTags, "spinUpYoungBiomasses"),
                         omitArgs = c("userTags", "paths", "modules"))

          ## method using modelBiomass
          ## -- deprecated, as it overestimates B for young ages at the moment
          # young <- Cache(updateYoungBiomasses,
          #                young = youngNoAgeEqZero,
          #                modelBiomass = modelBiomass,
          #                userTags = c(cacheTags, "updateYoungBiomasses"),
          #                omitArgs = c("userTags"))

          if (length(setdiff(colnames(young), colnames(pixelCohortData))) > 0) {
            set(young, NULL, setdiff(colnames(young), colnames(pixelCohortData)), NULL)
          }

          young <- rbindlist(list(young, youngWAgeEqZero), use.names = TRUE)
        } else {
          message(blue("No pixels found with ages needing age replacement with last fire year"))
        }

        lengthUniquePixelIndices <- length(unique(pixelCohortData$pixelIndex))
        pixelCohortData <- rbindlist(list(pixelCohortData[youngRows == FALSE],
                                          pixelCohortData[which(youngRows == TRUE)[!youngRows2]],
                                          young), use.names = TRUE)
        assertthat::assert_that(lengthUniquePixelIndices == length(unique(pixelCohortData$pixelIndex)))

        sim$imputedPixID <- unique(c(sim$imputedPixID, young$pixelIndex))

        ## TODO: reassess 2.8x multiplier; it's high, but needed in RoF_shield
        assertthat::assert_that(
          all(inRange(na.omit(young$B), 0, 2.8 * maxRawB / min(sim$species$longevity/maxAgeHighQualityData)))
        ) # /4 is too strong -- 25 years is a lot of time
      } else {
        ## return maxAgeHighQualityData to -1
        message(blue("Simulation start year is lower than oldest fire."))
        message(blue("B values will NOT be re-estimated inside fire perimeters"))
        maxAgeHighQualityData <- -1
      }
    }
  }

  assertthat::assert_that(all(inRange(na.omit(pixelCohortData$B), 0, round(maxRawB, -2)))) # should they all be below the initial biomass map?

  # Fill in any remaining B values that are still NA -- the previous chunk filled in B for young cohorts only
  if (anyNA(pixelCohortData$B)) {
    theNAsBiomass <- is.na(pixelCohortData$B)
    message(blue(" -- ", sum(theNAsBiomass),"cohort(s) has NA for Biomass: being replaced with model-derived estimates"))
    set(pixelCohortData, which(theNAsBiomass), "B",
        asInteger(predict(modelBiomass$mod, newdata = pixelCohortData[theNAsBiomass],
                          allow.new.levels = TRUE)))
    sim$imputedPixID <- unique(c(sim$imputedPixID, pixelCohortData[theNAsBiomass, pixelIndex]))
  }

  ## make cohortDataFiles: pixelCohortData (rm unnecessary cols, subset pixels with B>0,
  ## generate pixelGroups, add ecoregionGroup and totalBiomass) and cohortData
  cohortDataFiles <- Cache(makeCohortDataFiles,
                           pixelCohortData = pixelCohortData,
                           columnsForPixelGroups = sim$columnsForPixelGroups, # Par$cohortDefinitionCols
                           speciesEcoregion = speciesEcoregion,
                           pixelGroupBiomassClass = P(sim)$pixelGroupBiomassClass,
                           pixelGroupAgeClass = P(sim)$pixelGroupAgeClass,
                           minAgeForGrouping = maxAgeHighQualityData,
                           rmImputedPix = P(sim)$rmImputedPix,
                           imputedPixID = sim$imputedPixID,
                           pixelFateDT = pixelFateDT,
                           userTags = c(cacheTags, "makeCohortData"))

  sim$cohortData <- cohortDataFiles$cohortData
  pixelCohortData <- cohortDataFiles$pixelCohortData
  pixelFateDT <- cohortDataFiles$pixelFateDT

  ## Need to rerun this because we may have lost an Ecoregion_Group in the spinup
  if (is(P(sim)$minRelativeBFunction, "call")) {
    sim$minRelativeB <- eval(P(sim)$minRelativeBFunction)
  } else {
    stop("minRelativeBFunction should be a quoted function expression, using `pixelCohortData`, e.g.:\n",
         "    quote(LandR::makeMinRelativeB(pixelCohortData))")
  }

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

  # if (!is.na(P(sim)$.plotInitialTime)) {
  #   seStacks <- Cache(LandR::speciesEcoregionStack,
  #                     ecoregionMap = sim$ecoregionMap,
  #                     speciesEcoregion = sim$speciesEcoregion,
  #                     columns = c("establishmentprob", "maxB", "maxANPP"),
  #                     stackFilenames = NULL)
  #
  #   sim$ggSpeciesEcoregion <- Map(stk = seStacks, type = names(seStacks),
  #                                    function(stk, type) {
  #                                      ggSpeciesEcoregion <-
  #                                        gplot(stk, maxpixels = 2e6) +
  #                                        geom_tile(aes(fill = value)) +
  #                                        facet_wrap(~ variable) +
  #                                        scale_fill_gradient(low = 'light grey', high = 'blue', na.value = "white") +
  #                                        theme_bw() +
  #                                        coord_equal() +
  #                                        ggtitle(type)
  #                                    }
  #   )
  #   plotList <- lapply(unstack(seStacks$establishprob), plotFunction, studyArea = sim$studyArea)
  #   ggarrange(plotlist = plotList)
  #
  #   Map(gg = sim$ggSpeciesEcoregion, nam = names(sim$ggSpeciesEcoregion),
  #       function(gg, nam) {
  #         png(file.path(outputPath(sim), paste(nam, ".png")), width = 1600, height = 1200)
  #         print(gg)
  #         dev.off()
  #       })
  # }
  #
  sim$pixelGroupMap <- makePixelGroupMap(pixelCohortData, sim$rasterToMatch)

  ## make sure speciesLayers match RTM (since that's what is used downstream in simulations)
  message(blue("Writing sim$speciesLayers to disk as they are likely no longer needed in RAM"))

  # useTerra <- getOption("reproducible.useTerra") ## TODO: reproducible#242
  # options(reproducible.useTerra = FALSE) ## TODO: reproducible#242
  sim$speciesLayers <- Cache(postProcessTo,
                             sim$speciesLayers,
                             to = sim$rasterToMatch,
                             writeTo = .suffix(file.path(outputPath(sim), "speciesLayers.tif"),
                                               paste0("_", P(sim)$dataYear,
                                                      "_", P(sim)$.studyAreaName)),
                             overwrite = TRUE,
                             userTags = c(cacheTags, "speciesLayersRTM", P(sim)$dataYear),
                             quick = "writeTo", # don't digest the file content, just filename
                             # Cache reads file content if it is a file, so it is
                             #    reading content of writeTo, which is an output
                             omitArgs = c("userTags"))
  # options(reproducible.useTerra = useTerra) ## TODO: reproducible#242

  ## double check these rasters all match RTM
  .compareRas(sim$biomassMap, sim$ecoregionMap, sim$pixelGroupMap,
              sim$rasterToMatch, sim$speciesLayers, res = TRUE)

  ## rm ecoregions that may not be present in rasterToMatch
  ## make ecoregionGroup a factor and export speciesEcoregion to sim
  onMatch <- c("ecoregionGroup", "speciesCode")
  toRm <- speciesEcoregion[!sim$cohortData, on = onMatch]
  speciesEcoregion <- speciesEcoregion[!toRm, on = onMatch]
  sim$speciesEcoregion <- speciesEcoregion
  sim$speciesEcoregion$ecoregionGroup <- factor(as.character(sim$speciesEcoregion$ecoregionGroup))

  ## do assertions
  message(blue("Create pixelGroups based on: ", paste(sim$columnsForPixelGroups, collapse = ", ")),
          "\n", blue("Resulted in "), magenta(length(unique(sim$cohortData$pixelGroup))),
          " unique pixelGroup values")
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

plottingFn <- function(sim) {
  ## Step 1 make data
  seStacks <- Cache(LandR:::speciesEcoregionStack,
                    ecoregionMap = sim$ecoregionMap,
                    speciesEcoregion = sim$speciesEcoregion,
                    columns = c("establishprob", "maxB", "maxANPP"),
                    userTags = c("speciesEcoregionStks", P(sim)$.studyAreaName),
                    .cacheExtra = P(sim)$.studyAreaName)

  ## Step 2 make plots -- in this case up to 4 plots -- uses .plotInitialTime, .plots
  if (!is.null(mod$plotWindow)) {
    dev(mod$plotWindow)
  }
  Map(stk = seStacks, SEtype = names(seStacks),
      function(stk, SEtype) {
        Plots(stk,
              fn = plotFn_speciesEcoregion,
              types = P(sim)$.plots,
              filename = paste0("speciesEcoregion", "_", time(sim), "_", SEtype),
              SEtype = SEtype)
      }
  )
}

Save <- function(sim) {
  sim <- saveFiles(sim)
  return(invisible(sim))
}

## see other helper functions in R/ subdirectory

.inputObjects <- function(sim) {
  cacheTags <- c(currentModule(sim), "otherFunctions:.inputObjects")
  dPath <- asPath(inputPath(sim), 1)
  message(currentModule(sim), ": using dataPath '", dPath, "'.")

  # 1. test if all input objects are already present (e.g., from inputs, objects or another module)
  a <- depends(sim)
  whThisMod <- which(unlist(lapply(a@dependencies, function(x) x@name)) == "Biomass_borealDataPrep")
  objNames <- a@dependencies[[whThisMod]]@inputObjects$objectName
  objExists <- !unlist(lapply(objNames, function(x) is.null(sim[[x]])))
  names(objExists) <- objNames

  # for backwards compatibility -- change from parameter to object
  if (!suppliedElsewhere("cloudFolderID", sim)) {
    if (!is.null(P(sim)$cloudFolderID))
      sim$cloudFolderID <- P(sim)$cloudFolderID
  }
  ## Study area(s) ------------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    stop("Please provide a 'studyArea' polygon")
    # message("'studyArea' was not provided by user. Using a polygon (6250000 m^2) in southwestern Alberta, Canada")
    # sim$studyArea <- randomStudyArea(seed = 1234, size = (250^2)*100)  # Jan 2021 we agreed to force user to provide a SA/SAL
  }

  if (!suppliedElsewhere("studyAreaLarge", sim)) {
    stop("Please provide a 'studyAreaLarge' polygon.
         If parameterisation is to be done on the same area as 'studyArea'
         provide the same polygon to 'studyAreaLarge'")
    # message("'studyAreaLarge' was not provided by user. Using the same as 'studyArea'")
    # sim <- objectSynonyms(sim, list(c("studyAreaLarge", "studyArea"))) # Jan 2021 we agreed to force user to provide a SA/SAL
  }

  if (!.compareCRS(sim$studyArea, sim$studyAreaLarge)) {
    warning("studyArea and studyAreaLarge have different projections.\n
            studyAreaLarge will be projected to match crs(studyArea)")
    sim$studyAreaLarge <- projectTo(sim$studyAreaLarge, crs(sim$studyArea))
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

  ## this is necessary if studyArea and studyAreaLarge are multipolygon objects
  if (nrow(studyArea) > 1) {
    studyArea <- st_buffer(studyArea, 0) |> st_union()
  }

  if (nrow(studyAreaLarge) > 1) {
    studyAreaLarge <- st_buffer(studyAreaLarge, 0) |> st_union()
  }

  if (length(st_within(studyArea, studyAreaLarge))[[1]] == 0) {
    stop("studyArea is not fully within studyAreaLarge.
         Please check the aligment, projection and shapes of these polygons")
  }
  rm(studyArea, studyAreaLarge)

  ## Raster(s) to match ------------------------------------------------
  needRTML <- needRTM <- FALSE
  if (is.null(sim$rasterToMatch) || is.null(sim$rasterToMatchLarge)) {
    if (!suppliedElsewhere("rasterToMatch", sim)) {
      needRTM <- TRUE
      message("There is no rasterToMatch supplied; will attempt to use rawBiomassMap")
    }
    if (!suppliedElsewhere("rasterToMatchLarge", sim)) { ## Eliot changed this -- case where RTM was supplied, this broke that --> NOT TRUE --> if one is not provided, redo both (safer?)
      needRTML <- TRUE
      message("There is no rasterToMatchLarge supplied; will use rasterToMatch")
    }
  }

  ## biomass map
  if (!suppliedElsewhere("rawBiomassMap", sim)) {
    if (P(sim)$dataYear %in% c(2001, 2011)) {
      biomassURL <- extractURL("rawBiomassMap") |> gsub("2001", P(sim)$dataYear, x = _)
    } else {
      stop("'P(sim)$dataYear' must be one of 2001 or 2011")
    }

    sim$rawBiomassMap <- prepRawBiomassMap(
      url = biomassURL,
      studyAreaName = P(sim)$.studyAreaName,
      cacheTags = cacheTags,
      to = if (!needRTML) sim$rasterToMatchLarge else if (!needRTM) sim$rasterToMatch else sim$studyAreaLarge,
      projectTo = if (!needRTML) sim$rasterToMatchLarge else if (!needRTM) sim$rasterToMatch else NA, ## don't project to SA if RTMs not present
      destinationPath = dPath)
  }

  if (needRTML || needRTM) {
    if (!is.null(sim$rawBiomassMap)) {
      if (!.compareCRS(sim$rawBiomassMap, sim$studyAreaLarge)) {
        ## note that extents may never align if the resolution and projection do not allow for it
        # opt <- options("reproducible.useTerra" = TRUE) # Too many times this was failing with non-Terra # Eliot March 8, 2022
        # on.exit(options(opt), add = TRUE)
        sim$rawBiomassMap <- Cache(postProcess,
                                   sim$rawBiomassMap,
                                   method = "bilinear",
                                   to = sim$studyAreaLarge,
                                   projectTo = NA,  ## don't project to SA
                                   overwrite = TRUE)
        # options(opt)
      }
    }
  }

  if (needRTML || needRTM) {
    RTMs <- prepRasterToMatch(studyArea = sim$studyArea,
                              studyAreaLarge = sim$studyAreaLarge,
                              rasterToMatch = if (needRTM) NULL else sim$rasterToMatch,
                              rasterToMatchLarge = if (needRTML) NULL else sim$rasterToMatchLarge,
                              destinationPath = dPath,
                              templateRas = sim$rawBiomassMap,
                              studyAreaName = P(sim)$.studyAreaName,
                              cacheTags = cacheTags)
    sim$rasterToMatch <- RTMs$rasterToMatch
    sim$rasterToMatchLarge <- RTMs$rasterToMatchLarge
    rm(RTMs)
  }

  if (!.compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning(paste0("studyArea and rasterToMatch projections differ.\n",
                   "studyArea will be projected to match rasterToMatch"))
    sim$studyArea <- projectInputs(sim$studyArea, crs(sim$rasterToMatch))
    sim$studyArea <- fixErrors(sim$studyArea)
  }
  if (!.compareCRS(sim$studyAreaLarge, sim$rasterToMatchLarge)) {
    warning(paste0("studyAreaLarge and rasterToMatchLarge projections differ.\n",
                   "studyAreaLarge will be projected to match rasterToMatchLarge"))
    sim$studyAreaLarge <- projectInputs(sim$studyAreaLarge, crs(sim$rasterToMatchLarge))
    sim$studyAreaLarge <- fixErrors(sim$studyAreaLarge)
  }

  ## Land cover raster ------------------------------------------------
  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- Cache(prepInputs_NTEMS_LCC_FAO,
                        year = P(sim)$dataYear,
                        maskTo = sim$studyAreaLarge,
                        cropTo = sim$rasterToMatchLarge,
                        projectTo = sim$rasterToMatchLarge,
                        disturbedCode = 240,
                        destinationPath = dPath,
                        overwrite = TRUE,
                        # writeTo = .suffix("rstLCC.tif", paste0("_", P(sim)$.studyAreaName, "_", P(sim)$dataYear)),
                        userTags = c("rstLCC", currentModule(sim),
                                     P(sim)$.studyAreaName, P(sim)$dataYear))
  }

  ## Ecodistrict ------------------------------------------------
  if (!suppliedElsewhere("ecoregionLayer", sim)) {
    ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
    sim$ecoregionLayer <- Cache(prepInputs(targetFile = "ecodistricts.shp",
                                           archive = asPath("ecodistrict_shp.zip"),
                                           url = extractURL("ecoregionLayer", sim),
                                           alsoExtract = "similar",
                                           destinationPath = dPath,
                                           writeTo = NULL,
                                           to = sim$studyAreaLarge,
                                           fun = getOption("reproducible.shapefileRead"),
                                           overwrite = TRUE),
                                .functionName = "prepInputs_forEcoregionLayer",
                                userTags = c("prepInputsEcoDistrict_SA", currentModule(sim), cacheTags))
  }

  if (P(sim)$overrideAgeInFires) {
    if (!suppliedElsewhere("firePerimeters", sim)) {
      # opt <- options("reproducible.useTerra" = TRUE) # Too many times this was failing with non-Terra # Eliot March 8, 2022
      # on.exit(options(opt), add = TRUE)
      sa <- if (is(sim$studyAreaLarge, "sf")) {
        aggregate(sim$studyAreaLarge, list(rep(1, nrow(sim$studyAreaLarge))),
                  FUN = function(x) x)
      } else {
        aggregate(sim$studyAreaLarge)
      }
      sim$firePerimeters <- Cache(
        prepInputsFireYear,
        destinationPath = dPath,
        studyArea = sa,
        rasterToMatch = sim$rasterToMatchLarge,
        overwrite = TRUE,
        url = extractURL("firePerimeters"),
        fireField = "YEAR",
        omitArgs = "destinationPath",
        userTags = c(cacheTags, "firePerimeters")
      )
    }
  }

  ## Stand age map ------------------------------------------------
  if (!suppliedElsewhere("standAgeMap", sim)) {
    if (P(sim)$dataYear == 2001) {
      ageURL <- extractURL("standAgeMap")
    } else if (P(sim)$dataYear == 2011) {
      ageURL <- extractURL("standAgeMap") |> gsub("2001", "2011", x = _)
    } else {
      stop("'P(sim)$dataYear' must be 2001 OR 2011")
    }
    ## Ceres Sep 3rd 2022 -- this option caused failure when previously set to FALSE at project level.
    # opt <- options("reproducible.useTerra" = TRUE) # Too many times this was failing with non-Terra # Eliot March 8, 2022
    # on.exit(options(opt), add = TRUE)
    sa <- if (is(sim$studyAreaLarge, "sf")) {
      aggregate(sim$studyAreaLarge, list(rep(1, nrow(sim$studyAreaLarge))),
                FUN = function(x) x)
    } else {
      aggregate(sim$studyAreaLarge)
    }
    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      sim$standAgeMap <- Cache(LandR::prepInputsStandAgeMap,
                               ageFun = getOption("reproducible.rasterRead", "terra::rast"), # the backwards compatible default
                               destinationPath = dPath,
                               ageURL = ageURL,
                               studyArea = sa,
                               rasterToMatch = sim$rasterToMatchLarge,
                               # writeTo = .suffix("standAgeMap.tif", paste0("_", P(sim)$.studyAreaName)),
                               overwrite = TRUE,
                               useCache = FALSE, ### for now due to attributes being lost on retrieval
                               firePerimeters = if (P(sim)$overrideAgeInFires) sim$firePerimeters else NULL,
                               fireURL = if (P(sim)$overrideAgeInFires) extractURL("firePerimeters") else NULL,
                               startTime = start(sim),
                               userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
                               omitArgs = c("destinationPath", "targetFile", "overwrite",
                                            "alsoExtract", "userTags"))
    })
    # options(opt)
  }

  LandR::assertStandAgeMapAttr(sim$standAgeMap)
  sim$imputedPixID <- attr(sim$standAgeMap, "imputedPixID")

  ## check parameter consistency across modules
  paramCheckOtherMods(sim, "dataYear", ifSetButDifferent = "warning")
  paramCheckOtherMods(sim, "minCoverThreshold", ifSetButDifferent = "warning")

  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")
  paramCheckOtherMods(sim, "vegLeadingProportion", ifSetButDifferent = "error")

  ## Species equivalencies table and associated columns ----------------------------
  ## make sppEquiv table and associated columns, vectors
  ## do not use suppliedElsewhere here as we need the tables to exist (or not)
  ## already (rather than potentially being supplied by a downstream module)
  ## the function checks whether the tables exist internally.

  sppOuts <- sppHarmonize(sim$sppEquiv, sim$sppNameVector, P(sim)$sppEquivCol,
                          sim$sppColorVect, P(sim)$vegLeadingProportion, sim$studyAreaLarge)
  ## the following may, or may not change inputs
  sim$sppEquiv <- sppOuts$sppEquiv
  sim$sppNameVector <- sppOuts$sppNameVector
  P(sim, module = currentModule(sim))$sppEquivCol <- sppOuts$sppEquivCol
  sim$sppColorVect <- sppOuts$sppColorVect

  ## check again
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")

  ## Species raster layers -------------------------------------------
  if (!suppliedElsewhere("speciesLayers", sim)) {
    #opt <- options(reproducible.useCache = "overwrite")
    # opt <- options("reproducible.useTerra" = TRUE) # Too many times this was failing with non-Terra # Eliot March 8, 2022
    # on.exit(options(opt), add = TRUE)
    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      sim$speciesLayers <- Cache(prepSpeciesLayers_KNN,
                                 destinationPath = dPath, # this is generic files (preProcess)
                                 outputPath = dPath,
                                 studyArea = sim$studyAreaLarge,
                                 studyAreaName = P(sim)$.studyAreaName,
                                 rasterToMatch = sim$rasterToMatchLarge,
                                 sppEquiv = sim$sppEquiv,
                                 sppEquivCol = P(sim)$sppEquivCol,
                                 thresh = 10,
                                 year = P(sim)$dataYear,
                                 userTags = c(cacheTags, "speciesLayers"),
                                 omitArgs = c("userTags"))
    })
    # options(opt)

    ## make sure empty pixels inside study area have 0 cover, instead of NAs.
    ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
    sim$speciesLayers <- NAcover2zero(sim$speciesLayers, sim$rasterToMatchLarge)
  }

  # 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath = dPath, cacheTags = c(cacheTags, "speciesTable"))
  }

  if (!suppliedElsewhere("columnsForPixelGroups", sim)) {
    sim$columnsForPixelGroups <- LandR::columnsForPixelGroups()
  }

  return(invisible(sim))
}

