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
  version = list(Biomass_borealDataPrep = "1.5.9.9000"),
  timeframe = as.POSIXlt(c(NA, NA)),
  timeunit = "year",
  citation = list("citation.bib"),
  documentation = list("README.txt", "Biomass_borealDataPrep.Rmd"),
  loadOrder = list(after = c("Biomass_speciesData"),
                   before = c("Biomass_core")),
  reqdPkgs = list(
    "archive", "assertthat", "crayon", "data.table", "dplyr", "fasterize",  "ggplot2", "httr2",
    "merTools", "plyr", "rasterVis", "sf", "terra",
    "reproducible (>= 2.1.0)", "SpaDES.core (>= 2.1.0)", "SpaDES.tools (>= 2.0.0)",
    "PredictiveEcology/LandR@development (>= 1.1.5.9058)",
    "PredictiveEcology/SpaDES.project@development (>= 0.0.8.9026)",
    "PredictiveEcology/pemisc@development"
  ),
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
    defineParameter("earliestFireYear", "integer", 1950L, NA, NA,
                    paste("if using fires to impute stand age and biomass, the earliest year for which",
                          "fire data should be obtained")),
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
                    paste(
                      "Model to estimate the relationship between % cover and % biomass, referred to as",
                      "`P(sim)$fitDeciduousCoverDiscount`. It is a number between 0 and 1 that translates % cover,",
                      "as provided in several databases, to % biomass. It is assumed that all hardwoods",
                      "are equivalent and all softwoods are equivalent and that % cover of hardwoods will",
                      "be an overestimate of the % biomass of hardwoods. E.g., 30% cover of hardwoods",
                      "might translate to 20% biomass of hardwoods. The reason this discount exists is",
                      "because hardwoods in Canada have a much wider canopy than softwoods."
                    )),
    defineParameter("deciduousCoverDiscount", "numeric", 0.8418911, NA, NA,
                    paste("This was estimated with data from NWT on March 18, 2020 and may or may not be universal.",
                          "Will not be used if `P(sim)$fitDeciduousCoverDiscount == TRUE`")),
    defineParameter("fitDeciduousCoverDiscount", "logical", FALSE, NA, NA,
                    paste("If TRUE, this will re-estimate `P(sim)$fitDeciduousCoverDiscount` This may be unstable and",
                          "is not recommended currently. If `FALSE`, will use the current default")),
    ## -------------------------------------------------------------------------------------------
    defineParameter("dataSource", "character", "SCANFI", NA, NA,
                    paste(
                      "Source for species cover, biomass, age, and landcover data used to initialize cohorts.",
                      "Currently, only kNN (2001, 2011) and SCANFI (2020) provide all necesarry layers.",
                      "Mixing multiple datasets requires additonal raster geoprocessing and is not recommended."
                    )),
    defineParameter("dataYear", "numeric", 2020, NA, NA,
                    paste(
                      "the year for which SCANFI data wil be fetched for use with the module.",
                      "One of 2000, 2010, or 2020, but note that only 2020 is currently supported." ## TODO
                    )),
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
                    paste("The classes in the `rstLCC` layer that are 'treed' and will therefore be run in `Biomass_core`.",
                          "Defaults to forested classes in NTEMS map (210 conif, 220 deciduous, 230 mixed) plus",
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
                      "and applies them to all ecolocations (`ecoregionGroup` codes)."
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
                    paste("When assigning `pixelGroup` membership, this defines the resolution of biomass that will be considered",
                          "'the same `pixelGroup`', e.g., if it is 100, then 5160 and 5240 will be the same")),
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
    defineParameter("sppEquivCol", "character", "LandR", NA, NA,
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
                    "a number that defines whether a species is leading for a given pixel"),
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
                    paste(
                      "Named list of seeds to use for each event (names). E.g., `list('init' = 123)` will `set.seed(123)`",
                      "at the start of the init event and unset it at the end. Defaults to `NULL`, meaning that",
                      "no seeds will be set"
                    )),
    defineParameter(".sslVerify", "integer", as.integer(unname(curl::curl_options("^ssl_verifypeer$"))), NA_integer_, NA_integer_,
                    paste(
                      "Passed to `httr::config(ssl_verifypeer = P(sim)$.sslVerify)` when downloading NFI datasets.",
                      "Set to 0L if necessary to bypass checking the SSL certificate",
                      "(this may be necessary when NFI's website SSL certificate is not correctly configured)."
                    )),
    defineParameter(".studyAreaName", "character", NA, NA, NA,
                    "Human-readable name for the study area used. If `NA`, a hash of studyArea will be used."),
    defineParameter(".useCache", "character", c(".inputObjects", "init"), NA, NA,
                    "Internal. Can be names of events or the whole module name; these will be cached by SpaDES")
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
                              "It must have same extent and crs as `studyArea_biomassParam`.",
                              "It is superseded by `sim$ecoregionRst` if that object is supplied by the user"),
                 sourceURL = "https://sis.agr.gc.ca/cansis/nsdb/ecostrat/district/ecodistrict_shp.zip"),
    expectsInput("ecoregionRst", "SpatRaster",
                 desc = paste("A raster that characterizes the unique ecological regions used to",
                              "parameterize the biomass, cover, and species establishment probability models.",
                              "If this object is provided, it will supercede `sim$ecoregionLayer`.",
                              "It will be overlaid with landcover to generate classes for every ecoregion/LCC combination.",
                              "It must have same extent and crs as `rasterToMatch_biomassParam` if supplied by user - use `reproducible::postProcess`.",
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
                 paste(
                   "A land classification map in study area. It must be 'corrected', in the sense that:\n",
                   "1) Every class must not conflict with any other map in this module\n",
                   "    (e.g., `speciesLayers` should not have data in LCC classes that are non-treed);\n",
                   "2) It can have treed and non-treed classes. The non-treed will be removed within this\n",
                   "    module if `P(sim)$omitNonTreedPixels` is `TRUE`;\n",
                   "3) It can have transient pixels, such as 'young fire'. These will be converted to a\n",
                   "    the nearest non-transient class, probabilistically if there is more than 1 nearest\n",
                   "    neighbour class, based on `P(sim)$LCCClassesToReplaceNN`.\n",
                   "The default layer used, if not supplied, is SCANFI-derived data product for 2020 updated to use NTEMS land cover codes.",
                   "See <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc> for SCANFI metadata.",
                   "The metadata (res, proj, ext, origin) need to match `rasterToMatch_biomassParam`."),
                 sourceURL = NA), ## uses P(sim)$rstLCCYear and LandR::prepInputsLCC() defaults
    expectsInput("rasterToMatch", "SpatRaster",
                 desc = paste("A raster of the `studyArea` in the same resolution and projection as `rawBiomassMap`.",
                              "This is the scale used for all *outputs* for use in the simulation.",
                              "If not supplied will be forced to match the *default* `rawBiomassMap`.")),
    expectsInput("rasterToMatch_biomassParam", "SpatRaster",
                 desc = paste("A raster of the `studyArea_biomassParam` in the same resolution and projection as `rawBiomassMap`.",
                              "This is the scale used for all *inputs* for use in the simulation.",
                              "If not supplied will be forced to match the *default* `rawBiomassMap`.")),
    expectsInput("rawBiomassMap", "SpatRaster",
                 paste(
                   "total biomass raster layer in study area. Defaults to the Canadian Forestry",
                   "Service, National Forest Inventory, SCANFI-derived total aboveground biomass map",
                   "from 2020 (in tonnes/ha), unless `dataYear != 2020`.",
                   "See <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc> for metadata."
                 )),
    expectsInput("speciesLayers", "SpatRaster",
                 paste(
                   "cover percentage raster layers by species in Canada species map.",
                   "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                   "SCANFI-derived species cover maps from 2020 using a cover threshold of 10 -",
                   "see <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc> for metadata"
                 )),
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
                 paste(
                   "stand age map in study area. Must have a 'imputedPixID' attribute (a  vector of pixel IDs)",
                   "indicating which pixels suffered age imputation. If no pixel ages were imputed, please set",
                   "this attribute to `integer(0)`.",
                   "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                   "SCANFI-derived biomass map from 2020, unless `dataYear != 2020`.",
                   "See <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc> for metadata."
                  )),
    expectsInput("studyArea", "sf",
                 desc = paste("`sf` polygon or terra `SpatVector` to use as the study area - `nrow` must be one")),
    expectsInput("studyArea_biomassParam", "sf",
                 desc = paste("Polygon to use as the parametrisation study area. Must be provided by the user.",
                              "Note that `studyArea_biomassParam` is only used for parameter estimation, and",
                              "can be larger than the actual study area used for LandR simulations",
                              "(e.g., larger than `studyArea` in LandR `Biomass_core`)."))
  ),
  outputObjects = bindrows(
    createsOutput("biomassMap", "SpatRaster",
                  paste("total biomass raster layer in study area,",
                        "filtered for pixels covered by `cohortData`. Units in $g/m^2$")),
    createsOutput("cohortData", "data.table",
                  paste("initial community table, containing corrected biomass ($g/m^2$), age and",
                        "species cover data, as well as ecolocation and `pixelGroup` information. This table defines",
                        "the initial community composition and structure used by `Biomass_core`")),
    createsOutput("ecoregion", "data.table",
                  paste("`ecoregionGroup` look up table")),
    createsOutput("ecoregionMap", "SpatRaster",
                  paste("`ecoregionGroup` map that has mapcodes match `ecoregion` table and `speciesEcoregion` table")),
    createsOutput("firePerimeters", "SpatRaster",
                  paste("As the input object `firePerimeters`, but potentially cropped/masked/projected to match `rasterToMatch_biomassParam`")),
    createsOutput("imputedPixID", "integer",
                  paste("A vector of pixel IDs - matching `rasterMatch` IDs - that suffered data imputation.",
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
    #                     "Service, National Forest Inventory, SCANFI-derived total aboveground biomass map",
    #                     "(in tonnes/ha) from 2020, unless `dataYear != 2020`.",
    #                     "See <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc>",
    #                     "for metadata")),
    createsOutput("rstLCC", "SpatRaster",
                  paste("As the input object `rstLCC`, but potentially cropped/projected/masked",
                        "to match `rasterToMatch_biomassParam`")),
    createsOutput("species", "data.table",
                  paste("Table that of invariant species traits.",
                        "Will have the same traits as the input `speciesTable`,",
                        "with values adjusted where necessary.")),
    createsOutput("speciesLayers", "SpatRaster",
                  paste("cover percentage raster layers by species in Canada species map.",
                        "Defaults to the Canadian Forestry Service, National Forest Inventory,",
                        "SCANFI-derived species cover maps from 2020 using a cover threshold of 10 -",
                        "see <https://open.canada.ca/data/en/dataset/18e6a919-53fd-41ce-b4e2-44a9707c52dc>",
                        "for metadata.")),
    createsOutput("speciesEcoregion", "data.table",
                  paste("table of spatially-varying species traits (`maxB`, `maxANPP`, `establishprob`),",
                        "defined by species and `ecoregionGroup` (i.e. ecolocation)")),
    createsOutput("standAgeMap", "SpatRaster",
                  paste("As the input object `standAgeMap`, but potentially cropped, projected,",
                        "masked to match `rasterToMatch_biomassParam`.")),
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
  if (!.compareRas(sim$standAgeMap, sim$rasterToMatch_biomassParam, res = TRUE)) {
    ## note that extents may never align if the resolution and projection do not allow for it
    ## this is not working, need to use projectRaster
    sim$standAgeMap <- postProcess(
      sim$standAgeMap,
      to = sim$rasterToMatch_biomassParam,
      overwrite = TRUE
    ) |>
      Cache(.functionName = "postProcessStandAgeMap")
    attr(sim$standAgeMap, "imputedPixID") <- sim$imputedPixID
  }

  if (!.compareRas(sim$rstLCC, sim$rasterToMatch_biomassParam, res = TRUE)) {
    sim$rstLCC <- postProcess(sim$rstLCC,
                              to = sim$rasterToMatch_biomassParam,
                              overwrite = TRUE) |> Cache(.functionName = "postProcessRstLCC")
  }

  if (P(sim)$overrideAgeInFires) {
    sim$firePerimeters <- postProcess(sim$firePerimeters,
                                      to = sim$rasterToMatch_biomassParam,
                                      overwrite = TRUE) |> Cache(.functionName = "postProcessFirePerimeters")
  }
  # options(opt)
  if (!.compareRas(sim$speciesLayers, sim$rasterToMatch_biomassParam, res = TRUE)) {
    sim$speciesLayers <- postProcessTerra(sim$speciesLayers,
                                          to = sim$rasterToMatch_biomassParam,
                                          overwrite = TRUE) |> Cache(.functionName = "postProcessSpeciesLayers")
  }

  if (!.compareRas(sim$rasterToMatch_biomassParam, sim$rawBiomassMap, sim$rstLCC,
                   sim$speciesLayers, sim$standAgeMap, res = TRUE)) {
    stop(paste("sim$rasterToMatch_biomassParam, sim$rawBiomassMap, sim$rstLCC",
               "sim$speciesLayers, sim$standAgeMap properties do not match"))
  }

  ## species traits inputs ---------------------------------------
  message(blue("Prepare 'species' table, i.e., species level traits", Sys.time()))

  sim$species <- prepSpeciesTable(speciesTable = sim$speciesTable,
                                  sppEquiv = sim$sppEquiv,
                                  areas = P(sim)$speciesTableAreas,
                                  sppEquivCol = P(sim)$sppEquivCol) |> Cache()

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
  sim$sufficientLight <- data.frame(
    speciesshadetolerance = 1:5,
    X0 = c(rep(1, 4), 0),
    X1 = c(0, rep(1, 3), 0),
    X2 = c(0, 0, rep(1, 3)),
    X3 = c(rep(0, 3), rep(1, 2)),
    X4 = c(rep(0, 4), 1),
    X5 = c(rep(0, 4), 1)
  )

  ## initialEcoregionMap -----------------------------------------
  if (!.compareCRS(sim$studyArea, sim$rasterToMatch)) {
    warning("studyArea and rasterToMatch projections differ")
  }

  ## Clean pixels for veg. succession model
  ## remove pixels with no species data or non-forested LCC
  ## ELIOT: updated May 7, 2025. The way pixelFateDT works is that it separates the NAs
  ##        from the nonForestedPixels. This next line (nonForestedPixels) does both, confounded.
  ##        Need to separate these steps (NA in speciesLayers and nonForest classes)
  ##        Now it is run twice here, then after the printed line about NAs removed
  pixelsToRmDueToNAs <- nonForestedPixels(sim$speciesLayers, omitNonTreedPixels = FALSE)
  pixelFateDT <- pixelFate(fate = "Total number pixels", runningPixelTotal = ncell(sim$speciesLayers))
  pixelFateDT <- pixelFate(pixelFateDT, "NAs on sim$speciesLayers", sum(pixelsToRmDueToNAs))
  if (P(sim)$omitNonTreedPixels) {
    checkNonforest <- sum(!(as.vector(sim$rstLCC[]) %in% P(sim)$forestedLCCClasses)) -
      tail(pixelFateDT$pixelsRemoved, 1)
    if (checkNonforest < 0) browser() ## TODO: remove browser
    pixelFateDT <- pixelFate(pixelFateDT, "Non forested pixels (based on LCC classes)", checkNonforest)
  }
  pixelsToRmDueToNAsAndNonForest <- nonForestedPixels(
    sim$speciesLayers,
    P(sim)$omitNonTreedPixels,
    P(sim)$forestedLCCClasses,
    sim$rstLCC
  )

  ## The next function will remove the "zero" class on sim$ecoregionRst
  pixelFateDT <- pixelFate(pixelFateDT, "Removing 0 class in sim$ecoregionRst",
                           sum(as.vector(sim$ecoregionRst[])[!pixelsToRmDueToNAsAndNonForest] == 0, na.rm = TRUE))

  rstLCCAdj <- sim$rstLCC
  rstLCCAdj[pixelsToRmDueToNAsAndNonForest] <- NA

  ## make initial ecoregionFiles - some of these may have LCC that get replaced
  ecoregionFiles <- prepEcoregions(
    ecoregionRst = sim$ecoregionRst,
    ecoregionLayer = sim$ecoregionLayer,
    ecoregionLayerField = P(sim)$ecoregionLayerField,
    rasterToMatchLarge = sim$rasterToMatch_biomassParam,
    rstLCCAdj = rstLCCAdj,
    pixelsToRm = pixelsToRmDueToNAsAndNonForest,
    cacheTags = c(cacheTags, "prepEcoregionFiles")
  ) |>
    Cache()

  ## create pixelTable object ------------------------------------
  ##  Round age to pixelGroupAgeClass
  ##  Internal data.table is changed; using memoise here causes the internal changes to
  ##  come out to the pixelTable, which is not desired. Turn off memoising for one step
  opt <- options("reproducible.useMemoise" = FALSE)
  on.exit(try(options(opt), silent = TRUE), add = TRUE)

  pixelTable <- makePixelTable(
    speciesLayers = sim$speciesLayers,
    standAgeMap = sim$standAgeMap,
    ecoregionFiles = ecoregionFiles,
    biomassMap = sim$rawBiomassMap,
    rasterToMatch = sim$rasterToMatch_biomassParam,
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
    Cache(userTags = c(cacheTags, "pixelCohortData"))
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

  # Cache here, uses the previously digested object that was used to create the pixelCohortData; it hasn't
  #   changed in the code above since its creation just above
  pixelCohortData <- partitionBiomass(x = P(sim)$deciduousCoverDiscount, pixelCohortData) |>
    Cache(omitArgs = "pixelCohortData", .cacheExtra = attr(pixelCohortData, "tags"))

  set(pixelCohortData, NULL, "B", asInteger(pixelCohortData$B/P(sim)$pixelGroupBiomassClass) *
        P(sim)$pixelGroupBiomassClass)
  set(pixelCohortData, NULL, "cover", asInteger(pixelCohortData$cover))

  ## replace unwanted LCC classes to a neighbour class *that exists*.------------------------------------
  ## Originally 34/36 (hence the name) values from 2005 LCC, which were burns and cities.
  ## We need to have a spatial estimate of maxBiomass everywhere there is forest; can't have gaps.
  ## The LCC that are unwanted are places for which we don't want
  ## maxBiomass associated with their LCC ... i.e., we don't want a maximum
  ## biomass associated with disturbed forest because those classes are transient.
  ## They will transition to another class before they arrive at a tree maximum biomass.
  ## However, we need to give them a "forest class" that they might "become"
  ## The ecoregion map must be updated to reflect this new class.

  NTEMSlcc <- c(0, 20, 30, 31, 32, 33, 40, 50, 80, 81, 100, 210, 220, 230, P(sim)$LCCClassesToReplaceNN)
  ## unclassified, water, snow/ice, rock/rubble, exposed/barren,
  ## bryoids, shrubs, wetland, wetland-treed, herbs, coniferous,
  ## broadleaf, mixedwood, disturbed)
  if (length(P(sim)$LCCClassesToReplaceNN) && all(na.omit(as.vector(sim$rstLCC)) %in% NTEMSlcc)) {
    uwc <- P(sim)$LCCClassesToReplaceNN
    message("Replace ", paste(uwc, collapse = ", "), " values to a neighbour class *that exists*")
    availableCombinations <- unique(pixelCohortData[, .(speciesCode, initialEcoregionCode, pixelIndex)])

    freqsUpdates <- startFinishLCC <- list()
    lastYrOnNTEMS <- NTEMSfinalYearForLCC(timeout = 10) |> Cache()

    for (yr in c(lastYrOnNTEMS, 1984)) {
      freqs <- freq(rstLCCAdj)
      num2replace <- freqs$count[freqs$value %in% P(sim)$LCCClassesToReplaceNN]
      if ((length(num2replace) > 0) && (num2replace > 1000)) {
        yrChar <- as.character(yr)
        startFinishLCC[[yrChar]] <-
          prepInputs_NTEMS_LCC_FAO(
            year = yr,
            to = sim$rstLCC,
            disturbedCode = 240,
            destinationPath = inputPath(sim),
            overwrite = TRUE
          ) |>
          Cache(userTags = c("rstLCC", yr, "_", currentModule(sim),
                             P(sim)$.studyAreaName, P(sim)$dataYear))
        pixelTable <- copy(pixelTable) ## avoid super annoying warning
        cellsToUpdate <- which(rstLCCAdj[] == P(sim)$LCCClassesToReplaceNN)
        rstLCCAdj[cellsToUpdate] <- startFinishLCC[[yrChar]][cellsToUpdate]
        whUpdate <- match(cellsToUpdate, pixelTable$pixelIndex)
        pixelTable[whUpdate, newLcc := startFinishLCC[[yrChar]][cellsToUpdate]]
        pixelTable[whUpdate, initialEcoregionCode2 := gsub("_.+", "", initialEcoregionCode)]
        pixelsToRm2 <- nonForestedPixels(rstLCCAdj, P(sim)$omitNonTreedPixels,
                                         P(sim)$forestedLCCClasses, rstLCCAdj)
        pixelsToRm3 <- which(pixelsToRm2)
        pixelsToRm4 <- na.omit(match(pixelsToRm3, pixelTable$pixelIndex))
        if (length(pixelsToRm4)) {
          ## are there any new ones that are not forestedLCCClasses?
          pixelTable <- pixelTable[-pixelsToRm4]
        }

        ncharToPad <- max(nchar(pixelTable$lcc))

        ## Eliot added this after many failed assertions WAY below: Sep 5, 2025
        ##   assert_that(all(is.na(values(mat = FALSE, sim$ecoregionMap)) == is.na(values(mat = FALSE, sim$pixelGroupMap))))
        ##   The newLcc
        pixelTable <- pixelTable[!newLcc %in% 0] # These are pixels that turned to zero i.e., need to be removed

        pixelTable[!is.na(newLcc), lcc := newLcc]
        pixelTable[!is.na(newLcc),  initialEcoregionCode :=
                     paste0(initialEcoregionCode2, "_",
                            paddedFloatToChar(newLcc, ncharToPad))]

        set(pixelTable, NULL, c("newLcc", "initialEcoregionCode2"), NULL)
        rstLCCAdj[pixelsToRm2] <- NA
        rm(pixelsToRm2, pixelsToRm3, pixelsToRm4)
      }
    }

    ## create initial pixelCohortData table ---------------
    ## Might already have cover. in the names
    coverColNames <- colnames(pixelTable)[match(sim$species$species,
                                                 gsub("cover.(.+)", "\\1", colnames(pixelTable)))]
    # coverColNames <- paste0("cover.", coverColNames)

    pixelCohortData <- makeAndCleanInitialCohortData(
      inputDataTable = pixelTable,
      sppColumns = coverColNames,
      imputeBadAgeModel = P(sim)$imputeBadAgeModel,
      minCoverThreshold = P(sim)$minCoverThreshold,
      doSubset = P(sim)$subsetDataAgeModel
    ) |>
      Cache(userTags = c(cacheTags, "pixelCohortData"), omitArgs = c("userTags"))
    assertCohortDataAttr(pixelCohortData)
    pixelCohortData <- Cache(partitionBiomass(x = P(sim)$deciduousCoverDiscount, pixelCohortData))
    set(pixelCohortData, NULL, "B", asInteger(pixelCohortData$B / P(sim)$pixelGroupBiomassClass) *
          P(sim)$pixelGroupBiomassClass)
    set(pixelCohortData, NULL, "cover", asInteger(pixelCohortData$cover))
    pixelCohortData <- pixelCohortData[!is.na(pixelCohortData$lcc)]
    availableCombinations2 <- unique(pixelCohortData[, .(speciesCode, initialEcoregionCode, pixelIndex)])

    newLCCClasses <- convertUnwantedLCC(
      classesToReplace = P(sim)$LCCClassesToReplaceNN,
      rstLCC = rstLCCAdj,
      availableERC_by_Sp = availableCombinations2
    ) |>
      Cache(userTags = c(cacheTags, "newLCCClasses", "stable"))

    ## adjust rstLCCAdj so that ecoregionMap will contain the last set of updated LCCClassesToReplaceNN
    if (nrow(newLCCClasses)) {
      if (!is.null(newLCCClasses$newPossLCC)) {
        ## LandR versions prior to 1.1.5.9045 will not have this
        rstLCCAdj[newLCCClasses$pixelIndex] <- newLCCClasses$newPossLCC
      }
    }
  } else {
    newLCCClasses <- data.table(pixelIndex = numeric(), ecoregionGroup = numeric())
  }

  sim$imputedPixID <- unique(c(sim$imputedPixID, newLCCClasses$pixelIndex))
  ## split pixelCohortData into 2 parts -- one with the former P(sim)$LCCClassesToReplaceNN pixels, one without
  ##   The one without P(sim)$LCCClassesToReplaceNN can be used for statistical estimation, but not the one with
  cohortDataOnlyNonForestLCC <- pixelCohortData[pixelIndex %in% newLCCClasses$pixelIndex]
  cohortDataOnlyNonForestLCC <- merge(newLCCClasses, cohortDataOnlyNonForestLCC, all.x = TRUE,
                                      all.y = FALSE, by = "pixelIndex")
  cohortDataOnlyForestLCC <- pixelCohortData[!pixelIndex %in% newLCCClasses$pixelIndex]
  if (!length(P(sim)$LCCClassesToReplaceNN)) {
    if (!all.equal(cohortDataOnlyForestLCC, pixelCohortData, check.attributes = FALSE))
      stop("No LCC classes were listed for replacement, but some pixels may have been lost")
  }
  setnames(cohortDataOnlyForestLCC, "initialEcoregionCode", "ecoregionGroup")
  rmZeroBiomassQuote <- quote(totalBiomass > 0)
  cohortDataOnlyForestLCCBiomass <- cohortDataOnlyForestLCC[eval(rmZeroBiomassQuote),
                                                  .(B, logAge, speciesCode, ecoregionGroup, lcc, cover)]
  cohortDataOnlyForestLCCBiomass <- unique(cohortDataOnlyForestLCCBiomass)

  ## make sure ecoregionGroups match
  ## remember to match rmZeroBiomassQuote the rule used to filter `availableCombinations` (NULL if none)
  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert1(cohortDataOnlyNonForestLCC, pixelCohortData, rmZeroBiomassQuote = NULL,
            classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(cohortDataOnlyForestLCC, classesToReplace = P(sim)$LCCClassesToReplaceNN)
  }

  ## Statistical estimation of establishprob, maxB and maxANPP ----------------------
  cohortDataShort <- cohortDataOnlyForestLCC[, list(coverPres = sum(cover > 0)),
                                        by = c("ecoregionGroup", "speciesCode")]
  ## find coverNum for each known class
  ## add new ecoregions to pixelTable, before calc. table
  cohortDataShortNoCover <-
    (function(x) {
      tempDT <- rbind(cohortDataOnlyNonForestLCC[, .(pixelIndex, ecoregionGroup)],
                      cohortDataOnlyForestLCC[, .(pixelIndex, ecoregionGroup)])
      pixelTable <- tempDT[pixelTable, on = .(pixelIndex)]

      aa <- table(as.character(pixelTable$ecoregionGroup)) ## as.character avoids counting levels that don't exist anymore

      dt1 <- data.table(ecoregionGroup = factor(names(aa)), coverNum = as.integer(unname(aa)))
      allCombos <- expand.grid(ecoregionGroup = dt1$ecoregionGroup, speciesCode = unique(cohortDataShort$speciesCode))
      setDT(allCombos)
      dt1 <- dt1[allCombos, on = "ecoregionGroup", nomatch = 0]
      cohortDataShortNoCover <- cohortDataShort[dt1, on = c("ecoregionGroup", "speciesCode"), nomatch = NA]
    })() |>
    Cache(
      .functionName = "cohortDataShortNoCover",
      .cacheExtra = list(
        a = cohortDataOnlyNonForestLCC[, .(pixelIndex, ecoregionGroup)],
        b = cohortDataOnlyForestLCC[, .(pixelIndex, ecoregionGroup)],
        d = pixelTable,
        e = cohortDataShort
      )
    )

  cohortDataShort <- cohortDataShortNoCover[coverPres > 0] ## remove places where there is 0 cover
  cohortDataShortNoCover <- cohortDataShortNoCover[is.na(coverPres)][, coverPres := 0]
  ##  will be added back as establishprob = 0

  if (length(P(sim)$LCCClassesToReplaceNN)) {
    assert2(cohortDataShort, classesToReplace = P(sim)$LCCClassesToReplaceNN)
    assert2(cohortDataShortNoCover, classesToReplace = P(sim)$LCCClassesToReplaceNN)

    ## rebuild ecoregionFiles with updated rstLCC
    ecoregionFiles <- prepEcoregions(
      ecoregionRst = sim$ecoregionRst,
      ecoregionLayer = sim$ecoregionLayer,
      ecoregionLayerField = P(sim)$ecoregionLayerField,
      rasterToMatchLarge = sim$rasterToMatch_biomassParam,
      rstLCCAdj = rstLCCAdj,
      pixelsToRm = pixelsToRmDueToNAsAndNonForest,
      cacheTags = c(cacheTags, "prepEcoregionFiles")
    ) |>
      Cache()
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
    .cacheExtra = levels(cohortDataShort$speciesCode), ## in case sppEquivCol changes
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
  totalBiomass <- sum(cohortDataOnlyForestLCCBiomass$B, na.rm = TRUE)

  ## There are several reasons why the modelBiomass can fail;
  ##   1) inappropriate sub-sample
  ##   2) fit algorithm
  ## Run two nested loops to do both of these things
  modelBiomassTags <- c(cacheTags, "modelBiomass",
                        paste0("subsetSize:", P(sim)$subsetDataBiomassModel))
  maxDataSubsetTries <- ifelse(isTRUE(P(sim)$subsetDataBiomassModel > 0),
                               P(sim)$subsetDataAttempts, 1)
  for (tryBiomassDataSubset in 1:maxDataSubsetTries) {
    cohortDataOnlyForestLCCBiomassSubset <- subsetDT(cohortDataOnlyForestLCCBiomass,
                                                     by = c("ecoregionGroup", "speciesCode"),
                                                     doSubset = P(sim)$subsetDataBiomassModel)

    ## For Cache: doesn't need to cache all columns in the data.table; only the ones in the model.
    ## force parameter values to avoid more checks;
    ## If using mixed effect model, see here for good discussion of
    ##  shrinkage https://www.tjmahr.com/plotting-partial-pooling-in-mixed-effects-models/
    message(blue("Estimating biomass using P(sim)$biomassModel as:"), "\n",
            magenta(paste0(format(P(sim)$biomassModel, appendLF = FALSE), collapse = "")))

    ## NOTE: we are NOT using logB because the relationship between B~age should be hump-shaped
    ## (or at least capped at high age values). Ideally, we would want a non-linear model

    ## Default values of args to modelBiomass -- prior to any attempts to fix
    ueg <- .sortDotsUnderscoreFirst(as.character(unique(cohortDataOnlyForestLCCBiomassSubset$ecoregionGroup)))
    specDat <- cohortDataOnlyForestLCCBiomassSubset
    modelFn <- P(sim)$biomassModel
    sumResponse <- c(totalBiomass)

    fixModelBiomass <- P(sim)$fixModelBiomass
    timePriorToFit <- Sys.time()
    cohortDataOnlyForestLCCBiomassSubset2 <- copy(cohortDataOnlyForestLCCBiomassSubset)

    tryControl <- FALSE
    needRescaleModelB <- FALSE
    scaledVarsModelB <- NULL
    for (tryBiomassModel in 1:3) { ## try thrice -- default, then once to rescale, once to refit
      modelBiomass <- Cache(
        statsModel,
        modelFn = modelFn,
        uniqueEcoregionGroups = ueg,
        .cacheExtra = sumResponse, ## only digest on this
        .specialData = specDat,
        useCloud = useCloud,
        cloudFolderID = sim$cloudFolderID,
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
          logAge_sc <- scale(cohortDataOnlyForestLCCBiomassSubset$logAge)
          cover_sc <- scale(cohortDataOnlyForestLCCBiomassSubset$cover)

          scaledVarsModelB <- list(logAge = logAge_sc, cover = cover_sc)
          ## remove attributes with as.numeric
          ## don't change the original data
          cohortDataOnlyForestLCCBiomassSubset2[, `:=`(logAge = as.numeric(logAge_sc),
                                                  cover = as.numeric(cover_sc))]
          needRescaleModelB <- TRUE
          ueg <- .sortDotsUnderscoreFirst(as.character(unique(cohortDataOnlyForestLCCBiomassSubset2$ecoregionGroup)))
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
        suppressMessages(clearCache(userTags = userTagsToClear, # after = timePriorToFit,
                                    ask = FALSE))
        specDat <- cohortDataOnlyForestLCCBiomassSubset2
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

  if (any(P(sim)$exportModels %in% c("all", "biomassModel"))) {
    sim$modelBiomass <- modelBiomass
  }

  ## create speciesEcoregion ---------------------------------------------
  ## a single line for each combination of ecoregionGroup & speciesCode;
  ## doesn't include combinations with B = 0 because those places can't have the species/ecoregion combo
  ## cohortDataOnlyForestLCCBiomassSubset ends up determining which ecoregion combinations end up in
  ## species ecoregion, thus removing converted/masked classes present cohortDataShortNoCover
  message(blue("Create speciesEcoregion using modelCover and modelBiomass to estimate species traits"))
  speciesEcoregion <- makeSpeciesEcoregion(cohortDataBiomass = cohortDataOnlyForestLCCBiomassSubset,
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

  if (ncell(sim$rasterToMatch_biomassParam) > 3e7) replicate(3, gc())

  ## Create initial communities, i.e., pixelGroups -----------------------
  ## Rejoin back the pixels that were P(sim)$LCCClassesToReplaceNN
  set(cohortDataOnlyNonForestLCC, NULL, "initialEcoregionCode", NULL)
  pixelCohortData <- rbindlist(list(cohortDataOnlyNonForestLCC, cohortDataOnlyForestLCC),
                               use.names = TRUE, fill = TRUE)

  ## "Downsize" to studyArea after estimating parameters on studyArea_biomassParam --------------
  ## 1. Subset pixels (IDs) on rasterToMatch_biomassParam, using rasterToMatch
  ## 2. Subset data.tables using the pixel IDs / ecoregion/species combinations
  ##    that are common across the two rasters
  ## 3. Re-do pixel ID numbering so that it matches the final rasterToMatch
  ## Note: if SA and SALarge are the same, no subsetting will take place.
  if (sum(is.na(as.vector(values(sim$rasterToMatch)))) != sum(is.na(as.vector(values(sim$rasterToMatch_biomassParam))))) {
    message(blue("Subsetting to studyArea"))
    rasterToMatch_biomassParam <- sim$rasterToMatch_biomassParam
    rasterToMatch_biomassParam <- setValues(rasterToMatch_biomassParam, seq(ncell(rasterToMatch_biomassParam)))

    opt <- options(reproducible.gdalwarp = FALSE) ## gdalwarp will reproject even if same CRS, duplicating indices
    on.exit(options(opt), add = TRUE)
    rasterToMatch_biomassParamCropped <- Cache(postProcess,
                                       x = rasterToMatch_biomassParam,
                                       to = sim$rasterToMatch,
                                       datatype = assessDataType(rasterToMatch_biomassParam),
                                       method = "near",
                                       userTags = c(cacheTags, "rasterToMatch_biomassParamCropped"),
                                       omitArgs = c("userTags"))
    options(opt)

    rtmlc_int <- LandR::asInt(rasterToMatch_biomassParamCropped)
    assertthat::assert_that(all(na.omit(as.vector(rasterToMatch_biomassParamCropped - rtmlc_int)) == 0))
    rm(rtmlc_int)
    assertthat::assert_that(sum(is.na(as.vector(rasterToMatch_biomassParamCropped))) < ncell(rasterToMatch_biomassParamCropped))
    ## i.e., not all NA

    if (!.compareRas(rasterToMatch_biomassParamCropped, sim$rasterToMatch)) {
      stop("Downsizing to rasterToMatch after estimating parameters didn't work.",
           "Please debug Biomass_borealDataPrep::createBiomass_coreInputs().")
    }

    ## subset pixels that are in studyArea/rasterToMatch only
    pixToKeep <- na.omit(as.vector(values(rasterToMatch_biomassParamCropped))) ## these are the old indices of RTML
    pixelCohortData <- pixelCohortData[pixelIndex %in% pixToKeep]

    ## re-do pixelIndex (it now needs to match rasterToMatch)
    newPixelIndexDT <- data.table(pixelIndex = as.vector(values(rasterToMatch_biomassParamCropped)),
                                  newPixelIndex = as.integer(1:ncell(rasterToMatch_biomassParamCropped))) |>
      na.omit()

    pixelCohortData <- newPixelIndexDT[pixelCohortData, on = "pixelIndex"]
    pixelCohortData[, pixelIndex := NULL]
    setnames(pixelCohortData, old = "newPixelIndex", new = "pixelIndex")

    assertthat::assert_that(NROW(pixelCohortData) > 0)

    ## now convert imputedPixID to RTM
    sim$imputedPixID <- newPixelIndexDT[pixelIndex %in% sim$imputedPixID, newPixelIndex]

    rm(pixToKeep, rasterToMatch_biomassParamCropped, newPixelIndexDT)
    if (ncell(sim$rasterToMatch) > 3e7) replicate(3, gc())
  }
  ## subset ecoregionFiles$ecoregionMap to smaller area.

  ecoregionFiles$ecoregionMap <- postProcess(
    x = ecoregionFiles$ecoregionMap,
    to = sim$rasterToMatch,
    writeTo = NULL
  ) |>
    Cache(
      .functionName = "postProcessEcoregionMap",
      userTags = c(cacheTags, "ecoregionMap"),
      omitArgs = c("userTags")
    )

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

      firstFireYear <- P(sim)$earliestFireYear
      ## this is not necessary when using min(),
      ## but will be kept in case we use something else in the future
      maxAgeHighQualityData <- P(sim)$dataYear - firstFireYear
      ## if maxAgeHighQualityData is lower than 0, it means it's prior to the first fire Year
      ## or not following calendar year

      if (isTRUE(maxAgeHighQualityData >= 0)) {
        ## identify young in the pixelCohortData
        youngRows <- pixelCohortData$age <= maxAgeHighQualityData
        young <- pixelCohortData[youngRows == TRUE]

        youngRows2 <- !is.na(sim$firePerimeters[][young$pixelIndex])
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
        ) ## /4 is too strong -- 25 years is a lot of time
      } else {
        ## return maxAgeHighQualityData to -1
        message(blue("Simulation start year is lower than oldest fire."))
        message(blue("B values will NOT be re-estimated inside fire perimeters"))
        maxAgeHighQualityData <- -1
      }
    }
  }

  assertthat::assert_that(all(inRange(na.omit(pixelCohortData$B), 0, round(maxRawB, -2)))) # should they all be below the initial biomass map?

  ## Fill in any remaining B values that are still NA -- the previous chunk filled in B for young cohorts only
  if (anyNA(pixelCohortData$B)) {

    theNAsBiomass <- is.na(pixelCohortData$B)
    message(blue(" -- ", sum(theNAsBiomass),"cohort(s) has NA for Biomass: being replaced with model-derived estimates"))
    set(pixelCohortData, which(theNAsBiomass), "B",
        pmax(0, asInteger(predict(modelBiomass$mod, newdata = pixelCohortData[theNAsBiomass],
                          allow.new.levels = TRUE))))
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
  #initialize with disturbed (i.e. empty) pixels as pixelGroup 0
  sim$pixelGroupMap[is.na(sim$pixelGroupMap[]) & !is.na(sim$ecoregionMap[])] <- 0 #
  assert_that(all(is.na(values(mat = FALSE, sim$ecoregionMap)) == is.na(values(mat = FALSE, sim$pixelGroupMap))))


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
  # out <- lapply(capture.output(sim$pixelFateDT), function(x) message(blue(x)))

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

  rtm_res <- 240 ## SCANFI is 30m resolution and would be aggregated to this

  ## Study area(s) ------------------------------------------------
  if (!suppliedElsewhere("studyArea", sim)) {
    ## Jan 2021 we agreed to force user to provide a SA/SAL
    sim$studyArea <- randomStudyArea(seed = 1234, size = (rtm_res^2)*100)
  }

  if (!suppliedElsewhere("studyArea_biomassParam", sim)) {
    if (is.null(sim$studyAreaLarge)) {
      sim$studyArea_biomassParam <- sim$studyArea
    } else {
      warning("please replace studyAreaLarge with studyArea_biomassParam")
      sim$studyArea_biomassParam <- sim$studyAreaLarge
    }
  }

  if (is.na(P(sim)$.studyAreaName)) {
    params(sim)[[currentModule(sim)]][[".studyAreaName"]] <- reproducible::studyAreaName(sim$studyArea_biomassParam)
    message("The .studyAreaName is not supplied; derived name from sim$studyArea_biomassParam: ",
            params(sim)[[currentModule(sim)]][[".studyAreaName"]])
  }

  studyArea <- sf::st_as_sf(sim$studyArea)
  studyArea_biomassParam <- sf::st_as_sf(sim$studyArea_biomassParam)

  ## this is necessary if studyArea and studyArea_biomassParam are multipolygon objects
  if (nrow(studyArea) > 1) {
    stop("please provide a study area that is not a multipolygon",
         "which will incorrectly segment ecoregions. Try `terra::aggregate`")
  }

  if (length(st_within(studyArea, studyArea_biomassParam))[[1]] == 0) {
    stop("studyArea is not fully within studyArea_biomassParam.
         Please check the aligment, projection and shapes of these polygons")
  }
  rm(studyArea, studyArea_biomassParam)

  if (!suppliedElsewhere("rasterToMatch", sim)) {
    studyArea <- sim$studyArea
    if (!inherits(studyArea, "SpatVector")) {
      studyArea <- vect(studyArea)
    }
    if (terra::is.lonlat(studyArea)) {
      ## use SCANFI projection - LandR requires projected rasters for dispersal
      studyArea <- project(studyArea,
                           paste("+proj=lcc +lat_0=0 +lon_0=-95 +lat_1=49 +lat_2=77",
                                 "+x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs"))
    }
    sim$rasterToMatch <- rast(studyArea, res = c(rtm_res, rtm_res), vals = 1) |>
      mask(mask = studyArea)
  }

  if (!suppliedElsewhere("rasterToMatch_biomassParam", sim)) {
    if (!is.null(sim$rasterToMatchLarge)) {
      warning("please use rasterToMatch_biomassParam in place of rasterToMatchLarge")
      sim$rasterToMatch_biomassParam <- sim$rasterToMatchLarge
    } else {
      sim$rasterToMatch_biomassParam <- sim$rasterToMatch
    }
  }

  ## biomass map
  if (!suppliedElsewhere("rawBiomassMap", sim)) {
    stopifnot(
      "dataYear must be one of 2000, 2010, 2020" = P(sim)$dataYear %in% c(2000, 2010, 2020)
    )

    sim$rawBiomassMap <- prepRawBiomassMap(
      dataSource = P(sim)$dataSource,
      dataYear = P(sim)$dataYear,
      to =  sim$rasterToMatch_biomassParam,
      destinationPath = dPath,
      writeTo = paste0("_", P(sim)$.studyAreaName, "_", P(sim)$dataYear, "_", P(sim)$dataSource) |>
        .suffix("biomass.tif", suffix = _)
    )
  }

  ## Land cover raster ------------------------------------------------
  if (!suppliedElsewhere("rstLCC", sim)) {
    sim$rstLCC <- Cache(
      prepInputs_SCANFI_LCC_FAO,
      year = P(sim)$dataYear,
      maskTo = sim$studyArea_biomassParam,
      cropTo = sim$rasterToMatch_biomassParam,
      projectTo = sim$rasterToMatch_biomassParam,
      disturbedCode = 240,
      destinationPath = dPath,
      overwrite = TRUE,
      writeTo = .suffix("rstLCC.tif", paste0("_", P(sim)$.studyAreaName, "_", P(sim)$dataYear)),
      userTags = c("rstLCC", currentModule(sim), P(sim)$.studyAreaName, P(sim)$dataYear)
    )
  }

  ## Ecodistrict ------------------------------------------------
  if (!suppliedElsewhere("ecoregionLayer", sim)) {
    ## Ceres: makePixel table needs same no. pixels for this, RTM rawBiomassMap, LCC.. etc
    sim$ecoregionLayer <- Cache(
      prepInputs(targetFile = "ecodistricts.shp",
                 archive = asPath("ecodistrict_shp.zip"),
                 url = extractURL("ecoregionLayer", sim),
                 alsoExtract = "similar",
                 destinationPath = dPath,
                 writeTo = NULL,
                 to = sim$studyArea_biomassParam,
                 fun = getOption("reproducible.shapefileRead"),
                 overwrite = TRUE),
      .functionName = "prepInputs_forEcoregionLayer",
      userTags = c("prepInputsEcoDistrict_SA", currentModule(sim), cacheTags)
    )
  }

  if (P(sim)$overrideAgeInFires) {
    if (!suppliedElsewhere("firePerimeters", sim)) {
      sa <- if (is(sim$studyArea_biomassParam, "sf")) {
        aggregate(sim$studyArea_biomassParam, list(rep(1, nrow(sim$studyArea_biomassParam))),
                  FUN = function(x) x)
      } else {
        aggregate(sim$studyArea_biomassParam)
      }
      sim$firePerimeters <- Cache(
        prepInputsFireYear,
        destinationPath = dPath,
        studyArea = sa,
        rasterToMatch = sim$rasterToMatch_biomassParam,
        overwrite = TRUE,
        url = extractURL("firePerimeters"),
        fireField = "YEAR",
        omitArgs = "destinationPath",
        userTags = c(cacheTags, "firePerimeters")
      )
      whichFiresTooOld <- which(as.vector(sim$firePerimeters[]) < P(sim)$earliestFireYear)

      if (length(whichFiresTooOld)) {
        message("There were fires in the database older than ", P(sim)$earliestFireYear, ";",
                " The data from these will not be used")
        sim$firePerimeters[whichFiresTooOld] <- NA
      }
    }
  }

  ## Stand age map ------------------------------------------------
  if (!suppliedElsewhere("standAgeMap", sim)) {
    sa <- if (is(sim$studyArea_biomassParam, "sf")) {
      aggregate(
        sim$studyArea_biomassParam,
        list(rep(1, nrow(sim$studyArea_biomassParam))),
        FUN = function(x) x
      )
    } else {
      aggregate(sim$studyArea_biomassParam)
    }

    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      sim$standAgeMap <- Cache(
        LandR::prepInputsStandAgeMap,
        dataSource = P(sim)$dataSource,
        dataYear = P(sim)$dataYear,
        ageFun = getOption("reproducible.rasterRead", "terra::rast"), ## backwards compatible default
        destinationPath = dPath,
        rasterToMatch = sim$rasterToMatch_biomassParam,
        # writeTo = .suffix("standAgeMap.tif", paste0("_", P(sim)$.studyAreaName)),
        overwrite = TRUE,
        useCache = FALSE, ## TODO: temporary FALSE due to attributes being lost on retrieval
        firePerimeters = if (P(sim)$overrideAgeInFires) sim$firePerimeters else NULL,
        fireURL = if (P(sim)$overrideAgeInFires) extractURL("firePerimeters") else NULL,
        startTime = P(sim)$dataYear,
        userTags = c("prepInputsStandAge_rtm", currentModule(sim), cacheTags),
        omitArgs = c("destinationPath", "targetFile", "overwrite",
                     "alsoExtract", "userTags")
      )
    })
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
                          sim$sppColorVect, P(sim)$vegLeadingProportion, sim$studyArea_biomassParam)
  ## the following may, or may not change inputs
  sim$sppEquiv <- sppOuts$sppEquiv
  sim$sppNameVector <- sppOuts$sppNameVector
  P(sim, module = currentModule(sim))$sppEquivCol <- sppOuts$sppEquivCol
  sim$sppColorVect <- sppOuts$sppColorVect

  ## check again
  paramCheckOtherMods(sim, "sppEquivCol", ifSetButDifferent = "error")

  ## Species raster layers -------------------------------------------
  if (!suppliedElsewhere("speciesLayers", sim)) {
    httr::with_config(config = httr::config(ssl_verifypeer = P(sim)$.sslVerify), {
      sim$speciesLayers <- Cache(prepSpeciesLayers_SCANFI,
                                 destinationPath = dPath,
                                 outputPath = dPath,
                                 studyArea = sim$studyArea_biomassParam,
                                 studyAreaName = P(sim)$.studyAreaName,
                                 rasterToMatch = sim$rasterToMatch_biomassParam,
                                 sppEquiv = sim$sppEquiv,
                                 sppEquivCol = P(sim)$sppEquivCol,
                                 thresh = 10,
                                 year = P(sim)$dataYear,
                                 userTags = c(cacheTags, "speciesLayers"),
                                 omitArgs = c("userTags"))
    })

    ## make sure empty pixels inside study area have 0 cover, instead of NAs.
    ## this can happen when data has NAs instead of 0s and is not merged/overlayed (e.g. CASFRI)
    sim$speciesLayers <- NAcover2zero(sim$speciesLayers, sim$rasterToMatch_biomassParam)
  }

  ## 3. species maps
  if (!suppliedElsewhere("speciesTable", sim)) {
    sim$speciesTable <- getSpeciesTable(dPath = dPath, cacheTags = c(cacheTags, "speciesTable"))
  }

  if (!suppliedElsewhere("columnsForPixelGroups", sim)) {
    sim$columnsForPixelGroups <- LandR::columnsForPixelGroups()
  }

  return(invisible(sim))
}

#' Probe NTEMS NFI web page to find the final year available
#'
#' Starts searching
#' `paste0("https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_", lastYrOnNTEMS, ".zip")`
#' at current year (`Sys.Date()`), and subtract one year, try, subtract a year, try etc.
#'
#' @param timeout Numeric, in seconds, for how long to allow a download to happen
#'   before interrupting it and declaring, "that worked, use that year".
NTEMSfinalYearForLCC <- function(timeout = 5) {
  resp <- ""
  lastYrOnNTEMS <- as.integer(format(Sys.Date(), "%Y")) + 1

  while (!is(resp, "try-error")) {
    lastYrOnNTEMS <- lastYrOnNTEMS - 1
    url <- paste0(
      "https://opendata.nfis.org/downloads/forest_change/CA_forest_VLCE2_",
      lastYrOnNTEMS,
      ".zip"
    )
    req <- httr2::request(url) |> httr2::req_timeout(timeout)
    resp <- try(httr2::req_perform(req), silent = TRUE)
  }

  lastYrOnNTEMS
}
