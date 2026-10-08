## Estimate how much biomass a unit of deciduous cover carries, relative to a unit of conifer
## cover -- the weight `LandR::partitionBiomass()` applies when it splits a pixel's total biomass
## among its cohorts:
##
##     cover2 = cover * (w if deciduous, 1 otherwise);  B = totalBiomass * cover2 / sum(cover2)
##
## so w is (biomass per unit cover, deciduous) / (biomass per unit cover, conifer). Broadleaf
## crowns are wider per unit of wood, which is why it was ever expected to be below 1.
##
## `partitionBiomass()` on its own does not identify w: it is a partition, and a partition says
## nothing about the total it divides. The identifying statement is about the total,
##
##     totalBiomass_p = k_p * (C_p + w * D_p)
##
## with C and D the pixel's conifer and deciduous cover. What k is decides what is estimated:
##
##  * k constant within an ecoregion -- the obvious choice -- makes w absorb site quality. In the
##    boreal, broadleaf sits on richer, better-drained ground than black spruce, so the estimate
##    comes back ABOVE 1 (1.7 on 100 km of Alberta mixedwood): it has credited the site to the
##    species. That is not what this parameter means.
##  * k a function of the pixel's own structure -- SCANFI's canopy height and closure -- compares
##    composition between stands that carry the SAME amount of structure. Site quality shows up as
##    height and closure, so conditioning on them removes it, and what is left of the composition
##    signal is the crown-architecture difference the parameter is for. The same Alberta window
##    then gives 0.90.
##
## Cost: f and the stratum effects do not depend on w, so the design matrix is built once and each
## candidate w costs one `qr.resid()`. ~5 s on 158,000 pixels, against ~200 s for the AIC search
## this replaced -- cheap enough to run every time rather than carry a hardcoded number.
##
## CAVEAT on which structure layers are used. SCANFI's biomass and its height/closure come out of
## the same imputation, and it shows: on 100 km of Alberta mixedwood, SCANFI biomass regressed on
## SCANFI height, closure, age and ecoregion has R2 = 0.957 and a residual SD of 0.09 in log space.
## Conditioning on structure therefore leaves very little variance for composition to explain, and
## what is left is partly the internal structure of that one product. Feeding the same fit an
## independent biomass layer -- NTEMS 2015, same window, same SCANFI structure controls, R2 = 0.746
## and residual SD 0.38 -- gives 0.65 rather than 0.97. Both are below 1, so the sign is robust, but
## the size is not settled, and the SCANFI-on-SCANFI number is the more circular of the two.
##
## IDENTIFIABILITY. The model contributes `log1p((w - 1) * d)` per pixel, `d` the deciduous share
## of cover. If `d` is nearly the same in every pixel that term is nearly constant, the intercept
## absorbs it, and `w` is unidentified however many pixels there are -- a uniformly deciduous
## landscape is as uninformative as one with no deciduous at all. What identifies `w` is how much
## `d` varies between pixels, `sd(d)`, not how much deciduous there is. Neither the number of
## pixels nor a bootstrap catches this: across seven ELFs and windows (2026-09-18) the landscapes
## with sd(d) of 0.033-0.036 returned 0.63, 0.85 and 1.12 with bootstrap SDs of 0.001-0.012,
## i.e. precisely wrong, while those with sd(d) of 0.30-0.38 all returned 0.90-1.00.
## `minDeciduousShareSD = 0.1` sits in that gap (on a log scale, midway between 0.036 and 0.30).
##
## Returns `NA_real_` when the landscape cannot identify it (too few deciduous pixels, too little
## between-pixel variation in deciduous share, or an estimate on the search bound), and the caller
## keeps `P(sim)$deciduousCoverWeight`.

#' @importFrom data.table as.data.table
#' @importFrom stats as.formula model.matrix optimize sd
#' @importFrom utils data
deciduousCoverWeightFn <- function(pixelCohortData, canopyHeight, canopyClosure,
                                   interval = c(0.05, 5),
                                   minDeciduousPixels = 500L,
                                   minDeciduousCover = 0.02,
                                   minDeciduousShareSD = 0.1,
                                   minHeight = 2, minClosure = 5,
                                   dfStruct = 5L) {
  pcd <- as.data.table(pixelCohortData)
  ## `decid` exactly as LandR::partitionBiomass() decides it -- the packaged table, not
  ## `sim$sppEquiv` -- because that is the function that consumes the answer.
  sppEq <- get(data("sppEquivalencies_CA", package = "LandR", envir = environment()), inherits = FALSE)
  colName <- LandR::equivalentNameColumn(as.character(unique(pcd$speciesCode)), sppEq)
  decidSp <- LandR::equivalentName(sppEq[sppEq$Broadleaf == TRUE, ][["LandR"]], sppEq, colName)
  decidSp <- decidSp[nzchar(decidSp)]

  pix <- pcd[, list(totalBiomass = totalBiomass[1L], age = age[1L],
                    eco = as.character(initialEcoregionCode)[1L],
                    D = sum(cover[as.character(speciesCode) %in% decidSp]),
                    C = sum(cover[!as.character(speciesCode) %in% decidSp])),
             by = "pixelIndex"]
  pix[, eco := sub("_.*", "", eco)]   ## ecoregion, not ecoregion x land cover: land cover is itself
                                      ## a statement about composition, so conditioning on it would
                                      ## block part of the signal being estimated.
  pix[, height := as.vector(canopyHeight[])[pixelIndex]]
  pix[, closure := as.vector(canopyClosure[])[pixelIndex]]

  pix <- pix[totalBiomass > 0 & (C + D) > 0 & !is.na(height) & !is.na(closure) &
               height >= minHeight & closure >= minClosure]
  pix[, d := D / (C + D)]

  nDecid <- sum(pix$D > 0)
  if (NROW(pix) < minDeciduousPixels || nDecid < minDeciduousPixels ||
      mean(pix$d) < minDeciduousCover) {
    message(cli::col_yellow(
      "  not enough deciduous cover to estimate the deciduous cover weight (",
      nDecid, " pixels with any, mean cover share ", round(mean(pix$d) * 100, 2), "%)"))
    return(NA_real_)
  }
  sdD <- stats::sd(pix$d)
  if (sdD < minDeciduousShareSD) {
    message(cli::col_yellow(
      "  deciduous share of cover barely varies between pixels (sd ", round(sdD, 3), " < ",
      minDeciduousShareSD, ", mean ", round(mean(pix$d) * 100, 1), "%), so the deciduous cover ",
      "weight cannot be identified here"))
    return(NA_real_)
  }

  ## `splines::ns`, not `ns`: a SpaDES module's functions are sourced into an environment where
  ## `splines` is not attached, so a bare `ns` is only found once the module has been rendered as
  ## a package (which is how the tests run it).
  rhs <- paste("factor(eco) + splines::ns(height, dfStruct) * splines::ns(closure, dfStruct) +",
               "splines::ns(log1p(age), dfStruct)")
  qrX <- qr(model.matrix(as.formula(paste("~", rhs)), data = pix))
  y <- log(pix$totalBiomass)
  dd <- pix$d
  rss <- function(w) sum(qr.resid(qrX, y - log1p((w - 1) * dd))^2)
  out <- optimize(rss, interval = interval)

  atBound <- isTRUE(all.equal(out$minimum, interval[1L], tolerance = 1e-3)) ||
    isTRUE(all.equal(out$minimum, interval[2L], tolerance = 1e-3))
  if (atBound) {
    message(cli::col_yellow("  deciduous cover weight hit the search bound (",
                            round(out$minimum, 3), "); treat it as unidentified"))
    return(NA_real_)
  }
  message(cli::col_blue(
    "  deciduous cover weight estimated at ", round(out$minimum, 4),
    " from ", format(NROW(pix), big.mark = ","), " pixels (", nDecid, " with deciduous), ",
    "relative fit at 1.0 = ", round(rss(1) / out$objective, 4)))
  out$minimum
}
