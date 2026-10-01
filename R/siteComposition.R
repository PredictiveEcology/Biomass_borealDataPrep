## Stratum codes for grouping pixels when parameters are estimated.
##
## `ecoregionGroup` is "<ecoregion>_<code>". Two stratifications are offered:
##
## - "landcover": the code is the land-cover class, with the NTEMS wetland classes added from
##   the site layer -- one axis, as NTEMS itself has it. Treed wetland is 81 whatever it is
##   made of.
## - "siteComposition": site and composition are crossed, so black spruce on upland and black
##   spruce on peatland are different strata with their own maxB, maxANPP and establishment
##   probability. Neither axis is read off the other.
##
## Every code is three digits. LandR's convertUnwantedLCC() and the module's
## paddedFloatToChar() calls size the code field from the widest code present, so a four-digit
## code would re-pad every existing one, and the `_240` assertions would stop matching.
##
##   upland   210 conifer   220 broadleaf   230 mixed
##   wet      810 conifer   820 broadleaf   830 mixed      (upland + 600)
##   pooled   290 upland, any composition   890 wet, any composition   990 any site
##   no species cover to decide composition: 240 upland, 840 wet
.treedComposition <- c(210L, 220L, 230L)
.wetOffset <- 600L
.pooledUpland <- 290L
.pooledWet <- 890L
.pooledAny <- 990L

## Codes for each pixel from its composition and its site.
##
## `comp` is the land-cover class with class-240 pixels already resolved to 210/220/230 where
## species cover allowed (unresolved ones are still `unresolved`). `wet` is the site layer:
## non-zero is wet; 0 or NA is not.
siteCompositionCodes <- function(comp, wet, stratumType = c("landcover", "siteComposition"),
                                 unresolved = 240L) {
  stratumType <- match.arg(stratumType)
  if (length(comp) != length(wet)) {
    stop("`comp` and `wet` must have the same length (", length(comp), " vs ", length(wet), ")")
  }
  isWet <- !is.na(wet) & wet != 0

  if (identical(stratumType, "landcover")) {
    ## wet and treed -> 81 (240 included: forest land not currently stocked); wet otherwise -> 80
    return(LandR::wetlandToLCC(comp, as.integer(isWet),
                               treedClasses = c(.treedComposition, unresolved)))
  }

  out <- comp
  wetTreed <- isWet & comp %in% c(.treedComposition, unresolved)
  out[wetTreed] <- comp[wetTreed] + .wetOffset
  wetOther <- isWet & !wetTreed & !is.na(comp) & !comp %in% c(20, 80, 81)
  out[wetOther] <- 80
  out
}

## The classes that still need a neighbour class after composition has been resolved.
unresolvedClasses <- function(classesToReplace, stratumType = c("landcover", "siteComposition")) {
  stratumType <- match.arg(stratumType)
  if (identical(stratumType, "siteComposition")) {
    unique(c(classesToReplace, classesToReplace + .wetOffset))
  } else {
    classesToReplace
  }
}

## Pool strata that have too few pixels to estimate from.
##
## `counts` has one row per ecoregion x stratum with `N`, the number of *estimation* pixels
## (inferred pixels excluded). A thin stratum loses its composition first -- it joins the
## site's pooled class (290 or 890) -- and, if that pooled class is still thin, its site too
## (990). Strata with enough pixels keep their own code, so a pooled class holds only the
## thin ones: it is a residual, not an average over the site.
##
## Returns `counts` with `newStratum`.
collapseThinStrata <- function(counts, minN) {
  m <- data.table::copy(data.table::as.data.table(counts))
  stopifnot(all(c("eco", "stratum", "N") %in% names(m)))
  m[, newStratum := as.integer(stratum)]
  if (!nrow(m) || !is.finite(minN) || minN <= 0) {
    return(m[])
  }
  m[N < minN & stratum %in% .treedComposition, newStratum := .pooledUpland]
  m[N < minN & stratum %in% (.treedComposition + .wetOffset), newStratum := .pooledWet]
  lvl2 <- m[, list(N2 = sum(N)), by = c("eco", "newStratum")]
  m[lvl2, N2 := i.N2, on = c("eco", "newStratum")]
  m[newStratum %in% c(.pooledUpland, .pooledWet) & N2 < minN, newStratum := .pooledAny]
  m[, N2 := NULL]
  m[]
}

## Final stratum for pixels whose class was inferred.
##
## Inferred pixels are not used for estimation, so each needs a stratum that estimation
## pixels populate. Its own stratum is used if present (as pooled by collapseThinStrata());
## otherwise, when `allowPooling`, its site's pooled class, then 990; otherwise it goes back to
## the unresolved code for its site, for convertUnwantedLCC() to settle from its neighbours.
stratumForInferred <- function(eco, stratum, mapping, allowPooling = TRUE, unresolved = 240L) {
  if (length(eco) != length(stratum)) stop("`eco` and `stratum` must have the same length")
  out <- rep(NA_integer_, length(stratum))
  if (!length(stratum)) return(out)

  hit <- match(paste(eco, stratum), paste(mapping$eco, mapping$stratum))
  out[!is.na(hit)] <- mapping$newStratum[hit[!is.na(hit)]]

  isWet <- stratum >= 800 | stratum == 81
  if (allowPooling) {
    have <- unique(paste(mapping$eco, mapping$newStratum))
    pooled <- ifelse(isWet, .pooledWet, .pooledUpland)
    use2 <- is.na(out) & paste(eco, pooled) %in% have
    out[use2] <- pooled[use2]
    use1 <- is.na(out) & paste(eco, .pooledAny) %in% have
    out[use1] <- .pooledAny
  }
  back <- is.na(out)
  out[back] <- if (allowPooling) {
    ifelse(isWet[back], unresolved + .wetOffset, unresolved)
  } else {
    unresolved
  }
  out
}
