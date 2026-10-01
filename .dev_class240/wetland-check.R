#!/usr/bin/env Rscript
## Does SCANFI land cover + CWIM3A, with the 80/81 rule, reproduce NTEMS' wetland classes?
##
##   derived = wetlandToLCC(SCANFI 2020 land cover, CWIM wet)
##     wet & SCANFI in {210, 220, 230, 240} -> 81;  wet & other non-water -> 80
##   reference = NTEMS VLCE2 2020
##
## Everything on the harness's 240 m grid for each window, same year (2020; CWIM3A imagery is
## 2016-2020). Also reports the raw wet/not-wet agreement, and how much agreement a one-cell
## shift of the CWIM layer would buy -- a registration check.
##
## Usage: Rscript wetland-check.R <landrPath> [area ...]

suppressPackageStartupMessages({ library(terra); library(data.table) })
.args <- commandArgs(trailingOnly = TRUE)
landrPath <- .args[[1]]
areas <- if (length(.args) > 1) .args[-1] else c("big_s2020", "small_s2020")
ROOT <- Sys.getenv("CLASS240_ROOT", "/mnt/fast/class240-regression")
suppressPackageStartupMessages(pkgload::load_all(landrPath, quiet = TRUE, export_all = FALSE))
terraOptions(memmax = 4, todisk = TRUE, progress = 0)

scanfiFile <- file.path(ROOT, "inputs", "SCANFI_att_nfiLandcover_CanadaLCCclassCodes_2020_v2_20260119.tif")
vlceFile <- file.path(ROOT, "inputs", "CA_forest_VLCE2_2020.tif")
stopifnot(file.exists(scanfiFile), file.exists(vlceFile))

onGrid <- function(file, rtm, method) {
  r <- rast(file)
  win <- ext(project(as.polygons(ext(rtm), crs = crs(rtm)), crs(r)))
  win <- extend(win, 4 * max(res(r)))
  x <- crop(r, win)
  levels(x) <- NULL
  mask(project(x, rtm, method = method), rtm)
}
cls3 <- function(v) fifelse(v %in% 81, "81 treed wetland", fifelse(v %in% 80, "80 wetland", "other"))

res <- list()
for (area in areas) {
  rtm <- rast(file.path(ROOT, "fixtures", area, "rasterToMatch.tif"))
  scanfi <- onGrid(scanfiFile, rtm, "near")
  vlce   <- onGrid(vlceFile, rtm, "mode")
  wet    <- prepInputs_CWIM(rtm)
  derived <- wetlandToLCC(scanfi, wet)

  dt <- data.table(s = values(scanfi, mat = FALSE), v = values(vlce, mat = FALSE),
                   w = values(wet, mat = FALSE), d = values(derived, mat = FALSE))
  dt <- dt[!is.na(v) & !is.na(d)]
  n <- nrow(dt)

  vWet <- dt$v %in% c(80, 81); cWet <- dt$w %in% 1
  both <- sum(vWet & cWet)
  shiftAgree <- vapply(list(c(0, 0), c(1, 0), c(-1, 0), c(0, 1), c(0, -1)), function(sh) {
    ws <- terra::shift(wet, dx = sh[1] * res(wet)[1], dy = sh[2] * res(wet)[2])
    ws <- resample(ws, rtm, method = "near")
    wv <- values(ws, mat = FALSE)
    ok <- !is.na(values(vlce, mat = FALSE)) & !is.na(wv)
    mean((values(vlce, mat = FALSE)[ok] %in% c(80, 81)) == (wv[ok] %in% 1))
  }, numeric(1))

  conf <- dt[, .N, by = .(derived = cls3(d), NTEMS = cls3(v))]
  confW <- dcast(conf, derived ~ NTEMS, value.var = "N", fill = 0)

  treedAgree <- dt[v %in% 81, mean(d %in% 81)]
  treedPrec  <- dt[d %in% 81, mean(v %in% 81)]
  wetAgree   <- dt[v %in% 80, mean(d %in% 80)]
  wetPrec    <- dt[d %in% 80, mean(v %in% 80)]

  res[[area]] <- list(
    summary = data.table(
      area, cells = n,
      NTEMS_wet_pct = round(100 * mean(vWet), 1), CWIM_wet_pct = round(100 * mean(cWet), 1),
      wet_overall_agreement_pct = round(100 * mean(vWet == cWet), 1),
      NTEMS_wet_found_by_CWIM_pct = round(100 * both / max(sum(vWet), 1), 1),
      CWIM_wet_that_NTEMS_calls_wet_pct = round(100 * both / max(sum(cWet), 1), 1),
      NTEMS81_recalled_pct = round(100 * treedAgree, 1), derived81_precision_pct = round(100 * treedPrec, 1),
      NTEMS80_recalled_pct = round(100 * wetAgree, 1), derived80_precision_pct = round(100 * wetPrec, 1)
    ),
    confusion = confW,
    shift = data.table(area, shift = c("none", "+x", "-x", "+y", "-y"),
                       wet_agreement_pct = round(100 * shiftAgree, 2)),
    ntems81_underSCANFI = dt[v %in% 81, .N, by = .(SCANFI = s)][order(-N)],
    ntems81_notCWIMwet_underSCANFI = dt[v %in% 81 & !w %in% 1, .N, by = .(SCANFI = s)][order(-N)]
  )
}

f <- file.path(ROOT, "runs", "wetland-check.md")
con <- file(f, "w")
w <- function(...) { cat(..., "\n", sep = ""); cat(..., "\n", sep = "", file = con) }
w("# SCANFI + CWIM3A vs NTEMS 80/81 (2020, 240 m harness grids)\n")
w("```"); print(rbindlist(lapply(res, `[[`, "summary"))); capture.output(print(rbindlist(lapply(res, `[[`, "summary"))), file = con, append = TRUE); w("```\n")
for (a in names(res)) {
  w("## ", a, "\n")
  w("Confusion (rows: derived from SCANFI+CWIM; columns: NTEMS VLCE2 2020)\n")
  w("```"); print(res[[a]]$confusion); capture.output(print(res[[a]]$confusion), file = con, append = TRUE); w("```")
  w("Registration: wet/not-wet agreement with the CWIM layer shifted one cell\n")
  w("```"); print(res[[a]]$shift); capture.output(print(res[[a]]$shift), file = con, append = TRUE); w("```")
  w("What SCANFI calls the pixels NTEMS calls 81\n")
  w("```"); print(res[[a]]$ntems81_underSCANFI); capture.output(print(res[[a]]$ntems81_underSCANFI), file = con, append = TRUE); w("```")
  w("...of those, the ones CWIM does NOT call wet\n")
  w("```"); print(res[[a]]$ntems81_notCWIMwet_underSCANFI); capture.output(print(res[[a]]$ntems81_notCWIMwet_underSCANFI), file = con, append = TRUE); w("```\n")
}
close(con)
message("written: ", f)
