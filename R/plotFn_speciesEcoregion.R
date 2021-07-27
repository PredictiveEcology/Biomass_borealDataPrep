plotFn_speciesEcoregion <- function(stk, SEtype) {
  ggSpeciesEcoregion <-
    gplot(stk, maxpixels = 2e6) +
    geom_tile(aes(fill = value)) +
    facet_wrap(~ variable) +
    scale_fill_distiller(palette = "YlGnBu", na.value = "white", direction = 1) +
    theme_bw() +
    coord_equal() +
    ggtitle(SEtype)
}
