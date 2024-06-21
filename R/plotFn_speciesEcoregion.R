plotFn_speciesEcoregion <- function(stk, SEtype) {
  ggplot() +
    tidyterra::geom_spatraster(data = stk) +
    facet_wrap(~ lyr) +
    scale_fill_distiller(palette = "YlGnBu", direction = 1, na.value = "transparent") +
    theme_bw() +
    ggtitle(SEtype)
}
