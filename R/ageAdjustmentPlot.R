ageAdjustmentPlot <- function(df, longevity) {
  longevity$lty <- factor(1, levels = 1, labels = "Longevity")
  ggplot(df, aes(x = age, fill = processed)) +
    geom_histogram(alpha = 0.5, position = position_identity()) +
    geom_vline(data = longevity, aes(xintercept = longevity_orig, linetype = lty), linewidth = 1, col = "grey50") +
    geom_vline(data = longevity, aes(xintercept = longevity, linetype = lty), linewidth = 1, col = "#d95f02") +
    scale_fill_manual(values = c("#d95f02", "grey50")) +
    scale_linetype_manual(NULL, values = "dashed") +
    facet_wrap(~ speciesCode, ncol = 4, scales = "free") +
    theme_bw() +
    labs(y = "Number of cohorts", x = "Age", fill = NULL)
}
