# Packages needed:
# ----------------
# ggplot2
# dplyr
# tidyr
# pheatmap
# stats
# viridis


# Function to plot correlation matrices

plot_cormat <- function(rfm_df, cor_trans = NULL, variable = NULL) {
  rfm_df <- rfm_df %>%
    dplyr::select_if(is.numeric) %>%
    stats::cor(use = "pairwise.complete.obs") %>%
    # Transform to data frame
    tibble::as_tibble(rownames = "var_x") %>%
    # Prepare for plotting
    tidyr::pivot_longer(
      cols = -var_x,
      names_to = "var_y",
      values_to = "cor"
    ) %>%
    dplyr::mutate(
      var_x = factor(var_x, levels = unique(.[["var_x"]])),
      var_y = factor(var_y, levels = unique(.[["var_x"]]))
    )

  if (!is.null(variable)) {
    rfm_df <- rfm_df %>%
      dplyr::filter(var_x == variable)
  }

  # Transform correlation
  if (is.null(cor_trans)) {
    fill_limits <- c(-1, 1)
  } else {
    if (cor_trans == "abs") {
      rfm_df <- rfm_df %>%
        dplyr::mutate(cor = abs(cor))
      fill_limits <- c(0, 1)
    } else if (cor_trans == "squared") {
      rfm_df <- rfm_df %>%
        dplyr::mutate(cor = cor^2)
      fill_limits <- c(0, 1)
    }
  }

  # Plot
  gg <- rfm_df %>%
    ggplot() +
    aes(x = var_x, y = var_y, fill = cor) +
    geom_tile() +
    viridis::scale_fill_viridis(limits = fill_limits) +
    coord_fixed() +
    labs(
      x = NULL,
      y = NULL
    ) +
    geom_text(aes(label = round(cor, 2)), vjust = 0.5, size=3) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_blank())

  return(gg)
}

# Function to plot clustered corrlation matrix

plot_clustered_correlation <- function(data, cellsize, ncolor, ...) {
  data %>%
    stats::cor(method = "pearson") %>%
    pheatmap::pheatmap(
      mat = .,
      color = viridisLite::viridis(ncolor, begin = 0, end = 1),
      cellheight = cellsize,
      cellwidth = cellsize,
      border_color = NA,
      display_numbers = TRUE,
      breaks = seq(-1, 1, 2 / ncolor),
      ...
    )
}
