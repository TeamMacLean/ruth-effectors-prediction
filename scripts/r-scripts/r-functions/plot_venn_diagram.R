# Packages needed:
# ----------------
# dplyr
# limma
# tidyr


plot_venn_diagram_limma <- function(data) {
  # Store original column names
  model_names <- colnames(data)

  color_list <- c("#f35e5a", "#929705", "#18b56a", "#149ffe", "#de4af0", "#de4af0")

  venn_counts <- data %>%
    # Calculate counts
    dplyr::mutate_all(., as.logical) %>%
    limma::vennCounts() %>%
    `class<-`("matrix") %>%
    as.data.frame() %>%
    # Rename columns
    tidyr::unite(
      col = lookup_string,
      model_names,
      sep = "_",
      remove = FALSE
    ) %>%
    # Group 0000 and 1111 cases
    mutate(
      lookup_string = ifelse(
        lookup_string == paste0(rep("0", length(model_names)), collapse = "_"),
        paste0(rep("1", length(model_names)), collapse = "_"),
        lookup_string
      )
    ) %>%
    select(-model_names) %>%
    tidyr::separate(
      col = lookup_string,
      into = model_names,
      sep = "_",
      remove = FALSE
    ) %>%
    group_by_at(vars(model_names)) %>%
    summarise(
      Counts = sum(Counts)
    ) %>%
    # Transform back to VennCounts class
    ungroup() %>%
    rbind(
      data.frame(
        t(c(rep(0, length(model_names)),NA))
      ) %>%
        `colnames<-`(c(model_names, "Counts"))
    ) %>%
    tidyr::unite(
      col = row_num,
      model_names,
      sep = "",
      remove = FALSE
    ) %>%
    mutate(
      row_num = strtoi(row_num, base = 2)
    ) %>%
    arrange(row_num) %>%
    select(-row_num) %>%
    as.matrix() %>%
    `class<-`("VennCounts")

  venn_counts %>%
    limma::vennDiagram(
      circle.col = color_list[1:length(model_names)],
      include=c("up","down")
    )
}
