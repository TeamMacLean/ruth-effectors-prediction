data <- data.frame(
  cnn_lstm = rbinom(20, 1, 0.3),
  cnn_gru  = rbinom(20, 1, 0.25),
  gru_emb  = rbinom(20, 1, 0.1),
  lstm_emb = rbinom(20, 1, 0.4),
  ensemble = rbinom(20, 1, 0.4)
)

plot_venn_diagram_limma <- function(data) {
  # Store original column names
  model_names <- colnames(data)

  venn_counts <- data %>%
    # Calculate counts
    dplyr::mutate_all(., as.logical) %>%
    limma::vennCounts() %>%
    `class<-`("matrix") %>%
    as.data.frame() %>%
    # Rename columns
    `colnames<-`(c("A", "B", "C", "D", "E", "Counts")) %>%
    # Group 0000 and 1111 cases
    mutate(Counts = ifelse(A == B & A == C & A == D & A == E,
                           sum(Counts[A == B & A == C & A == D & A == E]),
                           Counts)) %>%
    slice(-1) %>%
    # Recover original names
    `colnames<-`(c(model_names, "Counts")) %>%
    # Transform back to VennCounts class
    as.matrix() %>%
    `class<-`("VennCounts")

  venn_counts %>%
    limma::vennDiagram(
      circle.col = c("#f35e5a", "#929705", "#18b56a", "#149ffe", "#de4af0")
    )
}

plot_venn_diagram_limma(data)
