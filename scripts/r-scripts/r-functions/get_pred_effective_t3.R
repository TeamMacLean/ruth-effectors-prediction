# Function
get_pred_effective_t3 <- function(fasta_path, module = "STD") {
  java_path <- "/Users/kristian/Documents/Workspace/Software/TTSS/"

  module_file <- Sys.glob(paste0(java_path, "module/*", module, "*.jar")) %>%
    stringr::str_split("/") %>%
    unlist() %>%
    .[[length(.)]]

  command <- glue::glue(
    "cd ", java_path, ";",
    "java -jar TTSS_GUI-1.0.1.jar -f ", fasta_path, " -m ", module_file, " -t selective -q"
    # "java -jar TTSS_GUI-1.0.1.jar -f ", fasta_path, " -m ", module_file, " -t cutoff=0.5 -q"
  )

  raw_data <- system(command, intern = TRUE, ignore.stderr = TRUE)

  col_names <- raw_data %>%
    # Get rid of header
    .[2] %>%
    stringr::str_split(";") %>%
    unlist() %>%
    stringr::str_trim() %>%
    stringr::str_to_lower() %>%
    stringr::str_replace_all(" ", "_")

  data <- raw_data %>%
    # Get rid of header
    .[-c(1, 2)] %>%
    stringr::str_split(";") %>%
    purrr::reduce(rbind) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    `colnames<-`(col_names) %>%
    mutate(
      score = as.numeric(score),
      is_effective = (is_effective == "true")
    ) %>%
    dplyr::as_tibble()

  return(data)
}


# This how this can be used
# bacteria_preds <- get_pred_effective_t3("/Users/kristian/Documents/Workspace/ruth-effectors-prediction/data/secreted_data/ready_to_process/fasta_files/bacteria_testing.fasta")
