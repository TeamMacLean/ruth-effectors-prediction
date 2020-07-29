get_formatted_date <- function(curr_date, to_date = NULL, range = FALSE, comment = NULL) {
  if (is.null(to_date)) {
    date_text <- curr_date %>%
      lubridate::ymd() %>%
      base::format("%d %B %Y (%A)")

    date_text <- paste(
      date_text,
      ifelse(
        is.null(comment),
        "",
        paste("–", stringr::str_to_title(comment))
      )
    )
  } else {
    days_list <- seq(as.Date(curr_date), as.Date(to_date), 1)
    if (!range) {
      days_text <- days_list %>%
        lubridate::ymd() %>%
        base::format("%d") %>%
        stringr::str_c(collapse = ", ")

      weekdays_text <- days_list %>%
        lubridate::ymd() %>%
        base::format("%A") %>%
        stringr::str_c(collapse = ", ")

      weekdays_text <- paste0(" (", weekdays_text, ") ")
    } else {
      days_text <- days_list %>%
        lubridate::ymd() %>%
        base::format("%d")

      days_text <- days_text[c(1, length(days_text))] %>%
        stringr::str_c(collapse = "–")

      weekdays_text <- ""
    }

    monthyear_text <- days_list %>%
      lubridate::ymd() %>%
      base::format("%B %Y") %>%
      unique()

    date_text <- paste0(
      days_text, " ",
      monthyear_text,
      weekdays_text,
      ifelse(
        is.null(comment),
        "",
        paste("–", stringr::str_to_title(comment))
      )
    )
  }

  return(date_text)
}

# Get total count of commits
get_git_commit_count_from_path_list <- function(git_path_list, ...) {
  git_count <- git_path_list %>%
    purrr::map(
      .f = function(path) {
        gitlogr::get_git_commit_count(path = path, ...)
      }
    ) %>%
    purrr::reduce(sum)

  return(git_count)
}
