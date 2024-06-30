# Fix excel dates ----
convert_date_columns <- function(data, date_columns, origin = "1899-12-30") {
  suppressWarnings({
    for (column_name in date_columns) {
      data$date_correct <- dmy(data[[column_name]])
      inds <- is.na(data$date_correct)
      data$date_correct[inds] <- as.Date(as.numeric(data[[column_name]][inds]), origin = origin)
      data[[column_name]] <- data$date_correct
      data$date_correct <- NULL
    }
  })
  return(data)
}


# Add dates of outcomes assessment columns ----
add_assessment_dates <- function(data) {
  data %>%
    mutate(clinical_assessment_date = visit_date,
           endoscopic_assessment_date = endoscopy_date) %>%
    mutate(across(c(clinical_assessment_date, endoscopic_assessment_date), ~ if_else(visit_timepoint %in% c("week_24", "week_12-16"), ., NA)))
}

# Format p value as required by CGH ----
format_p_value <- function(p_value) {
  if (p_value < 0.001) {
    return(bquote(italic(P) < ".001"))
  } else if (p_value < 0.1) {
    formatted_p_value <- formatC(p_value, format = "f", digits = 3)
  } else {
    formatted_p_value <- formatC(p_value, format = "f", digits = 2)
  }
    formatted_p_value <- sub("^0+", "", formatted_p_value)
    return(bquote(italic(P) == .(formatted_p_value)))
}

# Convert logits to probabilities ----
logit2prob <- function(logit){
  odds <- exp(logit)
  prob <- odds / (1 + odds)
  return(prob)
}

# Convert numeric age to Montreal categories ----
get_age_montreal <- function(data) {
  data %<>%
    mutate(age_montreal = case_when(
      age_at_diagnosis_yrs < 17 ~ "A1",
      age_at_diagnosis_yrs > 40 ~ "A3",
      .default = "A2"
    ), .after = age_at_diagnosis_yrs)
  return(data)
}

# Recode previous_advanced therapies ----
recode_previous_advancedrx <- function(data){
  data %<>% 
    mutate(previous_advancedrx_amount = as.character(previous_advancedrx_amount),
           previous_advancedrx_amount = if_else(previous_advancedrx_amount %in% c(3, 4, 5), "≥ 3", previous_advancedrx_amount))
}

# Recode steroids ----
recode_steroids <- function(data) {
  data %<>% 
    mutate(concomitant_cs = case_when(
      concomitant_cs == "0" ~ "no",
      str_detect(concomitant_cs, regex("medrol", ignore_case = TRUE)) ~ "sys",
      .default = "top"
    ))
}

# Get percentages of outcomes ----
outcome_pct <- function(data, outcomes) {
  data %>%
    summarise(across(any_of({{outcomes}}), 
                     list(n = ~sum(.), 
                          total = ~n(), 
                          percentage = ~round(sum(.) / n() * 100, 1)))) %>%
    pivot_longer(everything(), 
                 names_to = c("outcome", ".value"), 
                 names_pattern = "(.*_*)_(.*)") %>%
    mutate(report = str_glue("{n}/{total} ({percentage}%)")) %>%
    select(outcome, report)
}
