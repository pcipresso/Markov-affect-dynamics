library(dplyr)

# ============================================================
# Optional preprocessing example:
# from item-level PANAS data to state sequence, C, and P
# by Pietro Cipresso & Francesca Borghesi
# ============================================================

# ------------------------------------------------------------
# OPTION 1 (default): Simulated PANAS data
# ------------------------------------------------------------

# Comment this block if you want to use your own CSV data

set.seed(123)

n_time <- 180

panas_raw <- data.frame(
  subject = rep("S1", n_time),
  time = 1:n_time,
  
  interested   = sample(1:5, n_time, replace = TRUE),
  excited      = sample(1:5, n_time, replace = TRUE),
  strong       = sample(1:5, n_time, replace = TRUE),
  enthusiastic = sample(1:5, n_time, replace = TRUE),
  alert        = sample(1:5, n_time, replace = TRUE),
  
  distressed = sample(1:5, n_time, replace = TRUE),
  upset      = sample(1:5, n_time, replace = TRUE),
  guilty     = sample(1:5, n_time, replace = TRUE),
  scared     = sample(1:5, n_time, replace = TRUE),
  hostile    = sample(1:5, n_time, replace = TRUE)
)

# ------------------------------------------------------------
# OPTION 2 (alternative): Load your own CSV file
# ------------------------------------------------------------

# Uncomment and adapt this line to use real data.
# The CSV must contain one row per ESM prompt.
# Required columns: subject, time, and the PANAS item columns specified below.

# panas_raw <- read.csv("your_file.csv", stringsAsFactors = FALSE)

# ------------------------------------------------------------
# Define PANAS item names
# ------------------------------------------------------------

# Adapt these names to match the column names in your CSV file.

pa_items <- c("interested", "excited", "strong", "enthusiastic", "alert")
na_items <- c("distressed", "upset", "guilty", "scared", "hostile")

# ------------------------------------------------------------
# Input checks
# ------------------------------------------------------------

required_cols <- c("subject", "time", pa_items, na_items)
missing_cols <- setdiff(required_cols, names(panas_raw))

if (length(missing_cols) > 0) {
  stop(
    "The following required columns are missing from panas_raw: ",
    paste(missing_cols, collapse = ", ")
  )
}

if (any(is.na(panas_raw$subject))) {
  stop("The subject column contains missing values.")
}

if (any(is.na(panas_raw$time))) {
  stop("The time column contains missing values.")
}

# Convert PANAS item columns to numeric if needed.
# This is useful when item responses are imported from CSV as character strings.

panas_raw <- panas_raw %>%
  mutate(
    across(
      all_of(c(pa_items, na_items)),
      ~ as.numeric(.x)
    )
  )

# Check whether item conversion produced missing values.
# This may indicate non-numeric entries in the CSV file.

if (any(is.na(panas_raw[, c(pa_items, na_items)]))) {
  warning(
    "Some PANAS item values are missing or non-numeric. ",
    "PA and NA scores will be computed using available items only."
  )
}

# ------------------------------------------------------------
# Compute PA and NA scores
# ------------------------------------------------------------

panas_scored <- panas_raw %>%
  mutate(
    PA_score = rowMeans(select(., all_of(pa_items)), na.rm = TRUE),
    NA_score = rowMeans(select(., all_of(na_items)), na.rm = TRUE)
  )

# Remove rows where PA or NA cannot be computed

panas_scored <- panas_scored %>%
  filter(!is.na(PA_score), !is.na(NA_score))

if (nrow(panas_scored) == 0) {
  stop("No valid rows remain after computing PA_score and NA_score.")
}

# ------------------------------------------------------------
# Discretization into affective states (within-person)
# ------------------------------------------------------------

panas_states <- panas_scored %>%
  group_by(subject) %>%
  mutate(
    PA_cut = median(PA_score, na.rm = TRUE),
    NA_cut = median(NA_score, na.rm = TRUE),
    PA_level = ifelse(PA_score >= PA_cut, "high_PA", "low_PA"),
    NA_level = ifelse(NA_score >= NA_cut, "high_NA", "low_NA"),
    state = case_when(
      PA_level == "high_PA" & NA_level == "low_NA"  ~ "Positive-dominant",
      PA_level == "high_PA" & NA_level == "high_NA" ~ "Mixed-high",
      PA_level == "low_PA"  & NA_level == "high_NA" ~ "Negative-dominant",
      PA_level == "low_PA"  & NA_level == "low_NA"  ~ "Low-affect",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# ------------------------------------------------------------
# Build state sequence
# ------------------------------------------------------------

state_sequence <- panas_states %>%
  select(subject, time, PA_score, NA_score, state) %>%
  arrange(subject, time) %>%
  filter(!is.na(state))

# Check whether each subject has at least two valid observations

valid_counts <- state_sequence %>%
  count(subject, name = "n_valid_timepoints")

if (any(valid_counts$n_valid_timepoints < 2)) {
  warning(
    "At least one subject has fewer than two valid time points and cannot contribute transitions."
  )
}

# ------------------------------------------------------------
# Build transitions
# ------------------------------------------------------------

transitions <- state_sequence %>%
  group_by(subject) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    next_state = lead(state),
    next_time = lead(time)
  ) %>%
  ungroup() %>%
  filter(!is.na(next_state))

if (nrow(transitions) == 0) {
  stop("No valid transitions could be constructed.")
}

# ------------------------------------------------------------
# Define state order
# ------------------------------------------------------------

state_levels <- c(
  "Positive-dominant",
  "Mixed-high",
  "Negative-dominant",
  "Low-affect"
)

# ------------------------------------------------------------
# Transition count matrix C
# ------------------------------------------------------------

C <- table(
  factor(transitions$state, levels = state_levels),
  factor(transitions$next_state, levels = state_levels)
)

C <- as.matrix(C)

# ------------------------------------------------------------
# Transition probability matrix P
# ------------------------------------------------------------

row_totals <- rowSums(C)

if (any(row_totals == 0)) {
  warning(
    "At least one state has no outgoing transitions. ",
    "Rows with zero outgoing transitions will be set to NA in P."
  )
}

P <- matrix(
  NA_real_,
  nrow = nrow(C),
  ncol = ncol(C),
  dimnames = dimnames(C)
)

for (i in seq_len(nrow(C))) {
  if (row_totals[i] > 0) {
    P[i, ] <- C[i, ] / row_totals[i]
  }
}

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

cat("\n--- First rows of state sequence ---\n")
print(head(state_sequence, 12))

cat("\n--- Transition count matrix C ---\n")
print(C)

cat("\n--- Transition probability matrix P ---\n")
print(round(P, 3))