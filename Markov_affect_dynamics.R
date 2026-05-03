# ============================================================
# Discrete-time Markov chains for affect dynamics
# Worked example and reusable functions
# by Pietro Cipresso & Francesca Borghesi
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "igraph",
  "scales"
)

to_install <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)

library(dplyr)
library(tidyr)
library(ggplot2)
library(igraph)
library(scales)

# ============================================================
# Basic functions
# ============================================================

row_normalize <- function(count_matrix) {
  rs <- rowSums(count_matrix)
  P <- matrix(
    NA_real_,
    nrow = nrow(count_matrix),
    ncol = ncol(count_matrix),
    dimnames = dimnames(count_matrix)
  )
  
  for (i in seq_len(nrow(count_matrix))) {
    if (rs[i] > 0) {
      P[i, ] <- count_matrix[i, ] / rs[i]
    }
  }
  
  P
}

compute_stationary_distribution <- function(P, tol = 1e-10) {
  if (any(is.na(P))) return(rep(NA_real_, nrow(P)))
  
  ev <- eigen(t(P))
  idx <- which.min(abs(ev$values - 1))
  vec <- Re(ev$vectors[, idx])
  
  if (all(vec <= 0)) vec <- -vec
  vec[vec < 0] <- 0
  
  if (sum(vec) == 0) return(rep(NA_real_, nrow(P)))
  
  pi <- vec / sum(vec)
  names(pi) <- rownames(P)
  
  if (max(abs(as.numeric(pi %*% P - pi))) > sqrt(tol)) {
    warning("Stationary distribution may be unstable or not uniquely defined.")
  }
  
  pi
}

compute_mean_first_passage <- function(P, pi_vec = NULL) {
  if (any(is.na(P))) return(NULL)
  if (is.null(pi_vec)) pi_vec <- compute_stationary_distribution(P)
  if (any(is.na(pi_vec)) || any(pi_vec <= 0)) return(NULL)
  
  K <- nrow(P)
  I <- diag(K)
  One <- matrix(1, nrow = K, ncol = K)
  Pi <- matrix(rep(pi_vec, each = K), nrow = K, byrow = FALSE)
  
  Z <- solve(I - P + One * Pi)
  
  M <- matrix(NA_real_, nrow = K, ncol = K, dimnames = dimnames(P))
  
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i == j) {
        M[i, j] <- 1 / pi_vec[j]
      } else {
        M[i, j] <- (Z[j, j] - Z[i, j]) / pi_vec[j]
      }
    }
  }
  
  M
}

compute_structural_indices <- function(P, pi_vec = NULL, log_base = exp(1)) {
  if (any(is.na(P))) stop("Transition matrix contains NA values.")
  
  K <- nrow(P)
  diag_vals <- diag(P)
  
  if (is.null(pi_vec)) {
    pi_vec <- compute_stationary_distribution(P)
  }
  
  mean_persistence <- mean(diag_vals)
  
  state_volatility <- 1 - diag_vals
  names(state_volatility) <- rownames(P)
  
  pd_bar <- mean(diag_vals)
  differential_stability <- sqrt(mean((diag_vals - pd_bar)^2))
  
  asym_vals <- c()
  if (K > 1) {
    for (i in 1:(K - 1)) {
      for (j in (i + 1):K) {
        asym_vals <- c(asym_vals, abs(P[i, j] - P[j, i]))
      }
    }
  }
  
  reciprocal_asymmetry <- if (length(asym_vals) > 0) {
    (2 / (K * (K - 1))) * sum(asym_vals)
  } else {
    NA_real_
  }
  
  transition_dispersion <- apply(P, 1, function(r) {
    sqrt(mean((r - 1 / K)^2))
  })
  
  state_entropy <- apply(P, 1, function(r) {
    r_nonzero <- r[r > 0]
    -sum(r_nonzero * (log(r_nonzero) / log(log_base)))
  })
  
  entropy_rate <- NA_real_
  if (!any(is.na(pi_vec))) {
    entropy_rate <- sum(pi_vec * state_entropy)
  }
  
  attractor_strength <- colSums(P)
  names(attractor_strength) <- colnames(P)
  
  expected_residence_time <- sapply(diag_vals, function(x) {
    if (isTRUE(all.equal(x, 1))) Inf else 1 / (1 - x)
  })
  names(expected_residence_time) <- rownames(P)
  
  mean_recurrence_time <- rep(NA_real_, K)
  if (!any(is.na(pi_vec)) && all(pi_vec > 0)) {
    mean_recurrence_time <- 1 / pi_vec
    names(mean_recurrence_time) <- names(pi_vec)
  }
  
  mean_first_passage_time <- compute_mean_first_passage(P, pi_vec)
  
  list(
    K = K,
    mean_persistence = mean_persistence,
    state_persistence = setNames(diag_vals, rownames(P)),
    state_volatility = state_volatility,
    differential_stability = differential_stability,
    reciprocal_asymmetry = reciprocal_asymmetry,
    transition_dispersion = transition_dispersion,
    state_entropy = state_entropy,
    entropy_rate = entropy_rate,
    attractor_strength = attractor_strength,
    expected_residence_time = expected_residence_time,
    stationary_distribution = pi_vec,
    mean_recurrence_time = mean_recurrence_time,
    mean_first_passage_time = mean_first_passage_time
  )
}

compute_subspace_retention <- function(P, A) {
  if (length(A) == 0) return(NA_real_)
  mean(rowSums(P[A, A, drop = FALSE]))
}

compute_transition_mass <- function(P, A, B) {
  if (length(A) == 0 || length(B) == 0) return(NA_real_)
  mean(rowSums(P[A, B, drop = FALSE]))
}

compute_distance_weighted_jump <- function(P, distance_matrix) {
  if (is.null(distance_matrix)) return(NA_real_)
  
  if (!all(rownames(P) %in% rownames(distance_matrix)) ||
      !all(colnames(P) %in% colnames(distance_matrix))) {
    stop("Distance matrix row/column names must match transition matrix states.")
  }
  
  D <- distance_matrix[rownames(P), colnames(P), drop = FALSE]
  off_diag <- row(P) != col(P)
  
  sum(P[off_diag] * D[off_diag])
}

compute_theory_driven_indices <- function(P,
                                          subsets = list(),
                                          distance_matrix = NULL) {
  out <- list()
  
  if (!is.null(subsets$A)) {
    out$subspace_retention_A <- compute_subspace_retention(P, subsets$A)
  }
  
  if (!is.null(subsets$R) && !is.null(subsets$C)) {
    out$recovery_index <- compute_transition_mass(P, subsets$R, subsets$C)
  }
  
  if (!is.null(subsets$L) && !is.null(subsets$H)) {
    out$activation_index <- compute_transition_mass(P, subsets$L, subsets$H)
  }
  
  if (!is.null(subsets$N) && !is.null(subsets$P)) {
    out$valence_shift <- compute_transition_mass(P, subsets$N, subsets$P)
  }
  
  if (!is.null(distance_matrix)) {
    out$distance_weighted_jump <- compute_distance_weighted_jump(P, distance_matrix)
  }
  
  out
}

# ============================================================
# Worked example: transition count matrix from the paper
# ============================================================

state_levels <- c(
  "Positive-dominant",
  "Mixed-high",
  "Negative-dominant",
  "Low-affect"
)

count_matrix <- matrix(
  c(
    28, 8, 6, 8,
    7, 12, 7, 4,
    6, 5, 39, 8,
    7, 4, 8, 14
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(state_levels, state_levels)
)

P <- row_normalize(count_matrix)

# ============================================================
# Theory-driven structures used in the worked example
# ============================================================

theory_subsets <- list(
  P = c("Positive-dominant"),
  N = c("Negative-dominant"),
  L = c("Low-affect"),
  H = c("Mixed-high", "Negative-dominant"),
  R = c("Negative-dominant"),
  C = c("Positive-dominant", "Low-affect"),
  A = c("Positive-dominant", "Mixed-high")
)

distance_matrix <- matrix(
  c(
    0, 1, 2, 1,
    1, 0, 1, 2,
    2, 1, 0, 1,
    1, 2, 1, 0
  ),
  nrow = 4,
  byrow = TRUE,
  dimnames = list(state_levels, state_levels)
)

# ============================================================
# Compute indices
# ============================================================

structural <- compute_structural_indices(P)

theory <- compute_theory_driven_indices(
  P,
  subsets = theory_subsets,
  distance_matrix = distance_matrix
)

res_S1 <- list(
  subject_id = "S1",
  count_matrix = count_matrix,
  transition_matrix = P,
  structural_indices = structural,
  theory_driven_indices = theory
)

# ============================================================
# Tidy outputs
# ============================================================

tidy_structural_indices <- function(res) {
  s <- res$structural_indices
  
  data.frame(
    subject = res$subject_id,
    mean_persistence = s$mean_persistence,
    differential_stability = s$differential_stability,
    reciprocal_asymmetry = s$reciprocal_asymmetry,
    entropy_rate = s$entropy_rate,
    stringsAsFactors = FALSE
  )
}

tidy_state_level_indices <- function(res) {
  s <- res$structural_indices
  
  data.frame(
    subject = res$subject_id,
    state = names(s$state_persistence),
    persistence = as.numeric(s$state_persistence),
    volatility = as.numeric(s$state_volatility),
    transition_dispersion = as.numeric(s$transition_dispersion),
    entropy = as.numeric(s$state_entropy),
    attractor_strength = as.numeric(s$attractor_strength),
    expected_residence_time = as.numeric(s$expected_residence_time),
    stringsAsFactors = FALSE
  )
}

tidy_theory_indices <- function(res) {
  vals <- unlist(res$theory_driven_indices)
  
  data.frame(
    subject = res$subject_id,
    index = names(vals),
    value = as.numeric(vals),
    stringsAsFactors = FALSE
  )
}

structural_summary <- tidy_structural_indices(res_S1)
state_summary <- tidy_state_level_indices(res_S1)
theory_summary <- tidy_theory_indices(res_S1)

# ============================================================
# Print main results
# ============================================================

cat("\n--- Count matrix ---\n")
print(res_S1$count_matrix)

cat("\n--- Transition matrix ---\n")
print(round(res_S1$transition_matrix, 3))

cat("\n--- Structural summary ---\n")
print(structural_summary %>%
        mutate(across(where(is.numeric), ~ round(.x, 3))))

cat("\n--- State-level summary ---\n")
print(data.frame(
  subject = state_summary$subject,
  state = state_summary$state,
  persistence = round(state_summary$persistence, 3),
  volatility = round(state_summary$volatility, 3),
  transition_dispersion = round(state_summary$transition_dispersion, 3),
  entropy = round(state_summary$entropy, 3),
  attractor_strength = round(state_summary$attractor_strength, 3),
  expected_residence_time = round(state_summary$expected_residence_time, 3)
))

cat("\n--- Theory-driven summary ---\n")
print(data.frame(
  subject = theory_summary$subject,
  index = theory_summary$index,
  value = round(theory_summary$value, 3)
))

cat("\n--- Stationary distribution ---\n")
print(round(structural$stationary_distribution, 3))

cat("\n--- Mean recurrence time ---\n")
print(round(structural$mean_recurrence_time, 3))

cat("\n--- Mean first passage time ---\n")
print(round(structural$mean_first_passage_time, 3))

# ============================================================
# Optional visualization layer
# ============================================================
# This section generates graphical representations of the
# transition matrix and derived indices. It is not required
# for computing Markov-based features, but provides useful
# visual summaries of the dynamics.

plot_transition_heatmap <- function(P, title = "Transition matrix heatmap") {
  df <- as.data.frame(as.table(P))
  colnames(df) <- c("from", "to", "probability")
  
  ggplot(df, aes(x = to, y = from, fill = probability)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", probability)), size = 4) +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
    labs(
      title = title,
      x = "Next state",
      y = "Current state",
      fill = "P"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

plot_state_indices <- function(state_df,
                               index_col = "persistence",
                               title = NULL) {
  if (is.null(title)) title <- paste("State-level", index_col)
  
  ggplot(state_df, aes(x = state, y = .data[[index_col]])) +
    geom_col() +
    labs(
      title = title,
      x = "State",
      y = index_col
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1)
    )
}

plot_markov_chain <- function(P,
                              threshold = 0.05,
                              title = "Discrete-time Markov chain",
                              state_coordinates = NULL,
                              vertex_size = 35,
                              edge_curved = 0.15) {
  edges <- data.frame()
  
  for (i in seq_len(nrow(P))) {
    for (j in seq_len(ncol(P))) {
      if (!is.na(P[i, j]) && P[i, j] >= threshold) {
        edges <- rbind(
          edges,
          data.frame(
            from = rownames(P)[i],
            to = colnames(P)[j],
            weight = P[i, j]
          )
        )
      }
    }
  }
  
  g <- graph_from_data_frame(
    edges,
    directed = TRUE,
    vertices = data.frame(name = rownames(P))
  )
  
  if (is.null(state_coordinates)) {
    layout_mat <- layout_in_circle(g)
  } else {
    if (!all(V(g)$name %in% rownames(state_coordinates))) {
      stop("State coordinates must have rownames matching state names.")
    }
    layout_mat <- as.matrix(state_coordinates[V(g)$name, c("x", "y")])
  }
  
  E(g)$label <- sprintf("%.2f", E(g)$weight)
  E(g)$width <- 1 + 8 * E(g)$weight
  
  plot(
    g,
    layout = layout_mat,
    main = title,
    vertex.size = vertex_size,
    vertex.label.cex = 0.9,
    vertex.frame.color = "black",
    edge.arrow.size = 0.35,
    edge.curved = edge_curved,
    edge.label = E(g)$label,
    edge.label.cex = 0.75
  )
}

state_coordinates <- data.frame(
  x = c(1, 1, -1, -1),
  y = c(-1, 1, 1, -1),
  row.names = c(
    "Positive-dominant",
    "Mixed-high",
    "Negative-dominant",
    "Low-affect"
  )
)

p_heat <- plot_transition_heatmap(
  P,
  title = "Transition matrix heatmap"
)
print(p_heat)

p_persist <- plot_state_indices(
  state_summary,
  index_col = "persistence",
  title = "State persistence"
)
print(p_persist)

p_entropy <- plot_state_indices(
  state_summary,
  index_col = "entropy",
  title = "State entropy"
)
print(p_entropy)

plot_markov_chain(
  P,
  threshold = 0.05,
  title = "Discrete-time Markov chain",
  state_coordinates = state_coordinates
)

# ============================================================
# Generate plots
# ============================================================

# Transition matrix heatmap
p_heat <- plot_transition_heatmap(
  P,
  title = "Transition probability matrix (heatmap)"
)
print(p_heat)

# State-level persistence
p_persist <- plot_state_indices(
  state_summary,
  index_col = "persistence",
  title = "State-level persistence"
)
print(p_persist)

# State-level entropy
p_entropy <- plot_state_indices(
  state_summary,
  index_col = "entropy",
  title = "State-level entropy"
)
print(p_entropy)

# Markov chain graph
plot_markov_chain(
  P,
  threshold = 0.05,
  title = "Discrete-time Markov chain",
  state_coordinates = state_coordinates
)

# ============================================================
# Optional: save outputs
# ============================================================

# write.csv(structural_summary, "structural_summary.csv", row.names = FALSE)
# write.csv(state_summary, "state_level_summary.csv", row.names = FALSE)
# write.csv(theory_summary, "theory_driven_summary.csv", row.names = FALSE)
# ggsave("transition_heatmap.png", p_heat, width = 7, height = 5, dpi = 300)
# ggsave("state_persistence.png", p_persist, width = 7, height = 5, dpi = 300)
# ggsave("state_entropy.png", p_entropy, width = 7, height = 5, dpi = 300)


# ============================================================
# Optional: multi-subject workflow
# ============================================================

# This section applies the same Markov workflow to multiple subjects.
# It assumes that the input data contain at least three columns:
# subject, time, and state.
#
# Example input object:
# state_sequence_multi
#
# Required structure:
# subject | time | state
#
# The object state_sequence_multi can be obtained from the preprocessing
# script or loaded from a CSV file.

# Uncomment this line if you want to import a preprocessed CSV file:
# state_sequence_multi <- read.csv("state_sequence_multi.csv", stringsAsFactors = FALSE)

compute_subject_markov <- function(data,
                                   subject_id,
                                   state_levels) {
  
  dat_sub <- data |>
    dplyr::filter(subject == subject_id) |>
    dplyr::arrange(time)
  
  if (nrow(dat_sub) < 2) {
    warning("Subject ", subject_id, " has fewer than two valid observations.")
    return(NULL)
  }
  
  transitions_sub <- dat_sub |>
    dplyr::mutate(
      next_state = dplyr::lead(state),
      next_time = dplyr::lead(time)
    ) |>
    dplyr::filter(!is.na(next_state))
  
  if (nrow(transitions_sub) == 0) {
    warning("Subject ", subject_id, " has no valid transitions.")
    return(NULL)
  }
  
  C_sub <- table(
    factor(transitions_sub$state, levels = state_levels),
    factor(transitions_sub$next_state, levels = state_levels)
  )
  
  C_sub <- as.matrix(C_sub)
  
  P_sub <- row_normalize(C_sub)
  
  structural_sub <- compute_structural_indices(P_sub)
  
  theory_sub <- compute_theory_driven_indices(
    P_sub,
    subsets = theory_subsets,
    distance_matrix = distance_matrix
  )
  
  list(
    subject_id = subject_id,
    count_matrix = C_sub,
    transition_matrix = P_sub,
    structural_indices = structural_sub,
    theory_driven_indices = theory_sub
  )
}

analyze_multiple_subjects <- function(data,
                                      state_levels) {
  
  ids <- unique(data$subject)
  
  results <- lapply(
    ids,
    function(id) {
      compute_subject_markov(
        data = data,
        subject_id = id,
        state_levels = state_levels
      )
    }
  )
  
  names(results) <- ids
  results <- results[!sapply(results, is.null)]
  
  results
}

# ------------------------------------------------------------
# Example use
# ------------------------------------------------------------

# Uncomment the following lines when a multi-subject state sequence is available.

# results_all <- analyze_multiple_subjects(
#   data = state_sequence_multi,
#   state_levels = state_levels
# )

# structural_all <- do.call(
#   rbind,
#   lapply(results_all, tidy_structural_indices)
# )

# state_all <- do.call(
#   rbind,
#   lapply(results_all, tidy_state_level_indices)
# )

# theory_all <- do.call(
#   rbind,
#   lapply(results_all, function(x) {
#     y <- tidy_theory_indices(x)
#     if (is.null(y)) return(NULL)
#     y
#   })
# )

# print(structural_all)
# print(state_all)
# print(theory_all)

# Optional export:
# write.csv(structural_all, "multi_subject_structural_indices.csv", row.names = FALSE)
# write.csv(state_all, "multi_subject_state_level_indices.csv", row.names = FALSE)
# write.csv(theory_all, "multi_subject_theory_driven_indices.csv", row.names = FALSE)