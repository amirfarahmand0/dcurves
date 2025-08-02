library(broom)
library(tibble)
library(dplyr)
library(survival)


# Test that specifying weights of all 1 reproduces identical net benefit as unweighted analysis,
# ensuring that the weights pipeline does not alter results when weights are uniform.


test_that("weights of 1 reproduce unweighted results", {
  out_unweighted <- dca(cancer ~ cancerpredmarker, data = df_binary)

  expect_warning(
    out_weighted <- dca(
      cancer ~ cancerpredmarker,
      data = df_binary,
      weights = rep(1, nrow(df_binary))
    ),
    "All weights are equal; results are equivalent to unweighted analysis."
  )

  expect_equal(
    as_tibble(out_weighted)$net_benefit,
    as_tibble(out_unweighted)$net_benefit,
    tolerance = 1e-10
  )
})


test_that("scalar weight recycled correctly", {
  expect_warning(
    out_scalar <- dca(
      cancer ~ cancerpredmarker,
      data = df_binary,
      weights = 2
    ),
    "All weights are equal; results are equivalent to unweighted analysis."
  )

  expect_warning(
    out_vector <- dca(
      cancer ~ cancerpredmarker,
      data = df_binary,
      weights = rep(2, nrow(df_binary))
    ),
    "All weights are equal; results are equivalent to unweighted analysis."
  )

  out_unweighted <- dca(cancer ~ cancerpredmarker, data = df_binary)

  expect_equal(
    as_tibble(out_scalar)$net_benefit,
    as_tibble(out_vector)$net_benefit,
    tolerance = 1e-10
  )

  expect_equal(
    as_tibble(out_scalar)$net_benefit,
    as_tibble(out_unweighted)$net_benefit,
    tolerance = 1e-10
  )
})


# Test that invalid weight specifications are properly rejected:
# - Negative weights should trigger an error.
# - Weight vectors of incorrect length should trigger an error.
# This ensures robust and informative input validation for users.
test_that("invalid weights raise errors", {
  expect_error(
    dca(cancer ~ cancerpredmarker, data = df_binary,
        weights = rep(-1, nrow(df_binary))),
    "non-negative"
  )
  expect_error(
    dca(cancer ~ cancerpredmarker, data = df_binary,
        weights = rep(1, 5)),
    "length"
  )
})


# Test that weighted DCA for binary outcomes is equivalent to explicit replication:
# doubling the weight of the second half of the dataset should produce identical net benefit
# to duplicating the second half of the dataset explicitly and running unweighted analysis.
test_that("weighted binary DCA matches explicit replication approach", {
  n <- nrow(df_binary)
  weights <- rep(1, n)
  weights[(n/2 + 1):n] <- 2
  result_weighted <- dca(
    cancer ~ cancerpredmarker + famhistory,
    data = df_binary,
    thresholds = seq(0, 0.35, by = 0.05),
    weights = weights
  )
  nb_weighted <- as_tibble(result_weighted)$net_benefit
  df_replicated <- dplyr::bind_rows(
    df_binary[1:(n/2), ],
    df_binary[(n/2 + 1):n, ],
    df_binary[(n/2 + 1):n, ]
  )
  result_replicated <- dca(
    cancer ~ cancerpredmarker + famhistory,
    data = df_replicated,
    thresholds = seq(0, 0.35, by = 0.05)
  )
  nb_replicated <- as_tibble(result_replicated)$net_benefit
  expect_equal(nb_weighted, nb_replicated, tolerance = 1e-6)
})

# Test that manual calculation of weighted net benefit matches the output of dca:
# computes net benefit explicitly using the weighted formula,
# comparing it to the dca output to confirm the internal calculation logic is correct under arbitrary weights.
test_that("manual weighted net benefit matches explicit manual calculation", {
  set.seed(123)
  n <- 200
  toy <- data.frame(
    y = rbinom(n, 1, 0.4),
    risk = runif(n),
    w = sample(c(1, 2), n, replace = TRUE)
  )
  threshold <- 0.1
  result <- dca(y ~ risk, data = toy, thresholds = threshold, weights = toy$w)
  nb_df <- as_tibble(result)
  nb_model <- nb_df$net_benefit[nb_df$variable == "risk"]
  w_total <- sum(toy$w)
  predicted_positive <- toy$risk >= threshold
  TP <- sum(toy$w[predicted_positive & toy$y == 1])
  FP <- sum(toy$w[predicted_positive & toy$y == 0])
  p_threshold <- threshold / (1 - threshold)
  manual_nb <- TP / w_total - p_threshold * FP / w_total
  expect_equal(nb_model, manual_nb, tolerance = 1e-6)
})

# Test that weighted DCA for survival outcomes is equivalent to explicit replication:
# doubling the weight of the second half of the dataset should produce identical net benefit
# to duplicating the second half explicitly and running unweighted analysis,
# validating correctness of the weighted survival implementation.
test_that("weighted survival DCA matches explicit replication approach", {
  n <- nrow(df_surv)
  weights <- rep(1, n)
  weights[(n/2 + 1):n] <- 2
  result_weighted <- dca(
    Surv(ttcancer, cancer) ~ cancerpredmarker,
    data = df_surv,
    time = 1,
    thresholds = seq(0, 0.50, by = 0.05),
    weights = weights
  )
  nb_weighted <- as_tibble(result_weighted)$net_benefit
  df_replicated <- dplyr::bind_rows(
    df_surv[1:(n/2), ],
    df_surv[(n/2 + 1):n, ],
    df_surv[(n/2 + 1):n, ]
  )
  result_replicated <- dca(
    Surv(ttcancer, cancer) ~ cancerpredmarker,
    data = df_replicated,
    time = 1,
    thresholds = seq(0, 0.50, by = 0.05)
  )
  nb_replicated <- as_tibble(result_replicated)$net_benefit
  expect_equal(nb_weighted, nb_replicated, tolerance = 1e-6)
})

# Test that manual calculation of weighted survival net benefit matches the output of dca:
# This test generates a small synthetic survival dataset with arbitrary observation weights,
# runs the weighted dca() function at a specified time point and threshold,
# and manually computes the time-dependent weighted net benefit using the same
# survival-Kaplan-Meier-based formula that dca() uses internally for survival outcomes.
# The manually computed net benefit is then compared to the dca() output to confirm
# that the internal calculation logic produces correct results under arbitrary weights.

test_that("manual weighted survival net benefit matches dca output", {
  set.seed(123)
  n <- 200
  toy <- data.frame(
    time = rexp(n, 0.1),
    status = sample(c(0, 1), n, replace = TRUE),
    risk = runif(n),
    w = sample(c(1, 2), n, replace = TRUE)
  )

  t0 <- 3
  threshold <- 0.4
  w_total <- sum(toy$w)
  result <- dca(Surv(time, status) ~ risk,
                data = toy,
                time = t0,
                thresholds = threshold,
                weights = toy$w)

  nb_df <- as_tibble(result)
  nb_model <- nb_df$net_benefit[nb_df$variable == "risk"]
  predicted_positive <- toy$risk >= threshold
  sf_pos <- survfit(Surv(time, status) ~ 1,
                    data = toy[predicted_positive, ],
                    weights = toy$w[predicted_positive])
  S_pos_t0 <- summary(sf_pos, times = t0, extend = TRUE)$surv
  risk_rate_among_pos <- 1 - S_pos_t0

  test_pos_rate <- sum(toy$w[predicted_positive]) / w_total
  tp_rate <- risk_rate_among_pos * test_pos_rate
  fp_rate <- (1 - risk_rate_among_pos) * test_pos_rate

  p_threshold <- threshold / (1 - threshold)
  manual_nb <- tp_rate - p_threshold * fp_rate

  expect_equal(nb_model, manual_nb, tolerance = 1e-6)
})

# Test that weighted competing-risk DCA matches explicit row-replication:
# This test runs dca() on a competing-risk dataset using observation weights,
# then creates an equivalent unweighted dataset by replicating rows in proportion
# to their weights. It compares the net benefit results from the weighted and
# replicated analyses to confirm they are identical within tolerance.

test_that("weighted Competing-Risk DCA matches explicit replication approach", {
  n <- nrow(df_surv)
  weights <- rep(1, n)
  weights[(n/2 + 1):n] <- 2
  result_weighted <- dca(
    Surv(ttcancer, cancer_cr) ~ cancerpredmarker,
    data = df_surv,
    time = 1,
    thresholds = seq(0, 0.50, by = 0.05),
    weights = weights
  )
  nb_weighted <- as_tibble(result_weighted)$net_benefit
  df_replicated <- dplyr::bind_rows(
    df_surv[1:(n/2), ],
    df_surv[(n/2 + 1):n, ],
    df_surv[(n/2 + 1):n, ]
  )
  result_replicated <- dca(
    Surv(ttcancer, cancer_cr) ~ cancerpredmarker,
    data = df_replicated,
    time = 1,
    thresholds = seq(0, 0.50, by = 0.05)
  )
  nb_replicated <- as_tibble(result_replicated)$net_benefit
  expect_equal(nb_weighted, nb_replicated, tolerance = 1e-6)
})

# Test that manual weighted competing-risk net benefit matches dca output:
# This test creates a toy competing-risk dataset with "censor" as the lowest factor level
# and observation weights. It runs weighted dca() at a fixed time and threshold, then
# manually computes the net benefit. For the manual calculation, the cumulative incidence
# function (CIF) for the event of interest is extracted from survfit() output using
# broom::tidy(), filtering for the correct state and time <= t0, and taking the last value.
# The CIF is used to compute TP and FP rates, and the resulting manual NB is compared
# to the dca() output for equality within tolerance.

test_that("manual weighted competing-risk net benefit matches dca output", {
  set.seed(123)
  n <- 200
  toy <- data.frame(
    time = rexp(n, 0.1),
    status = sample(
      c("censor", "diagnosed with cancer", "dead other causes"),
      n,
      replace = TRUE,
      prob = c(0.6, 0.25, 0.15)
    ),
    risk = runif(n),
    w = sample(c(1, 2), n, replace = TRUE)
  )

  toy$status <- factor(
    toy$status,
    levels = c("censor", "diagnosed with cancer", "dead other causes")
  )

  t0 <- 3
  threshold <- 0.4
  w_total <- sum(toy$w)

  result <- dca(Surv(time, status) ~ risk,
                data = toy,
                time = t0,
                thresholds = threshold,
                weights = toy$w)

  nb_df <- as_tibble(result)
  nb_model <- nb_df$net_benefit[nb_df$variable == "risk"]

  predicted_positive <- toy$risk >= threshold

  if (any(predicted_positive)) {
    sf_pos <- survfit(Surv(time, status) ~ 1,
                      data = toy[predicted_positive, ],
                      weights = toy$w[predicted_positive])

    sf_tidy <- broom::tidy(sf_pos)

    cif_event <- sf_tidy %>%
      filter(state == "diagnosed with cancer", time <= t0) %>%
      arrange(time) %>%
      slice_tail(n = 1) %>%
      pull(estimate)

    risk_rate_among_pos <- cif_event
  } else {
    risk_rate_among_pos <- 0
  }

  test_pos_rate <- sum(toy$w[predicted_positive]) / w_total
  tp_rate <- risk_rate_among_pos * test_pos_rate
  fp_rate <- (1 - risk_rate_among_pos) * test_pos_rate

  p_threshold <- threshold / (1 - threshold)
  manual_nb <- tp_rate - p_threshold * fp_rate

  expect_equal(nb_model, manual_nb, tolerance = 1e-6)
})

