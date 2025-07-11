# Test that specifying weights of all 1 reproduces identical net benefit as unweighted analysis,
# ensuring that the weights pipeline does not alter results when weights are uniform.
test_that("weights of 1 reproduce unweighted results", {
  out_unweighted <- dca(cancer ~ cancerpredmarker, data = df_binary)
  out_weighted <- dca(cancer ~ cancerpredmarker, data = df_binary,
                      weights = rep(1, nrow(df_binary)))
  expect_equal(as_tibble(out_weighted)$net_benefit,
               as_tibble(out_unweighted)$net_benefit,
               tolerance = 1e-10)
})

# Test that passing a scalar weight (e.g., 2) behaves identically to passing a vector of that scalar,
# and that results match unweighted analysis up to a scaling factor, confirming consistent weight broadcasting.
test_that("scalar weight recycled correctly", {
  out_scalar <- dca(cancer ~ cancerpredmarker, data = df_binary, weights = 2)
  out_vector <- dca(cancer ~ cancerpredmarker, data = df_binary,
                    weights = rep(2, nrow(df_binary)))
  out_unweighted <- dca(cancer ~ cancerpredmarker, data = df_binary)
  expect_equal(as_tibble(out_scalar)$net_benefit,
               as_tibble(out_vector)$net_benefit,
               tolerance = 1e-10)
  expect_equal(as_tibble(out_scalar)$net_benefit,
               as_tibble(out_unweighted)$net_benefit,
               tolerance = 1e-10)
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

# Test that manual calculation of weighted net benefit matches the output of dca:
# computes net benefit explicitly using the weighted formula,
# comparing it to the dca output to confirm the internal calculation logic is correct under arbitrary weights.
test_that("manual weighted net benefit matches explicit manual calculation", {
  set.seed(123)
  n <- 20
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
