test_that("weights of 1 reproduce unweighted results", {
  out_unweighted <- dca(cancer ~ cancerpredmarker, data = df_binary)
  out_weighted <- dca(cancer ~ cancerpredmarker, data = df_binary,
                      weights = rep(1, nrow(df_binary)))
  expect_equal(as_tibble(out_weighted)$net_benefit,
               as_tibble(out_unweighted)$net_benefit,
               tolerance = 1e-10)
})

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

test_that("weighted calculations on toy data", {
  toy <- data.frame(
    y = c(1, 0, 1, 0),
    risk = c(0.9, 0.8, 0.3, 0.2),
    w = c(1, 0.5, 2, 0.5)
  )
  out <- dca(y ~ risk, data = toy, thresholds = 0.5, weights = toy$w)
  nb <- as_tibble(out)$net_benefit
  expect_true(is.numeric(nb))
  expect_length(nb, length(out$dca$threshold))
})

test_that("invalid weights raise errors", {
  expect_error(
    dca(cancer ~ cancerpredmarker, data = df_binary,
        weights = rep(-1, nrow(df_binary))),
    "non-negative"
  )
  expect_error(
    dca(cancer ~ cancerpredmarker, data = df_binary,
        weights = rep(1, 5)),  # wrong length
    "length"
  )
})

test_that("case-control data with prevalence and weights runs correctly", {
  set.seed(42)
  n_cases <- 100
  n_controls <- 100
  df_case_control <- data.frame(
    casecontrol = c(rep(1, n_cases), rep(0, n_controls)),
    cancerpredmarker = c(
      rbeta(n_cases, 2, 1),
      rbeta(n_controls, 1, 2)
    )
  )

  true_prevalence <- 0.15

  # For cases: weight = 1
  # For controls: weight = (p / (1-p)) * (n_cases / n_controls)
  weight_case <- 1
  weight_control <- (true_prevalence / (1 - true_prevalence)) * (n_cases / n_controls)
  weights <- c(rep(weight_case, n_cases), rep(weight_control, n_controls))

  result <- dca(casecontrol ~ cancerpredmarker,
                data = df_case_control,
                prevalence = true_prevalence,
                weights = weights)

  nb <- as_tibble(result)$net_benefit
  expect_true(is.numeric(nb))
  expect_length(nb, length(result$dca$threshold))
})

test_that("manual weighted net benefit matches expected value", {
  toy <- data.frame(
    y = c(1, 0, 1, 0),
    risk = c(0.9, 0.8, 0.3, 0.2),
    w = c(1.0, 0.5, 2.0, 0.5)
  )
  result <- dca(y ~ risk, data = toy, thresholds = 0.5, weights = toy$w)
  nb_df <- as_tibble(result)

  nb_model <- nb_df$net_benefit[nb_df$variable == "risk"]

  expect_equal(nb_model, 0.125, tolerance = 1e-4)
})





