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

