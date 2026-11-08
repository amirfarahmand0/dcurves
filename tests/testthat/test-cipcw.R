# cIPCW (`cipcw=`): ported from voipred, so the net benefit equals voipred's.

# Net benefit from voipred (0.0.2.9000), at thresholds seq(0.05, 0.5, by = 0.05)
voipred_ref <- list(
  cox = list(
    model = c(0.10955405317296737, 0.085702504196406223, 0.077744405785840909, 0.063715262270553075, 0.054884394850305163, 0.053290786752865166, 0.0486335570773576, 0.045745635382867446, 0.037950258367561801, 0.035543341208372221),
    all = c(0.1017372949389598, 0.051833811324457579, -0.0039406703623390549, -0.066686962259985216, -0.13779942641065088, -0.21907081401141168, -0.31284549201228962, -0.42224928301331366, -0.55154467237816018, -0.7066991396159763)),
  untruncated = list(
    model = c(0.10955405317296737, 0.085702504196406223, 0.077744405785840909, 0.063715262270553075, 0.054884394850305163, 0.053290786752865166, 0.0486335570773576, 0.045745635382867446, 0.037950258367561801, 0.035543341208372221),
    all = c(0.1017372949389598, 0.051833811324457579, -0.0039406703623390549, -0.066686962259985216, -0.13779942641065088, -0.21907081401141168, -0.31284549201228962, -0.42224928301331366, -0.55154467237816018, -0.7066991396159763)),
  own_model = list(
    model = c(0.10895152441285083, 0.085263919847311062, 0.077460597297347861, 0.063579492048155717, 0.054872454317016373, 0.053225494584654018, 0.048652569973001845, 0.045720664773655155, 0.038002775111367328, 0.035495540411639999),
    all = c(0.10102442559921722, 0.051081338132507065, -0.004737406683227835, -0.067533494600929556, -0.13870239424099151, -0.22003827954391952, -0.31388737797037491, -0.4233779928012395, -0.55277599214680651, -0.70805359136148738)),
  competing = list(
    model = c(0.10446508115232503, 0.081691985491407978, 0.074932085451985264, 0.061595185499501971, 0.055305693391277722, 0.053149387045219369, 0.048321915917822506, 0.044748369290513089, 0.036856068068763903, 0.034128973588479355),
    all = c(0.096405908147814107, 0.046206236378248222, -0.0098992791289136683, -0.073017984074470743, -0.14455251634610211, -0.22630626751368088, -0.32063751886088715, -0.43069064543262775, -0.56075343138104827, -0.71682877451915328)),
  weights_binary = list(
    model = c(0.11235406620090618, 0.096078370194763035, 0.086950812361683033, 0.073547367174188388, 0.069139582137422662, 0.058467059369117683, 0.056762834610468785, 0.043693916449021215, 0.033308912470370974, 0.026471974480612406),
    all = c(0.10057577662749465, 0.050607764217911008, -0.0052388378869177854, -0.0680662652548501, -0.13927068293850675, -0.22064716029125725, -0.3145430956982771, -0.42408835367313352, -0.55355093127978194, -0.70890602440776007)),
  weights_km = list(
    model = c(0.11273955014230705, 0.085335781789153589, 0.076358917546580973, 0.059342569913617103, 0.054503061365607813, 0.051035321962905042, 0.046818347568609385, 0.043485104779443097, 0.033667032970654653, 0.025667381868723277),
    all = c(0.10512946276624494, 0.055414432919925205, -0.00014942396713804595, -0.06265876296508413, -0.13350268049608974, -0.21446715767438188, -0.30788770826471901, -0.41687835062011225, -0.54568547340375861, -0.70025402074413456)),
  weights_cipcw = list(
    model = c(0.11427866902865197, 0.08841366159245849, 0.076835312565006308, 0.060036144232528191, 0.05181996258501409, 0.04875577225756883, 0.042466831813653641, 0.041224285110240427, 0.033784535721387343, 0.027210134235467406),
    all = c(0.1066575748118824, 0.057027440079209207, 0.0015584659662214728, -0.060844129910889638, -0.13156707190494893, -0.21239329132673102, -0.30565431373647961, -0.41445883988118626, -0.54304600714311213, -0.69735060785742342))
)

z <- seq(0.05, 0.5, by = 0.05)
nb <- function(x, v) x$dca$net_benefit[x$dca$variable == v]
quiet <- function(e) suppressMessages(suppressWarnings(e))
# P(remaining uncensored) from a coxph censoring model, for predict_function=
pf_cox <- function(model, newdata, times) {
  sf <- survival::survfit(model, newdata = newdata)
  vapply(seq_len(nrow(newdata)), function(i) {
    j <- findInterval(times[i], sf$time)
    if (j == 0) 1 else sf$surv[j, i]
  }, numeric(1))
}

test_that("cIPCW gives the net benefit of voipred", {
  cox <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   thresholds = z, cipcw = ~ age + famhistory))
  expect_equal(nb(cox, "cancerpredmarker"), voipred_ref$cox$model, tolerance = 1e-12)
  expect_equal(nb(cox, "all"), voipred_ref$cox$all, tolerance = 1e-12)
  expect_equal(nb(cox, "none"), rep(0, length(z)))

  untr <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                    thresholds = z, cipcw = cipcw_control(~ age + famhistory, min_prob = NULL)))
  expect_equal(nb(untr, "cancerpredmarker"), voipred_ref$untruncated$model, tolerance = 1e-12)

  cf <- survival::coxph(Surv(ttcancer, cancer == 0) ~ age + marker, data = df_surv)
  own <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   thresholds = z, cipcw = cipcw_control(model = cf, predict_function = pf_cox)))
  expect_equal(nb(own, "cancerpredmarker"), voipred_ref$own_model$model, tolerance = 1e-8)
  expect_equal(nb(own, "all"), voipred_ref$own_model$all, tolerance = 1e-8)

  cr <- quiet(dca(Surv(ttcancer, cancer_cr) ~ cancerpredmarker, data = df_surv, time = 1,
                  thresholds = z, cipcw = ~ age + famhistory))
  expect_equal(nb(cr, "cancerpredmarker"), voipred_ref$competing$model, tolerance = 1e-12)
  expect_equal(nb(cr, "all"), voipred_ref$competing$all, tolerance = 1e-12)
})

test_that("the ways of giving cIPCW agree", {
  f <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z, cipcw = ~ age + famhistory))
  ctl <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   thresholds = z, cipcw = cipcw_control(~ age + famhistory)))
  expect_equal(f$dca, ctl$dca)
  # the same Cox model, fitted yourself and given with predict_function
  cf <- survival::coxph(Surv(ttcancer, cancer == 0) ~ age + famhistory, data = df_surv)
  p <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z, cipcw = cipcw_control(model = cf, predict_function = pf_cox)))
  expect_equal(nb(p, "cancerpredmarker"), nb(f, "cancerpredmarker"), tolerance = 1e-8)
  # given weights: those of the Cox model, without further truncation
  y <- survival::Surv(df_surv$ttcancer, df_surv$cancer)
  cw <- quiet(dcurves:::.cipcw_weights(cipcw_control(~ age + famhistory), df_surv,
                                       y[, "time"], y[, "status"], horizon = 1))
  g <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z,
                 cipcw = cipcw_control(weight_event = cw$weight_event,
                                       weight_horizon = cw$weight_horizon, min_prob = NULL)))
  expect_equal(nb(g, "cancerpredmarker"), nb(f, "cancerpredmarker"), tolerance = 1e-12)
})

test_that("cIPCW rates: tp_rate = prevalence x Se, fp_rate = (1 - prevalence)(1 - Sp)", {
  x <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z, cipcw = ~ age + famhistory))
  d <- x$dca
  m <- d[d$variable == "cancerpredmarker", ]
  expect_equal(m$net_benefit, m$tp_rate - m$threshold / (1 - m$threshold) * m$fp_rate)
  a <- d[d$variable == "all", ]
  expect_equal(a$tp_rate, a$pos_rate)                         # everyone treated
  expect_equal(a$fp_rate, 1 - a$pos_rate)
  expect_true(all(m$tp_rate <= m$pos_rate + 1e-12))
  expect_equal(x$cipcw$min_prob, 0.02)
  # test_consequences() gives the same rates
  tc <- quiet(test_consequences(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv,
                                time = 1, thresholds = z, cipcw = ~ age + famhistory,
                                statistics = c("tp_rate", "fp_rate")))
  expect_equal(tc$tp_rate, m$tp_rate)
  expect_equal(tc$fp_rate, m$fp_rate)
  # net interventions avoided and standardized net benefit work on the result
  nia <- net_intervention_avoided(x)
  expect_true(all(is.finite(nia$dca$net_intervention_avoided[nia$dca$variable == "cancerpredmarker"])))
  expect_true(all(is.finite(standardized_net_benefit(x)$dca$standardized_net_benefit)))
})

test_that("truncation", {
  expect_warning(
    dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1, thresholds = z,
        cipcw = cipcw_control(~ age + famhistory, min_prob = 0.9)),
    "truncated at min_prob = 0.9")
  hi <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                  thresholds = z, cipcw = cipcw_control(~ age + famhistory, min_prob = 0.9)))
  lo <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                  thresholds = z, cipcw = cipcw_control(~ age + famhistory, min_prob = NULL)))
  expect_false(isTRUE(all.equal(nb(hi, "cancerpredmarker"), nb(lo, "cancerpredmarker"))))
})

test_that("rows dropped for missing values are left out of the weights too", {
  d <- df_surv
  d$cancerpredmarker[c(3, 10, 50)] <- NA
  x <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = d, time = 1,
                 thresholds = z, cipcw = ~ age + famhistory))
  ref <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = d[-c(3, 10, 50), ],
                   time = 1, thresholds = z, cipcw = ~ age + famhistory))
  expect_equal(x$dca, ref$dca)
  expect_equal(x$n, nrow(d) - 3)
})

test_that("cipcw_control() and cipcw= check their inputs", {
  expect_error(cipcw_control(), "Give the covariates of the censoring")
  expect_error(cipcw_control("age"), "`formula=` must be a formula")
  expect_error(cipcw_control(~ age, min_prob = 2), "`min_prob=` must be a number between 0 and 1")
  expect_error(cipcw_control(weight_event = 1:3), "must be given together")
  expect_error(cipcw_control(weight_event = c(1, -1), weight_horizon = c(1, 1)),
               "must be positive, finite")
  # your own model always needs predict_function, and goes in cipcw_control()
  cf <- survival::coxph(Surv(ttcancer, cancer == 0) ~ age, data = df_surv)
  expect_error(cipcw_control(model = cf), "needs `predict_function=`")
  expect_error(cipcw_control(cf), "A fitted censoring model goes in `model=`")
  expect_error(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   cipcw = cf), "Your own censoring model goes in")
  expect_error(cipcw_control(~ age, predict_function = identity),
               "is for your own censoring `model=`")
  expect_output(print(cipcw_control(~ age + famhistory)), "Cox censoring model on ~age \\+ famhistory; truncation at 0.02")
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary, cipcw = ~ age),
               "`cipcw=` is for survival outcomes")
  expect_error(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   cipcw = ~ age, prevalence = 0.2),
               "`prevalence=` cannot be combined with `cipcw=`")
  expect_error(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   cipcw = ~ age + nope),
               "not in `data=`: nope")
  d <- df_surv
  d$age[5] <- NA
  expect_error(quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = d, time = 1,
                         cipcw = ~ age)), "have missing values \\(age\\)")
  ev <- survival::coxph(Surv(ttcancer, cancer) ~ age, data = df_surv)
  expect_error(quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                         cipcw = cipcw_control(model = ev, predict_function = pf_cox))),
               "fitted to the event, not to censoring")
  expect_error(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                   cipcw = "age"), "`cipcw=` must be")
})
