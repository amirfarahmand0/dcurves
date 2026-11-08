# Case weights (`weights=`): weighted proportions for a binary outcome, weighted
# Kaplan-Meier (Aalen-Johansen) for a survival outcome; equal to voipred's.

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
set.seed(1)
wb <- round(runif(nrow(df_binary), 0.2, 3), 3)
set.seed(2)
ws <- round(runif(nrow(df_surv), 0.2, 3), 3)

test_that("weights give the net benefit of voipred", {
  b <- quiet(dca(cancer ~ cancerpredmarker, data = df_binary, thresholds = z, weights = wb))
  expect_equal(nb(b, "cancerpredmarker"), voipred_ref$weights_binary$model, tolerance = 1e-12)
  expect_equal(nb(b, "all"), voipred_ref$weights_binary$all, tolerance = 1e-12)
  s <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z, weights = ws))
  expect_equal(nb(s, "cancerpredmarker"), voipred_ref$weights_km$model, tolerance = 1e-12)
  expect_equal(nb(s, "all"), voipred_ref$weights_km$all, tolerance = 1e-12)
  sc <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                  thresholds = z, weights = ws, cipcw = ~ age + famhistory))
  expect_equal(nb(sc, "cancerpredmarker"), voipred_ref$weights_cipcw$model, tolerance = 1e-12)
  expect_equal(nb(sc, "all"), voipred_ref$weights_cipcw$all, tolerance = 1e-12)
  expect_true(b$weights)
})

test_that("equal weights change nothing", {
  for (v in list(1, 2.5)) {
    b0 <- quiet(dca(cancer ~ cancerpredmarker + marker, data = df_binary, thresholds = z,
                    as_probability = "marker"))
    b1 <- quiet(dca(cancer ~ cancerpredmarker + marker, data = df_binary, thresholds = z,
                    as_probability = "marker", weights = rep(v, nrow(df_binary))))
    # (the unweighted binary prevalence keeps a name from table(); the values are equal)
    expect_equal(b1$dca, b0$dca, tolerance = 1e-12, ignore_attr = TRUE)
    s0 <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                    thresholds = z))
    s1 <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                    thresholds = z, weights = rep(v, nrow(df_surv))))
    expect_equal(s1$dca, s0$dca, tolerance = 1e-12)
    c0 <- quiet(dca(Surv(ttcancer, cancer_cr) ~ cancerpredmarker, data = df_surv, time = 1,
                    thresholds = z))
    c1 <- quiet(dca(Surv(ttcancer, cancer_cr) ~ cancerpredmarker, data = df_surv, time = 1,
                    thresholds = z, weights = rep(v, nrow(df_surv))))
    expect_equal(c1$dca, c0$dca, tolerance = 1e-12)
  }
})

test_that("whole-number weights equal duplicated patients", {
  set.seed(3)
  k <- sample(0:3, nrow(df_binary), replace = TRUE)
  b <- quiet(dca(cancer ~ cancerpredmarker, data = df_binary, thresholds = z, weights = k))
  bd <- quiet(dca(cancer ~ cancerpredmarker, data = df_binary[rep(seq_len(nrow(df_binary)), k), ],
                  thresholds = z))
  expect_equal(b$dca$net_benefit, bd$dca$net_benefit, tolerance = 1e-12)
  expect_equal(b$dca$test_pos_rate, bd$dca$test_pos_rate, tolerance = 1e-12)
  k <- sample(0:3, nrow(df_surv), replace = TRUE)
  s <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
                 thresholds = z, weights = k))
  sd <- quiet(dca(Surv(ttcancer, cancer) ~ cancerpredmarker,
                  data = df_surv[rep(seq_len(nrow(df_surv)), k), ], time = 1, thresholds = z))
  expect_equal(s$dca$net_benefit, sd$dca$net_benefit, tolerance = 1e-12)
  cr <- quiet(dca(Surv(ttcancer, cancer_cr) ~ cancerpredmarker, data = df_surv, time = 1,
                  thresholds = z, weights = k))
  crd <- quiet(dca(Surv(ttcancer, cancer_cr) ~ cancerpredmarker,
                   data = df_surv[rep(seq_len(nrow(df_surv)), k), ], time = 1, thresholds = z))
  expect_equal(cr$dca$net_benefit, crd$dca$net_benefit, tolerance = 1e-12)
})

test_that("weights as a column, with prevalence, and in test_consequences()", {
  d <- df_binary
  d$w <- wb
  a <- quiet(dca(cancer ~ cancerpredmarker, data = d, thresholds = z, weights = "w"))
  b <- quiet(dca(cancer ~ cancerpredmarker, data = d, thresholds = z, weights = wb))
  expect_equal(a$dca, b$dca)
  # case-control: the given prevalence is used, with weighted sensitivity and specificity
  p <- quiet(dca(casecontrol ~ cancerpredmarker, data = df_case_control, thresholds = z,
                 prevalence = 0.2, weights = rep(1, nrow(df_case_control))))
  p0 <- quiet(dca(casecontrol ~ cancerpredmarker, data = df_case_control, thresholds = z,
                  prevalence = 0.2))
  expect_equal(p$dca$net_benefit, p0$dca$net_benefit, tolerance = 1e-12)
  tc <- quiet(test_consequences(cancer ~ cancerpredmarker, data = df_binary, thresholds = z,
                                weights = wb, statistics = c("tp_rate", "fp_rate")))
  m <- b$dca[b$dca$variable == "cancerpredmarker", ]
  expect_equal(tc$tp_rate, m$tp_rate)
  expect_equal(tc$fp_rate, m$fp_rate)
})

test_that("rows dropped for missing values are left out of the weights", {
  d <- df_binary
  d$cancerpredmarker[c(2, 7)] <- NA
  x <- quiet(dca(cancer ~ cancerpredmarker, data = d, thresholds = z, weights = wb))
  ref <- quiet(dca(cancer ~ cancerpredmarker, data = d[-c(2, 7), ], thresholds = z,
                   weights = wb[-c(2, 7)]))
  expect_equal(x$dca, ref$dca)
})

test_that("weights are checked", {
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary, weights = "nope"),
               "names a column \\(\"nope\"\\) that is not in `data=`")
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary, weights = 1:3),
               "one weight per row of `data=` \\(750\\)")
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary,
                   weights = c(-1, rep(1, 749))), "non-negative, finite")
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary,
                   weights = c(NA, rep(1, 749))), "non-negative, finite")
  expect_error(dca(cancer ~ cancerpredmarker, data = df_binary, weights = rep(0, 750)),
               "positive total")
})
