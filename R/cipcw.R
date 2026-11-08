#' Options for conditional inverse probability of censoring weighting (cIPCW)
#'
#' With a survival outcome, `dca(cipcw = )` estimates the net benefit with
#' conditional inverse probability of censoring weighting (cIPCW) instead of
#' Kaplan-Meier, which allows censoring to depend on covariates. Usually the
#' covariates of the censoring are all that is needed, e.g.
#' `dca(..., cipcw = ~ age + famhistory)`: a Cox model of the censoring times
#' on these covariates is fitted, and the weights are truncated at
#' `min_prob = 0.02` (i.e. weights of at most 50). `cipcw_control()` sets the
#' other options: another truncation, your own censoring model with a function
#' that predicts from it, or the weights themselves.
#'
#' @section cIPCW:
#' With \eqn{T_i} the event time and \eqn{C_i} the censoring time of patient
#' \eqn{i}, \eqn{\tilde T_i = \min(T_i, C_i)} the observed time, \eqn{X_i} the
#' covariates of the censoring model, \eqn{\tau} the time of the assessment
#' (`time=` of [dca()]) and \eqn{\hat G(t \mid X)} the estimated probability of
#' remaining uncensored up to time \eqn{t} given \eqn{X}, each patient gets a
#' weight \eqn{w_i = 1/\hat G(\tilde T_i- \mid X_i)} if an event (of any cause)
#' is observed by \eqn{\tau}, \eqn{w_i = 1/\hat G(\tau \mid X_i)} if the
#' patient is event-free at \eqn{\tau}, and \eqn{w_i = 0} if censored before
#' \eqn{\tau}. The prevalence \eqn{\rho} is the weighted share of patients
#' with the event of interest by \eqn{\tau}; the sensitivity and specificity
#' are the weighted proportions treated (risk at or above the threshold
#' \eqn{z}) among the patients with the event of interest by \eqn{\tau}, and
#' not treated among the others (event-free at \eqn{\tau}, or with a competing
#' event by \eqn{\tau}). The net benefit is
#' \deqn{\rho \, Se - (1 - \rho) (1 - Sp) \frac{z}{1-z},}
#' i.e. `tp_rate` is \eqn{\rho \, Se} and `fp_rate` is \eqn{(1 - \rho)(1 - Sp)}.
#' With case weights (`weights=` of [dca()]), every weighted sum is also
#' weighted by them, and the Cox censoring model is fitted with them.
#'
#' @param formula The covariates of the censoring model, e.g.
#'   `~ age + famhistory`: a Cox model of the censoring times on them is fitted
#'   on `data=` of [dca()] (with its case weights, if any). The same as giving
#'   the formula as `cipcw=` directly.
#' @param model,predict_function Your own censoring model, fitted on the same
#'   patients with the censored patients as the events, together with
#'   `predict_function`, a `function(model, newdata, times)` returning, for each
#'   row of `newdata`, the probability of remaining uncensored at the time of
#'   the same position in `times`. Both are needed, for any kind of model.
#' @param weight_event,weight_horizon Your own censoring weights, one per
#'   patient (row of `data=`), instead of a model: `weight_event` for the
#'   patients with an event by the time of the assessment
#'   (\eqn{1/\hat G(\tilde T_i- \mid X_i)}), and `weight_horizon` for the
#'   patients event-free at it (\eqn{1/\hat G(\tau \mid X_i)}); the other values
#'   are not used. Give both.
#' @param min_prob Truncation of the censoring probabilities at `min_prob`
#'   (equivalently, of the weights at `1 / min_prob`), with a warning saying how
#'   many were truncated. Default 0.02 (weights of at most 50). `NULL`: no
#'   truncation (an error for a censoring probability of 0, and a warning for
#'   probabilities below 0.01).
#' @param ties For a censoring model: `"left"` (default) evaluates the censoring
#'   probability just before each patient's own time, so a censoring step at
#'   exactly that time is not included (the usual convention, as for
#'   Kaplan-Meier: an event tied with a censoring comes first); `"right"`
#'   includes it.
#'
#' @return An object of class `cipcw_control`, for `cipcw=` of [dca()] and
#'   [test_consequences()].
#' @seealso [dca()]
#' @export
#' @examples
#' # the covariates of the censoring: a Cox censoring model is fitted
#' dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
#'     thresholds = seq(0, 0.5, 0.01), cipcw = ~ age + famhistory)
#'
#' # the same, truncated at 0.05 instead of 0.02
#' dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
#'     thresholds = seq(0, 0.5, 0.01),
#'     cipcw = cipcw_control(~ age + famhistory, min_prob = 0.05))
#'
#' # your own censoring model (censored patients as the events), with a function
#' # returning the probabilities of remaining uncensored; no truncation
#' cens_fit <- survival::survreg(Surv(ttcancer, cancer == 0) ~ age + famhistory,
#'                               data = df_surv, dist = "weibull")
#' uncensored <- function(model, newdata, times) {
#'   lp <- predict(model, newdata = newdata, type = "lp")
#'   1 - stats::pweibull(times, shape = 1 / model$scale, scale = exp(lp))
#' }
#' dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1,
#'     thresholds = seq(0, 0.5, 0.01),
#'     cipcw = cipcw_control(model = cens_fit, predict_function = uncensored,
#'                           min_prob = NULL))
cipcw_control <- function(formula = NULL, model = NULL, weight_event = NULL,
                          weight_horizon = NULL, min_prob = 0.02,
                          ties = c("left", "right"), predict_function = NULL) {
  ties <- match.arg(ties)
  if (!is.null(formula) && !inherits(formula, "formula"))
    stop("`formula=` must be a formula of the covariates of the censoring model, ",
         "e.g. ~ age + famhistory. A fitted censoring model goes in `model=`, ",
         "with `predict_function=`.", call. = FALSE)
  if (!is.null(formula)) {
    if (!is.null(model) || !is.null(weight_event) || !is.null(weight_horizon))
      stop("Give one of `formula=`, a censoring `model=` or the weights, not several.",
           call. = FALSE)
    if (!is.null(predict_function))
      stop("`predict_function=` is for your own censoring `model=`. With `formula=`, ",
           "the Cox censoring model is fitted and predicted from by dcurves.",
           call. = FALSE)
  }
  has_w <- !is.null(weight_event) || !is.null(weight_horizon)
  if (!is.null(model) && has_w)
    stop("Give either a censoring `model=` or the weights, not both.", call. = FALSE)
  if (!is.null(predict_function)) {
    if (!is.function(predict_function))
      stop("`predict_function=` must be a function(model, newdata, times).",
           call. = FALSE)
    if (is.null(model))
      stop("`predict_function=` needs the censoring `model=` it predicts from.",
           call. = FALSE)
  }
  if (xor(is.null(weight_event), is.null(weight_horizon)))
    stop("`weight_event=` and `weight_horizon=` must be given together.", call. = FALSE)
  for (nm in c("weight_event", "weight_horizon")) {
    v <- get(nm)
    if (!is.null(v) && (!is.numeric(v) || anyNA(v) || any(!is.finite(v)) || any(v <= 0)))
      stop("`", nm, "=` must be positive, finite numbers without missing values.",
           call. = FALSE)
  }
  if (has_w && length(weight_event) != length(weight_horizon))
    stop("`weight_event=` and `weight_horizon=` must have the same length.", call. = FALSE)
  if (!is.null(model) && (is.numeric(model) || is.data.frame(model)))
    stop("`model=` must be a fitted censoring model.", call. = FALSE)
  if (!is.null(model) && is.null(predict_function))
    stop("Your own censoring `model=` needs `predict_function=`, a function(model, ",
         "newdata, times) returning the probabilities of remaining uncensored.",
         call. = FALSE)
  if (is.null(formula) && is.null(model) && !has_w)
    stop("Give the covariates of the censoring (`formula=`), a censoring `model=`, ",
         "or the weights.", call. = FALSE)
  if (!is.null(min_prob) && (!is.numeric(min_prob) || length(min_prob) != 1 ||
                             is.na(min_prob) || min_prob <= 0 || min_prob >= 1))
    stop("`min_prob=` must be a number between 0 and 1 (e.g. 0.02), or NULL for no ",
         "truncation.", call. = FALSE)
  structure(list(formula = formula, model = model, weight_event = weight_event,
                 weight_horizon = weight_horizon, min_prob = min_prob, ties = ties,
                 predict_function = predict_function),
            class = "cipcw_control")
}

#' @export
print.cipcw_control <- function(x, ...) {
  what <- if (!is.null(x$formula)) {
    paste("Cox censoring model on", paste(deparse(x$formula), collapse = " "))
  } else if (!is.null(x$model)) {
    paste("censoring model of class", class(x$model)[1])
  } else {
    "given censoring weights"
  }
  cat("cIPCW: ", what, "; ",
      if (is.null(x$min_prob)) "no truncation" else paste("truncation at", x$min_prob),
      "\n", sep = "")
  invisible(x)
}

# `cipcw=` of dca(): NULL, the covariates of the censoring (a formula), or
# cipcw_control().
.as_cipcw <- function(x) {
  if (is.null(x) || inherits(x, "cipcw_control")) return(x)
  if (inherits(x, "formula")) return(cipcw_control(formula = x))
  stop("`cipcw=` must be the covariates of the censoring (e.g. ~ age + famhistory) or ",
       "cipcw_control(). Your own censoring model goes in ",
       "cipcw_control(model = , predict_function = ).", call. = FALSE)
}

# The censoring weights of every patient: list(weight_event, weight_horizon).
# `data` holds the rows used (aligned with `time` and `status`), `w` the case
# weights (or NULL).
.cipcw_weights <- function(cipcw, data, time, status, horizon, w = NULL) {
  n <- length(time)
  needs_event <- status != 0 & time <= horizon
  needs_horizon <- time > horizon

  if (!is.null(cipcw$weight_event)) {
    if (length(cipcw$weight_event) != n)
      stop("The cIPCW weights have length ", length(cipcw$weight_event),
           " but there are ", n, " patients.", call. = FALSE)
    we <- cipcw$weight_event
    wh <- cipcw$weight_horizon
    used <- c(we[needs_event], wh[needs_horizon])
    if (!is.null(cipcw$min_prob)) {
      cap <- 1 / cipcw$min_prob
      n_trunc <- sum(used > cap)
      if (n_trunc > 0)
        warning(n_trunc, " cIPCW weight", if (n_trunc > 1) "s were" else " was",
                " truncated at 1 / min_prob = ", signif(cap, 4), ".", call. = FALSE)
      we <- pmin(we, cap)
      wh <- pmin(wh, cap)
    } else if (any(used > 100)) {
      warning(sum(used > 100), " cIPCW weight", if (sum(used > 100) > 1) "s are" else " is",
              " above 100 (censoring probability below 0.01), a sign of near-",
              "violations of positivity. Consider truncation with `min_prob=`.",
              call. = FALSE)
    }
    return(list(weight_event = we, weight_horizon = wh))
  }

  model <- cipcw$model
  if (!is.null(cipcw$formula)) {
    model <- .fit_censoring_cox(cipcw$formula, data, time, status, w)
  } else {
    .check_censoring_model(model, status)
  }
  g_at <- function(rows, times) {
    nd <- data[rows, , drop = FALSE]
    # the Cox model fitted from `formula=`; a given model has predict_function
    if (is.null(cipcw$predict_function)) {
      return(.cox_survival(model, nd, times, cipcw$ties))
    }
    eps <- if (cipcw$ties == "left") 1e-9 * max(1, max(times)) else 0
    g <- cipcw$predict_function(model, nd, times - eps)
    if (is.matrix(g)) g <- g[cbind(seq_len(nrow(nd)), seq_len(nrow(nd)))]
    if (length(g) != nrow(nd))
      stop("`predict_function=` must return one probability per row of `newdata`.",
           call. = FALSE)
    as.numeric(g)
  }
  g_event <- rep(1, n)
  g_horizon <- rep(1, n)
  if (any(needs_event)) g_event[needs_event] <- g_at(which(needs_event), time[needs_event])
  if (any(needs_horizon))
    g_horizon[needs_horizon] <- g_at(which(needs_horizon), rep(horizon, sum(needs_horizon)))
  used <- c(g_event[needs_event], g_horizon[needs_horizon])
  if (anyNA(used) || any(used < 0 | used > 1))
    stop("The censoring model returned invalid censoring probabilities ",
         "(missing, or outside [0, 1]).", call. = FALSE)
  if (!is.null(cipcw$min_prob)) {
    n_trunc <- sum(used < cipcw$min_prob)
    if (n_trunc > 0)
      warning(n_trunc, " censoring probabilit", if (n_trunc > 1) "ies were" else "y was",
              " truncated at min_prob = ", cipcw$min_prob, ". This signals ",
              "near-violations of positivity.", call. = FALSE)
    g_event <- pmax(g_event, cipcw$min_prob)
    g_horizon <- pmax(g_horizon, cipcw$min_prob)
  } else {
    if (any(used == 0))
      stop(sum(used == 0), " censoring probabilit", if (sum(used == 0) > 1) "ies are" else "y is",
           " exactly 0, which gives an infinite weight (positivity is violated). ",
           "Set `min_prob=` to truncate them.", call. = FALSE)
    if (any(used < 0.01))
      warning(sum(used < 0.01), " censoring probabilit",
              if (sum(used < 0.01) > 1) "ies are" else "y is",
              " below 0.01, a sign of near-violations of positivity. ",
              "Consider truncation with `min_prob=`.", call. = FALSE)
  }
  list(weight_event = 1 / g_event, weight_horizon = 1 / g_horizon)
}

# The Cox model of the censoring times on the covariates of `formula`.
.fit_censoring_cox <- function(formula, data, time, status, w = NULL) {
  rhs <- formula[[length(formula)]]
  vars <- all.vars(rhs)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop("The censoring formula (`cipcw=`) uses ",
         if (length(miss) == 1) "a variable that is" else "variables that are",
         " not in `data=`: ", paste(miss, collapse = ", "), ".", call. = FALSE)
  na <- vars[vapply(vars, function(v) anyNA(data[[v]]), logical(1))]
  if (length(na))
    stop("The covariates of the censoring model have missing values (",
         paste(na, collapse = ", "), "), but every patient needs a censoring ",
         "weight. Impute them, or leave them out of `cipcw=`.", call. = FALSE)
  if (any(c(".dcurves_time", ".dcurves_censored") %in% names(data)))
    stop("`data=` has a column named .dcurves_time or .dcurves_censored, which ",
         "dcurves uses for the censoring model; rename it.", call. = FALSE)
  d <- data
  d$.dcurves_time <- time
  d$.dcurves_censored <- as.integer(status == 0)
  f <- stats::as.formula(call("~", quote(survival::Surv(.dcurves_time, .dcurves_censored)),
                              rhs), env = environment(formula))
  env <- new.env(parent = if (is.environment(environment(formula))) environment(formula)
                          else globalenv())
  env$.dcurves_data <- d
  cl <- as.call(c(list(quote(survival::coxph), formula = f, data = quote(.dcurves_data)),
                  if (!is.null(w)) list(weights = w)))
  tryCatch(eval(cl, env), error = function(e)
    stop("Could not fit the Cox censoring model ", paste(deparse(formula), collapse = " "),
         " (", conditionMessage(e), "). Give your own censoring model with ",
         "cipcw_control(model = , predict_function = ).", call. = FALSE))
}

# A given censoring model must be fitted to the censoring times: its events are
# the censored patients (checked for models that keep their outcome).
.check_censoring_model <- function(model, status) {
  y <- tryCatch(model[["y"]], error = function(e) NULL)
  if (!inherits(y, "Surv") || nrow(y) != length(status)) return(invisible(NULL))
  ev <- as.integer(y[, ncol(y)] != 0)
  refit <- " Refit it with only censoring as the event, e.g. Surv(time, status == 0)."
  if (identical(ev, as.integer(status == 0))) return(invisible(NULL))
  if (identical(ev, as.integer(status == 1)))
    stop("The censoring model appears to have been fitted to the event, not to ",
         "censoring.", refit, call. = FALSE)
  if (identical(ev, as.integer(status != 0)))
    stop("The censoring model appears to have been fitted to all events, not to ",
         "censoring.", refit, call. = FALSE)
  invisible(NULL)
}

# P(remaining uncensored up to `times`) for each row of `newdata`, from a coxph
# censoring model: just before each time (`side = "left"`) or at it.
.cox_survival <- function(model, newdata, times, side = "left") {
  if (.cox_simple(model)) {
    bh <- .cox_basehaz(model)
    lp <- stats::predict(model, newdata = newdata, type = "lp")
    i <- if (side == "left") findInterval(times, bh$time, left.open = TRUE)
         else findInterval(times, bh$time)
    h <- ifelse(i == 0, 0, bh$hazard[pmax(i, 1)])
    return(unname(exp(-h * exp(lp))))
  }
  # stratified or time-dependent: survfit, which needs the model's data
  sf <- tryCatch(survival::survfit(model, newdata = newdata),
                 error = function(e)
                   stop("Could not compute survival from this coxph censoring model: ",
                        conditionMessage(e), ".", call. = FALSE))
  vapply(seq_len(nrow(newdata)), function(r) {
    sr <- sf[r]
    i <- findInterval(times[r], sr$time, left.open = side == "left")
    if (i == 0) 1 else sr$surv[i]
  }, numeric(1))
}

.cox_simple <- function(model) {
  !is.null(model$y) && is.null(model$strata) &&
    is.null(attr(stats::terms(model), "specials")$strata) &&
    is.null(attr(stats::terms(model), "specials")$tt) &&
    !is.null(model$linear.predictors) && ncol(model$y) == 2
}

# The cumulative baseline hazard of a coxph model (Breslow or Efron), at the
# centred linear predictor, as used by predict(type = "lp").
.cox_basehaz <- function(model) {
  y <- model$y
  tt <- y[, 1]
  ev <- y[, 2]
  w <- if (is.null(model$weights)) rep(1, length(tt)) else model$weights
  r <- w * exp(model$linear.predictors)
  ut <- sort(unique(tt[ev == 1]))
  efron <- identical(model$method, "efron")
  o <- order(tt, decreasing = TRUE)
  risk_sum <- cumsum(r[o])
  at_risk <- function(t) risk_sum[sum(tt[o] >= t)]
  inc <- vapply(ut, function(t) {
    dd <- ev == 1 & tt == t
    R <- at_risk(t)
    if (!efron) return(sum(w[dd]) / R)
    D <- sum(r[dd])
    m <- sum(dd)
    sum((sum(w[dd]) / m) / (R - (0:(m - 1)) / m * D))
  }, numeric(1))
  data.frame(hazard = cumsum(inc), time = ut)
}

# The rates of .calculate_test_consequences() with cIPCW, at each threshold:
# pos_rate (the prevalence), test_pos_rate, tp_rate (prevalence x sensitivity)
# and fp_rate ((1 - prevalence) x (1 - specificity)). `k` is the status of the
# event of interest; other events by `horizon` are competing events.
.cipcw_rates <- function(time, status, risk, thresholds, horizon, cw, w = NULL, k = 1) {
  if (is.null(w)) w <- rep(1, length(time))
  event_k <- status == k & time <= horizon
  competing <- status != 0 & status != k & time <= horizon
  event_free <- time > horizon
  we <- ifelse(event_k | competing, cw$weight_event, 0)
  wh <- ifelse(event_free, cw$weight_horizon, 0)
  pos_w <- w * event_k * we                        # weighted events of interest
  neg_w <- w * (event_free * wh + competing * we)  # weighted non-events
  den_se <- sum(pos_w)
  den_sp <- sum(neg_w)
  prev <- den_se / sum(w)
  treated <- function(z) risk >= z
  tp_w <- vapply(thresholds, function(z) sum(pos_w[treated(z)]), numeric(1))
  fp_w <- vapply(thresholds, function(z) sum(neg_w[treated(z)]), numeric(1))
  tibble::tibble(
    threshold = thresholds,
    n = length(time),
    pos_rate = prev,
    test_pos_rate = vapply(thresholds, function(z) sum(w[treated(z)]) / sum(w), numeric(1)),
    tp_rate = tp_w / sum(w),
    fp_rate = if (den_sp == 0) 0 else (1 - prev) * fp_w / den_sp
  )
}
