# Case weights (`weights=`) and cIPCW (`cipcw=`) of dca() and test_consequences():
# checked and aligned with the rows of the model frame (rows dropped for missing
# values are left out). Returns list(weights, cipcw_weights, cipcw_info), with
# NULL for what is not used.
.weights_and_cipcw <- function(model_frame, data, outcome_name, outcome_type, weights,
                               cipcw, time, prevalence) {
  n_data <- nrow(data)
  dropped <- attr(model_frame, "na.action")
  kept <- if (is.null(dropped)) seq_len(n_data) else seq_len(n_data)[-dropped]

  # case weights
  w <- NULL
  if (!is.null(weights)) {
    if (is.character(weights) && length(weights) == 1L) {
      if (!weights %in% names(data))
        stop("`weights=` names a column (\"", weights, "\") that is not in `data=`.",
             call. = FALSE)
      weights <- data[[weights]]
    }
    if (!is.numeric(weights) || length(weights) != n_data)
      stop("`weights=` must be the name of a column of `data=`, or numbers with one ",
           "weight per row of `data=` (", n_data, ").", call. = FALSE)
    w <- as.numeric(weights)[kept]
    if (anyNA(w) || any(!is.finite(w)) || any(w < 0))
      stop("`weights=` must be non-negative, finite numbers without missing values.",
           call. = FALSE)
    if (sum(w) <= 0)
      stop("`weights=` must have a positive total.", call. = FALSE)
  }

  # cIPCW
  cipcw <- .as_cipcw(cipcw)
  cw <- NULL
  info <- NULL
  if (!is.null(cipcw)) {
    if (outcome_type != "survival")
      stop("`cipcw=` is for survival outcomes: a Surv() outcome with `time=`.",
           call. = FALSE)
    if (!is.null(prevalence))
      stop("`prevalence=` cannot be combined with `cipcw=`: with cIPCW the prevalence ",
           "is estimated from the data.", call. = FALSE)
    y <- model_frame[[outcome_name]]
    if (!attr(y, "type") %in% c("right", "mright"))
      stop("`cipcw=` needs right-censored survival times, e.g. Surv(time, status).",
           call. = FALSE)
    if (!is.null(cipcw$weight_event) && length(cipcw$weight_event) == n_data) {
      cipcw$weight_event <- cipcw$weight_event[kept]
      cipcw$weight_horizon <- cipcw$weight_horizon[kept]
    }
    cw <- .cipcw_weights(cipcw, data[kept, , drop = FALSE], time = y[, "time"],
                         status = y[, "status"], horizon = time, w = w)
    info <- list(min_prob = cipcw$min_prob)
  }
  list(weights = w, cipcw_weights = cw, cipcw_info = info)
}

# The rates of .calculate_test_consequences() with case weights and/or cIPCW:
# pos_rate (the prevalence), test_pos_rate (the share treated), tp_rate
# (prevalence x sensitivity) and fp_rate ((1 - prevalence) x (1 - specificity)).
.weighted_test_consequences <- function(outcome, risk, thresholds, outcome_type,
                                        prevalence, time, weights, cipcw_weights) {
  # cIPCW: the event of interest is status 1 (for competing risks, the first
  # event of the Surv() outcome, as for the Aalen-Johansen estimate)
  if (!is.null(cipcw_weights)) {
    return(.cipcw_rates(time = outcome[, "time"], status = outcome[, "status"],
                        risk = risk, thresholds = thresholds, horizon = time,
                        cw = cipcw_weights, w = weights, k = 1))
  }
  w <- weights
  treated <- function(z) risk >= z
  df <- tibble::tibble(threshold = thresholds, n = length(outcome))
  df$test_pos_rate <- vapply(thresholds, function(z) sum(w[treated(z)]) / sum(w), numeric(1))

  if (outcome_type == "binary") {
    y <- outcome == "TRUE"
    df$pos_rate <- if (!is.null(prevalence)) prevalence else sum(w[y]) / sum(w)
    sens <- vapply(thresholds, function(z) sum(w[y & treated(z)]) / sum(w[y]), numeric(1))
    fpr <- vapply(thresholds, function(z) sum(w[!y & treated(z)]) / sum(w[!y]), numeric(1))
    df$tp_rate <- sens * df$pos_rate
    df$fp_rate <- fpr * (1 - df$pos_rate)
  } else {
    # survival: weighted Kaplan-Meier (Aalen-Johansen for competing risks)
    if (!is.null(prevalence)) {
      df$pos_rate <- prevalence
    } else {
      prev <- .surv_to_risk(outcome ~ 1, time = time, quiet = TRUE, weights = w)
      if (is.na(prev))
        stop("Cannot calculate outcome prevalence at specified time, ",
             "likely due to no observed data beyond selected time.", call. = FALSE)
      df$pos_rate <- prev
    }
    risk_pos <- vapply(thresholds, function(z) {
      sel <- treated(z) & w > 0
      if (!any(sel)) return(0)
      tryCatch(.surv_to_risk(outcome[sel] ~ 1, time = time, weights = w[sel]),
               error = function(e) NA_real_)
    }, numeric(1))
    df$tp_rate <- risk_pos * df$test_pos_rate
    df$fp_rate <- (1 - risk_pos) * df$test_pos_rate
  }
  df[, c("threshold", "n", "pos_rate", "test_pos_rate", "tp_rate", "fp_rate")]
}
