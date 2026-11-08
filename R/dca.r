#' Perform Decision Curve Analysis
#'
#' Diagnostic and prognostic models are typically evaluated with measures of
#' accuracy that do not address clinical consequences.
#' Decision-analytic techniques allow assessment of clinical outcomes but often
#' require collection of additional information may be cumbersome to apply to
#' models that yield a continuous result. Decision curve analysis is a method
#' for evaluating and comparing prediction models that incorporates clinical
#' consequences, requires only the data set on which the models are tested,
#' and can be applied to models that have either continuous or dichotomous
#' results.
#' The dca function performs decision curve analysis for binary outcomes.
#' Review the
#' [DCA Vignette](https://www.danieldsjoberg.com/dcurves/articles/dca.html)
#' for a detailed walk-through of various applications.
#' Also, see [www.decisioncurveanalysis.org](https://www.mskcc.org/departments/epidemiology-biostatistics/biostatistics/decision-curve-analysis) for more information.
#'
#' @author Daniel D Sjoberg
#'
#' @param formula a formula with the outcome on the LHS and a sum of
#' markers/covariates to test on the RHS
#' @param data a data frame containing the variables in `formula=`.
#' @param thresholds vector of threshold probabilities between 0 and 1.
#' Default is `seq(0, 0.99, by = 0.01)`. Thresholds at zero are replaced
#' with 10e-10.
#' @param label named list of variable labels, e.g. `list(age = "Age, years")`
#' @param harm named list of harms associated with a test. Default is `NULL`
#' @param as_probability character vector including names of variables
#' that will be converted to a probability. Details below.
#' @param time if outcome is survival, `time=` specifies the time the
#' assessment is made
#' @param weights Optional case weights, one per patient: the name of a column
#' of `data=`, or a numeric vector with one weight per row of `data=`
#' (non-negative). Every proportion is then weighted (for a survival outcome,
#' the Kaplan-Meier or Aalen-Johansen estimates are weighted). Use them, e.g.,
#' for samples that do not represent the target population in equal
#' proportions (each patient weighted by the inverse of the probability of
#' being sampled). Default `NULL`: every patient counts once.
#' @param cipcw For a survival outcome, optional conditional inverse probability
#' of censoring weighting (cIPCW), which allows censoring to depend on
#' covariates, instead of Kaplan-Meier: the covariates of the censoring, e.g.
#' `cipcw = ~ age + famhistory` (a Cox model of the censoring times on them is
#' fitted, and the censoring probabilities are truncated at 0.02), or
#' [cipcw_control()] for the other options (another truncation, your own
#' censoring model with a function predicting from it, or given weights). Not
#' combined with `prevalence=`. Details in [cipcw_control()].
#' @param prevalence When `NULL`, the prevalence is estimated from `data=`.
#' If the data passed is a case-control set, the population prevalence
#' may be set with this argument.
#'
#' @section weights and cipcw arguments:
#' With case weights \eqn{v_i} (`weights=`), the prevalence and the rates are
#' weighted proportions; for a binary outcome \eqn{Y},
#' \deqn{\mathrm{tp\_rate}(z) = \frac{\sum_i v_i Y_i 1(\hat\pi_i \ge z)}{\sum_i v_i},
#' \quad \mathrm{fp\_rate}(z) = \frac{\sum_i v_i (1 - Y_i) 1(\hat\pi_i \ge z)}{\sum_i v_i},}
#' where \eqn{\hat\pi_i} is the predicted risk and \eqn{1(\hat\pi_i \ge z)} is 1
#' for the patients treated at threshold \eqn{z}. For a survival outcome, the
#' Kaplan-Meier (or Aalen-Johansen) risk by `time=`, overall and among the
#' patients treated, is the weighted estimate, and the share treated is
#' weighted. Without `weights=`, every \eqn{v_i = 1} and the results are
#' unchanged.
#'
#' With `cipcw=`, a survival outcome is estimated with conditional inverse
#' probability of censoring weights instead of Kaplan-Meier: each patient with
#' an event (of any cause) by `time=` \eqn{\tau} is weighted by
#' \eqn{1/\hat G(\tilde T_i- \mid X_i)}, and each patient event-free at
#' \eqn{\tau} by \eqn{1/\hat G(\tau \mid X_i)}, with \eqn{\tilde T_i} the observed
#' time, \eqn{X_i} the covariates of the censoring model and \eqn{\hat G(t \mid X)}
#' the estimated probability of remaining uncensored up to \eqn{t}; patients
#' censored before \eqn{\tau} get no weight. The prevalence is then the weighted
#' share of patients with the event of interest by \eqn{\tau}, and
#' `tp_rate` and `fp_rate` are the prevalence times the weighted sensitivity,
#' and one minus the prevalence times one minus the weighted specificity (a
#' competing event by \eqn{\tau} counts as no event of interest). Case weights
#' multiply the censoring weights. See [cipcw_control()].
#'
#' @section as_probability argument:
#' While the `as_probability=` argument can be used to convert a marker to the
#' probability scale, use the argument only when the consequences are fully
#' understood. For example, when the outcome is binary, logistic regression
#' is used to convert the marker to a probability. The logistic regression
#' model assumes linearity on the log-odds scale and can induce
#' miscalibration when this assumption is not true. Miscalibration in a
#' model will adversely affect performance on decision curve
#' analysis. Similarly, when the outcome is time-to-event, Cox Proportional
#' Hazards regression is used to convert the marker to a probability.
#' The Cox model also has a linearity assumption and additionally assumes
#' proportional hazards over the follow-up period. When these assumptions
#' are violated, important miscalibration may occur.
#'
#' Instead of using the `as_probability=` argument, it is suggested to perform
#' the regression modeling outside of the `dca()` function utilizing methods,
#' such as non-linear modeling, as appropriate.
#'
#' @return List including net benefit of each variable
#' @seealso [`net_intervention_avoided()`], [`standardized_net_benefit()`], [`plot.dca()`],
#' [`as_tibble.dca()`]
#' @export
#'
#' @examples
#' # calculate DCA with binary endpoint
#' dca(cancer ~ cancerpredmarker + marker,
#'     data = df_binary,
#'     as_probability = "marker",
#'     label = list(cancerpredmarker = "Prediction Model", marker = "Biomarker")) %>%
#'   # plot DCA curves with ggplot
#'   plot(smooth = TRUE) +
#'   # add ggplot formatting
#'   ggplot2::labs(x = "Treatment Threshold Probability")
#'
#' # calculate DCA with time to event endpoint
#' dca(Surv(ttcancer, cancer) ~ cancerpredmarker, data = df_surv, time = 1)
dca <- function(formula, data, thresholds = seq(0, 0.99, by = 0.01),
                label = NULL, harm = NULL, as_probability = character(),
                time = NULL, prevalence = NULL, weights = NULL, cipcw = NULL) {
  # checking inputs ------------------------------------------------------------
  if (!is.data.frame(data))
    stop("`data=` must be a data frame", call. = FALSE)
  if (!inherits(formula, "formula"))
    stop("`formula=` must be a formula", call. = FALSE)
  if (!is.null(label) && (!rlang::is_list(label) || !rlang::is_named(label)))
    stop("`label=` must be a named list or NULL.", call. = FALSE)
  if (!is.null(harm) && (!rlang::is_list(harm) || !rlang::is_named(harm)))
    stop("`harm=` must be a named list or NULL.", call. = FALSE)
  if (!is.character(as_probability))
    stop("`as_probability=` must be character.", call. = FALSE)

  # prepping data --------------------------------------------------------------
  thresholds <- thresholds[thresholds >= 0 & thresholds <= 1]

  label <-
    list(all = "Treat All", none = "Treat None") %>%
    purrr::list_modify(!!!label)
  model_frame <- stats::model.frame(formula, data)
  outcome_name <- names(model_frame)[1]
  if (any(c("all", "none") %in% names(model_frame))) {
    stop("Variables cannot be named 'all' or 'none': they are reserved.",
         call. = FALSE)
  }

  outcome_type <- .outcome_type(model_frame, outcome_name, time)

  # case weights and cIPCW ------------------------------------------------------
  wc <- .weights_and_cipcw(model_frame, data, outcome_name, outcome_type, weights,
                           cipcw, time, prevalence)

  # convert to probability if requested ----------------------------------------
  model_frame <- .as_probability(model_frame, outcome_name, outcome_type, as_probability,
                                 time, weights = wc$weights)

  # add treat all and treat none -----------------------------------------------
  model_frame <-
    model_frame %>%
    dplyr::mutate(
      all = 1,
      none = 0,
      dplyr::across(.cols = -dplyr::all_of(outcome_name),
                    .fns = ~dplyr::case_when(. == 0 ~ 0 - .Machine$double.eps,
                                             . == 1 ~ 1 + .Machine$double.eps,
                                             TRUE ~ .)),
      .after = all_of(outcome_name)
    ) %>%
    copy_labels_from(model_frame)

  # calculate net benefit ------------------------------------------------------
  dca_result <-
    test_consequences_data_frame(
      model_frame = model_frame,
      outcome_name = outcome_name,
      outcome_type = outcome_type,
      statistics = c("pos_rate", "tp_rate", "fp_rate", "harm", "net_benefit"),
      thresholds = thresholds,
      label = label,
      time = time,
      prevalence = prevalence,
      harm = harm,
      weights = wc$weights,
      cipcw_weights = wc$cipcw_weights
    )

  # return results -------------------------------------------------------------
  lst_result <-
    list(
      call = match.call(),
      y = outcome_name,
      n = dca_result$n[1],
      prevalence = dca_result$pos_rate[1],
      time = time,
      weights = if (!is.null(wc$weights)) TRUE,
      cipcw = wc$cipcw_info,
      dca = dca_result
    ) %>%
    purrr::compact()
  class(lst_result) <- "dca"
  lst_result
}


#' Calculate risks
#'
#' When users request, we can convert a marker to a probability with logistic, or Cox regression
#'
#' @param outcome outcome object
#' @param variable variable
#' @param outcome_type type of outcome
#' @param time time to calculate risk if Surv() outcome
#' @param prevalence specified prevalence (if cannot be estimated from data)
#'
#' @noRd
#' @keywords internal
.convert_to_risk <- function(outcome, variable, outcome_type,
                             time = NULL, prevalence = NULL, weights = NULL) {
  if (outcome_type == "binary" && !is.null(prevalence)) {
    stop("Cannot convert to risks in case-control setting.")
  }

  if (outcome_type == "binary") {
    # with case weights, quasibinomial gives the weighted logistic fit without
    # warnings for non-integer weights
    risk <-
      stats::glm(outcome ~ variable,
                 family = if (is.null(weights)) stats::binomial else stats::quasibinomial,
                 weights = weights) %>%
      stats::predict(type = "response")
  } else if (outcome_type == "survival") {
    # construct data frame
    df <- data.frame(outcome = outcome, variable = variable)
    new_df <- data.frame(outcome = outcome, variable = variable)
    new_df$outcome[, 1] <- time
    # build model, and get predictions for time point of interest
    risk <-
      survival::coxph(outcome ~ variable, data = df, weights = weights) %>%
      stats::predict(newdata = new_df, type = "expected") %>%
      {
        exp(-.)
      }
  }

  risk
}


#' Convert a `Surv()` object to a n-time risk estimate
#'
#' @param outcome `Surv()` object
#' @param time numeric time
#' @param quiet logical, default is TRUE
#'
#' @noRd
#' @keywords internal
.surv_to_risk <- function(outcome, time, quiet = TRUE, weights = NULL) {
  fit <- if (is.null(weights)) {
    survival::survfit(outcome)
  } else {
    # weighted Kaplan-Meier (Aalen-Johansen for competing risks)
    y <- eval(outcome[[2]], environment(outcome))
    survival::survfit(y ~ 1, weights = weights)
  }
  df_tidy <- broom::tidy(fit)

  # if multistate (i.e. competing risks) delete states not of interest
  if ("state" %in% names(df_tidy)) {
    state <- unique(df_tidy$state) %>%
      setdiff("(s0)") %>%
      purrr::pluck(1)
    df_tidy <- df_tidy %>% dplyr::filter(.data$state %in% .env$state)
    if (!isTRUE(quiet)) {
      glue::glue(
        "Multi-state model detected. Showing probabilities into state '{state}'") %>%
        message()
    }
  }
  # if regular survfit() model, convert survival to risk
  else {
    df_tidy <- dplyr::mutate(df_tidy, estimate = 1 - .data$estimate)
  }

  # if no observed times after specified time, return NA
  if (max(df_tidy$time) < time) {
    return(NA_real_)
  }

  df_tidy <-
    df_tidy %>%
    dplyr::filter(.data$time <= .env$time)

  # if no observed times after time point, return NA
  if (nrow(df_tidy) == 0L) {
    return(NA_real_)
  }

  df_tidy %>%
    dplyr::slice_tail() %>%
    dplyr::pull("estimate")
}


# convert marker to a probability
.as_probability <- function(model_frame, outcome_name, outcome_type, as_probability, time,
                            weights = NULL) {
  as_probability <-
    model_frame %>%
    dplyr::select(-dplyr::all_of(outcome_name)) %>%
    dplyr::select(dplyr::all_of(as_probability))
  for (v in names(as_probability)) {
    model_frame[[v]] <- .convert_to_risk(model_frame[[outcome_name]],
                                         model_frame[[v]],
                                         outcome_type = outcome_type,
                                         time = time,
                                         weights = weights
    )
  }
  .check_probability_range(model_frame, outcome_name)

  model_frame
}

# a preceding step may result in the loss of label attributes
# this function helps restore them
copy_labels_from <- function(x, y) {
  for (v in names(x)) {
    attr(x[[v]], "label") <- attr(y[[v]], "label") %||% attr(x[[v]], "label")
  }
  x
}

