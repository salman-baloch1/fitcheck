#' Does screening change the comparison of test forms? Invariance with and without flagged records
#'
#' Tests measurement invariance across test forms twice, once with all
#' records and once with the retained records, and reports whether screening
#' changes the conclusion. Records without real responses are often
#' concentrated in some forms, so they can create, hide, or change
#' differences between forms. The function shows which of these happens.
#'
#' For each data version, three multiple-group models are fitted with
#' [lavaan::cfa()]: a configural model (same pattern of loadings in every
#' form), a metric model (equal loadings), and a scalar model (equal loadings
#' and intercepts). Each step is compared with the previous one by the change
#' in robust CFI and robust RMSEA. A step is marked "not supported" when CFI
#' decreases by more than `cfi_cut` and RMSEA increases by more than
#' `rmsea_cut` (Chen, 2007).
#'
#' If `model` contains a second-order factor, the higher-order part is
#' removed and the first-order factors are allowed to correlate. With three
#' first-order factors this correlated model is equivalent to the
#' second-order model, and it avoids identification problems of the
#' higher-order factor in multiple-group models. With more or fewer than three
#' first-order factors the two models are not equivalent; the function then
#' warns, and the results refer to the correlated first-order model.
#'
#' Each model is also checked for admissibility (no negative variances and a
#' positive definite latent covariance matrix). A level whose solution is
#' improper is marked "improper solution" and is not compared.
#'
#' The rule is a threshold-based assessment, not proof of invariance, and the
#' decision depends on the cutoffs and the model.
#'
#' @param items Item scores (data frame or matrix), one row per record.
#' @param keep Logical vector: records to retain (for example
#'   `screen_records(...)$flags$keep`).
#' @param form Test form of each record. At least two forms are needed.
#' @param parcels Named list of item columns per parcel (indicator). If
#'   `NULL`, the item columns are used as indicators.
#' @param model lavaan syntax whose indicators are the parcel names. A
#'   second-order model is turned into the correlated first-order model.
#' @param estimator Estimator for [lavaan::cfa()]. Default `"MLR"`.
#' @param cfi_cut Largest acceptable decrease in robust CFI. Default `.010`.
#' @param rmsea_cut Largest acceptable increase in robust RMSEA. Default `.015`.
#' @return An object of class `fc_invariance`: a list with
#'   \describe{
#'     \item{table}{one row per data version and invariance level, with
#'       convergence, admissibility, robust fit indices, the changes from the
#'       previous level, and the decision.}
#'     \item{verdict}{the highest level supported in each data version.}
#'     \item{same_conclusion}{`TRUE` when screening does not change the
#'       highest level supported.}
#'     \item{model}{the first-order model that was fitted.}
#'   }
#' @references
#' Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of
#' measurement invariance. *Structural Equation Modeling, 14*(3), 464-504.
#' \doi{10.1080/10705510701301834}
#'
#' Rashid, S. (2026). *The impact of sample size and factor loadings on
#' goodness-of-fit indices in second-order confirmatory factor analysis: A
#' bootstrapping study of the ProTEFL test data* (Doctoral dissertation).
#' Universitas Negeri Yogyakarta.
#' @seealso [screening_by_form()], [screening_sensitivity()], [screen_records()]
#' @examples
#' \donttest{
#' data(sim_test)
#' X <- sim_test[, sim_test_design$items]
#' s <- screen_records(X, sim_test_design$sections, form = sim_test$form)
#' screening_invariance(X, s$flags$keep, sim_test$form,
#'                      sim_test_design$parcels, sim_test_design$model)
#' }
#' @export
screening_invariance <- function(items, keep, form, parcels = NULL, model, estimator = "MLR",
                                 cfi_cut = .010, rmsea_cut = .015) {
  X <- as_item_matrix(items)
  if (length(keep) != nrow(X)) stop("`keep` must have one value per record.", call. = FALSE)
  if (length(form) != nrow(X)) stop("`form` must have one value per record.", call. = FALSE)
  if (length(unique(form)) < 2) stop("At least two forms are needed.", call. = FALSE)
  D <- if (is.null(parcels)) as.data.frame(X) else parcel_scores(X, parcels)
  D$.form <- as.character(form)
  m1 <- first_order_model(model)
  n_first <- attr(m1, "n_first"); n_higher <- attr(m1, "n_higher")
  if (n_higher > 0 && n_first != 3)
    warning("The model has ", n_first, " first-order factors. The correlated first-order model is ",
            "equivalent to the second-order model only with three first-order factors; the results ",
            "refer to the correlated first-order model.", call. = FALSE)
  levels <- list(configural = character(0), metric = "loadings",
                 scalar = c("loadings", "intercepts"))
  one <- function(d, version) {
    rows <- lapply(names(levels), function(lev) {
      fit <- tryCatch(suppressWarnings(lavaan::cfa(m1, data = d, group = ".form", estimator = estimator,
                                                   group.equal = levels[[lev]])),
                      error = function(e) NULL)
      ok <- !is.null(fit) && isTRUE(lavaan::lavInspect(fit, "converged"))
      adm <- ok && isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
      fm <- if (ok) lavaan::fitMeasures(fit, c("chisq.scaled", "df", "cfi.robust", "rmsea.robust", "srmr"))
            else stats::setNames(rep(NA_real_, 5), c("chisq.scaled", "df", "cfi.robust", "rmsea.robust", "srmr"))
      data.frame(version = version, level = lev, N = nrow(d), converged = ok, admissible = adm,
                 chisq_scaled = unname(fm["chisq.scaled"]), df = unname(fm["df"]),
                 CFI = unname(fm["cfi.robust"]), RMSEA = unname(fm["rmsea.robust"]),
                 SRMR = unname(fm["srmr"]))
    })
    t <- do.call(rbind, rows)
    t$delta_CFI <- c(NA, diff(t$CFI)); t$delta_RMSEA <- c(NA, diff(t$RMSEA))
    t$decision <- ifelse(!t$converged, "not estimated",
                  ifelse(!t$admissible, "improper solution",
                         invariance_decision(t$CFI, t$RMSEA, cfi_cut, rmsea_cut)))
    t
  }
  tab <- rbind(one(D, "all records"), one(D[keep, , drop = FALSE], "retained"))
  verdict <- vapply(c("all records", "retained"), function(v) {
    t <- tab[tab$version == v, ]
    if (!t$converged[1]) return("configural model not estimated")
    if (!t$admissible[1]) return("configural model improper")
    hit <- 1
    for (k in 2:3) if (identical(t$decision[k], "supported")) hit <- k else break
    t$level[hit]
  }, character(1))
  rownames(tab) <- NULL
  structure(list(table = tab, verdict = verdict,
                 same_conclusion = identical(unname(verdict[1]), unname(verdict[2])),
                 model = m1, cuts = c(cfi = cfi_cut, rmsea = rmsea_cut)),
            class = "fc_invariance")
}

#' @export
print.fc_invariance <- function(x, digits = 3, ...) {
  y <- x$table; num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], function(v) round(v, digits))
  cat("Measurement invariance across forms: all records vs retained records\n")
  cat("Model: correlated first-order factors\n")
  cat("Rule (Chen, 2007): the metric or scalar step is not supported when the decrease in robust CFI\n",
      " exceeds", x$cuts[["cfi"]], "AND the increase in robust RMSEA exceeds", x$cuts[["rmsea"]], "\n")
  cat("This is a threshold-based assessment, not proof of invariance; the result depends on the cutoffs.\n\n")
  z <- y[, c("version", "level", "N", "admissible", "CFI", "RMSEA", "delta_CFI", "delta_RMSEA", "decision")]
  names(z)[names(z) %in% c("delta_CFI", "delta_RMSEA")] <- c("dCFI", "dRMSEA")
  print(z, row.names = FALSE)
  cat("\nHighest level supported: all records =", x$verdict[["all records"]],
      "| retained =", x$verdict[["retained"]], "\n")
  if (x$same_conclusion) {
    cat("Screening does not change the invariance conclusion.\n")
  } else {
    cat("Screening changes the invariance conclusion: report both, and check the flagged records",
        "against administration records.\n")
  }
  invisible(x)
}

#' Does screening change the evidence within each form?
#'
#' Repeats [screening_sensitivity()] separately for each test form and
#' collects the main fit and measurement indices in one table. Flagged
#' records are often concentrated in a few forms, so the effect of screening
#' can be large in one form and absent in another.
#'
#' @inheritParams screening_invariance
#' @param form Test form of each record.
#' @param model lavaan syntax for a second-order model whose indicators are
#'   the parcel names (as in [screening_sensitivity()]).
#' @param std.lv Fix latent variances to 1? Default `FALSE` (marker
#'   variables), as in [fitcheck()].
#' @param min_n Forms with fewer retained records than this are skipped.
#'   Default `200`.
#' @return A data frame of class `fc_by_form` with one row per form and index:
#'   `form`, `index`, `all_records`, `retained`, `difference`. Forms whose
#'   model could not be estimated are listed with a note.
#' @seealso [screening_sensitivity()], [screening_invariance()]
#' @examples
#' \donttest{
#' data(sim_test)
#' X <- sim_test[, sim_test_design$items]
#' s <- screen_records(X, sim_test_design$sections, form = sim_test$form)
#' screening_by_form(X, s$flags$keep, sim_test$form,
#'                   sim_test_design$parcels, sim_test_design$model)
#' }
#' @export
screening_by_form <- function(items, keep, form, parcels = NULL, model, estimator = "MLR",
                              std.lv = FALSE, min_n = 200) {
  X <- as_item_matrix(items)
  if (length(keep) != nrow(X)) stop("`keep` must have one value per record.", call. = FALSE)
  if (length(form) != nrow(X)) stop("`form` must have one value per record.", call. = FALSE)
  show <- c("N", "admissible", "CFI", "RMSEA", "SRMR", "mean_first_order_loading", "omega_h", "ECV")
  out <- lapply(sort(unique(as.character(form))), function(f) {
    i <- form == f
    if (sum(keep[i]) < min_n)
      return(data.frame(form = f, index = "note", all_records = NA, retained = NA, difference = NA,
                        note = paste("skipped: fewer than", min_n, "retained records")))
    s <- tryCatch(screening_sensitivity(X[i, , drop = FALSE], keep[i], parcels, model,
                                        estimator = estimator, std.lv = std.lv),
                  error = function(e) conditionMessage(e))
    if (is.character(s))
      return(data.frame(form = f, index = "note", all_records = NA, retained = NA, difference = NA,
                        note = s))
    s <- as.data.frame(unclass(s))
    s <- s[s$index %in% show, ]
    data.frame(form = f, s, note = "")
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  class(out) <- c("fc_by_form", "data.frame")
  out
}

#' @export
print.fc_by_form <- function(x, digits = 3, ...) {
  y <- as.data.frame(unclass(x)); num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], function(v) round(v, digits))
  cat("Screening sensitivity by form: all records vs retained records\n\n")
  wide <- y[y$index != "note", c("form", "index", "all_records", "retained", "difference")]
  print(wide, row.names = FALSE)
  notes <- y[y$index == "note", ]
  if (nrow(notes)) { cat("\n"); for (k in seq_len(nrow(notes))) cat("Form", notes$form[k], ":", notes$note[k], "\n") }
  meas <- wide[wide$index %in% c("mean_first_order_loading", "omega_h", "ECV"), ]
  if (nrow(meas)) {
    big <- tapply(abs(meas$difference), meas$form, max, na.rm = TRUE)
    cat("\nLargest change in loadings/omega_h/ECV by form:",
        paste(names(big), round(big, 3), sep = " = ", collapse = ", "), "\n")
  }
  invisible(x)
}

# Chen (2007) rule: a level is not supported when CFI drops by more than cfi_cut
# and RMSEA rises by more than rmsea_cut, relative to the previous level.
invariance_decision <- function(cfi, rmsea, cfi_cut, rmsea_cut) {
  d_cfi <- c(NA, diff(cfi)); d_rmsea <- c(NA, diff(rmsea))
  out <- ifelse(d_cfi < -cfi_cut & d_rmsea > rmsea_cut, "not supported", "supported")
  out[1] <- "baseline"
  out[is.na(out)] <- "not compared"
  out
}

# Turn a (possibly second-order) model into a correlated first-order model.
first_order_model <- function(model) {
  pt <- lavaan::lavaanify(model)
  lv <- unique(pt$lhs[pt$op == "=~"])
  hi <- unique(pt$lhs[pt$op == "=~" & pt$rhs %in% lv])
  first <- setdiff(lv, hi)
  if (!length(first)) stop("The model has no first-order factors.", call. = FALSE)
  lines <- vapply(first, function(f) {
    ind <- pt$rhs[pt$op == "=~" & pt$lhs == f]
    paste(f, "=~", paste(ind, collapse = " + "))
  }, character(1))
  structure(paste(lines, collapse = "\n"), n_first = length(first), n_higher = length(hi))
}
