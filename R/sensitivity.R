#' Does screening change the evidence? Fit the model with and without flagged records
#'
#' Fits the same model to all records and to the retained records, and
#' compares fit indices with measurement evidence (loadings, omega,
#' omega-hierarchical, explained common variance). In operational data,
#' records with no real responses can raise loadings and reliability while
#' leaving CFI and RMSEA almost unchanged, so fit indices alone do not reveal
#' the problem.
#'
#' @param items Item scores (data frame or matrix), one row per record.
#' @param keep Logical vector: records to retain (for example
#'   `screen_records(...)$flags$keep`).
#' @param parcels Named list of item columns per parcel (indicator). If
#'   `NULL`, the item columns are used as indicators.
#' @param model lavaan syntax for a second-order model whose indicators are the
#'   parcel names.
#' @param estimator Estimator for [lavaan::cfa()]. Default `"MLR"`.
#' @param std.lv Fix latent variances to 1? Default `TRUE`.
#' @return A data frame of class `fc_sensitivity` with one column per data
#'   version (all records, retained records) and one row per index. The row
#'   `admissible` is 1 for an admissible solution and 0 otherwise. Omega,
#'   omega-hierarchical and ECV are computed only from admissible solutions;
#'   for an improper solution they are `NA`.
#' @references
#' Rashid, S., Retnawati, H., & Hadi, S. (2026). When all-zero records stay in:
#' Data screening, measurement evidence, and fit-index behavior in an
#' operational English proficiency test. Unpublished manuscript.
#' @examples
#' data(sim_test)
#' X <- sim_test[, sim_test_design$items]
#' s <- screen_records(X, sim_test_design$sections)
#' screening_sensitivity(X, s$flags$keep, sim_test_design$parcels, sim_test_design$model)
#' @export
screening_sensitivity <- function(items, keep, parcels = NULL, model, estimator = "MLR", std.lv = TRUE) {
  X <- as_item_matrix(items)
  if (length(keep) != nrow(X)) stop("`keep` must have one value per record.", call. = FALSE)
  D <- if (is.null(parcels)) as.data.frame(X) else parcel_scores(X, parcels)
  one <- function(d) {
    fit <- tryCatch(suppressWarnings(lavaan::cfa(model, data = d, estimator = estimator, std.lv = std.lv)),
                    error = function(e) NULL)
    if (is.null(fit) || !isTRUE(lavaan::lavInspect(fit, "converged"))) return(NULL)
    ck <- four_checks(fit, verbose = FALSE)
    hix <- if (isTRUE(ck$admissible))
      tryCatch(hierarchical_indices(fit), error = function(e) NULL, warning = function(w) NULL) else NULL
    std <- lavaan::standardizedSolution(fit)
    pe <- lavaan::parameterEstimates(fit); lv <- lavaan::lavNames(fit, "lv")
    hi <- higher_order_factors(pe, lv)
    so <- std[std$op == "=~" & std$lhs %in% hi, ]
    v <- c(N = nrow(d), admissible = as.numeric(ck$admissible),
           CFI = ck$cfi, TLI = ck$tli, RMSEA = ck$rmsea, SRMR = ck$srmr, uSRMR = ck$usrmr,
           mean_first_order_loading = ck$mean_first_order,
           stats::setNames(so$est.std, paste0("second_order_", so$rhs)))
    if (!is.null(hix)) v <- c(v, stats::setNames(hix$omega, paste0("omega_", names(hix$omega))),
                             omega_h = hix$omega_h, ECV = hix$ecv)
    v
  }
  a <- one(D); b <- one(D[keep, , drop = FALSE])
  if (is.null(a) || is.null(b)) stop("The model did not converge in one of the data versions.", call. = FALSE)
  rn <- union(names(a), names(b))
  out <- data.frame(index = rn, all_records = unname(a[rn]), retained = unname(b[rn]))
  out$difference <- out$retained - out$all_records
  class(out) <- c("fc_sensitivity", "data.frame")
  out
}

#' @export
print.fc_sensitivity <- function(x, digits = 3, ...) {
  y <- as.data.frame(unclass(x)); num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], function(v) round(v, digits))
  cat("Screening sensitivity: all records vs retained records\n\n")
  print(y, row.names = FALSE)
  adm <- y[y$index == "admissible", c("all_records", "retained")]
  if (nrow(adm) && any(unlist(adm) == 0)) {
    cat("\nWarning: the solution is improper in",
        paste(c("all records", "retained records")[unlist(adm) == 0], collapse = " and "),
        "- omega and ECV are not computed there, and its fit indices should not be read.\n")
  }
  fit_rows <- y$index %in% c("CFI", "RMSEA")
  meas_rows <- grepl("^omega|^ECV|loading", y$index)
  if (any(fit_rows) && any(meas_rows)) {
    cat("\nLargest change in CFI/RMSEA:", max(abs(y$difference[fit_rows]), na.rm = TRUE),
        "| largest change in loadings/omega/ECV:", max(abs(y$difference[meas_rows]), na.rm = TRUE), "\n")
  }
  invisible(x)
}
