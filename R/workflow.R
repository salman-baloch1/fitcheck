#' The whole workflow: check the data, then the solution, then the fit, then the scores
#'
#' Runs the steps of the fitcheck workflow on item-level data from an
#' operational test and returns one report:
#' \enumerate{
#'   \item Check 0, data: [screen_records()] flags all-zero and whole-section
#'     zero records.
#'   \item Screening sensitivity: [screening_sensitivity()] fits the model
#'     with and without the flagged records.
#'   \item If records come from two or more forms: [screening_by_form()]
#'     repeats the comparison within each form, and [screening_invariance()]
#'     tests invariance across forms with and without the flagged records.
#'   \item Checks 1 to 4 on the retained records: [four_checks()] and
#'     [identification_check()].
#'   \item Optional sample-size behaviour: [boot_by_n()].
#'   \item Optional score reporting: [added_value()] and [profile_check()].
#' }
#'
#' @param items Dichotomous item scores (data frame or matrix), one row per record.
#' @param sections Named list of item columns per section.
#' @param parcels Named list of item columns per parcel (the indicators).
#' @param model lavaan syntax. If `NULL`, a second-order model is built with
#'   one first-order factor per section, measured by the parcels whose items
#'   belong to that section.
#' @param form Optional test form of each record.
#' @param subscores Run [added_value()] and [profile_check()]? Default `TRUE`.
#' @param B Bootstrap resamples for [added_value()].
#' @param boot_ns Sample sizes for [boot_by_n()]; `NULL` (default) skips it.
#' @param reps Replications per size for [boot_by_n()].
#' @param cores Cores for [boot_by_n()].
#' @param seed Random seed.
#' @param n_options,max_p_zero Passed to [screen_records()]: whole-section
#'   zeros are flagged only in sections where a zero score is unusual under
#'   random guessing.
#' @param invariance When `form` has two or more forms, run
#'   [screening_by_form()] and [screening_invariance()]? Default `TRUE`.
#' @return An object of class `fc_report`, a list with one element per step.
#' @examples
#' \donttest{
#' data(sim_test)
#' rep <- fitcheck(sim_test[, sim_test_design$items], sim_test_design$sections,
#'                 sim_test_design$parcels, form = sim_test$form, B = 200)
#' rep
#' }
#' @export
fitcheck <- function(items, sections, parcels, model = NULL, form = NULL, subscores = TRUE,
                     B = 1000, boot_ns = NULL, reps = 50, cores = 1, seed = 2025, invariance = TRUE,
                     n_options = 4, max_p_zero = .01) {
  X <- as_item_matrix(items)
  sec_idx <- resolve_sections(sections, colnames(X))
  par_idx <- resolve_sections(parcels, colnames(X))
  if (is.null(model)) model <- model_from_parcels(sec_idx, par_idx)
  out <- list(model = model)
  out$screen <- screen_records(X, sec_idx, form, n_options = n_options, max_p_zero = max_p_zero)
  keep <- out$screen$flags$keep
  out$sensitivity <- tryCatch(screening_sensitivity(X, keep, par_idx, model, std.lv = FALSE),
                              error = function(e) conditionMessage(e))
  if (invariance && !is.null(form) && length(unique(form)) > 1) {
    out$by_form <- tryCatch(screening_by_form(X, keep, form, par_idx, model),
                            error = function(e) conditionMessage(e))
    out$invariance <- tryCatch(screening_invariance(X, keep, form, par_idx, model),
                               error = function(e) conditionMessage(e))
  }
  d <- parcel_scores(X[keep, , drop = FALSE], par_idx)
  fit <- tryCatch(suppressWarnings(lavaan::cfa(model, data = d, estimator = "MLR")), error = function(e) NULL)
  out$fit <- fit
  out$checks <- if (!is.null(fit)) four_checks(fit, verbose = FALSE) else NULL
  out$identification <- identification_check(model, d)
  out$hierarchical <- if (!is.null(out$checks) && out$checks$admissible)
    tryCatch(hierarchical_indices(fit), error = function(e) NULL) else NULL
  if (!is.null(boot_ns)) out$boot <- boot_by_n(d, model, ns = boot_ns, reps = reps, seed = seed, cores = cores)
  if (subscores) {
    out$added_value <- added_value(X[keep, , drop = FALSE], sec_idx, B = B, seed = seed)
    out$profiles <- profile_check(X[keep, , drop = FALSE], sec_idx)
  }
  class(out) <- "fc_report"
  out
}

#' @export
print.fc_report <- function(x, ...) {
  rule <- function(t) cat("\n", strrep("=", 4), " ", t, " ", strrep("=", max(4, 60 - nchar(t))), "\n", sep = "")
  rule("CHECK 0  Data")
  bf <- x$screen$by_form
  cat("Records:", sum(bf$records), "| all-zero:", sum(bf$all_zero), "| whole-section zero:",
      sum(bf$section_zero), "| retained:", sum(x$screen$flags$keep), "\n")
  if (nrow(bf) > 1) print(bf, row.names = FALSE)
  rule("Screening sensitivity")
  if (is.character(x$sensitivity)) cat(x$sensitivity, "\n") else print(x$sensitivity)
  if (!is.null(x$by_form)) {
    rule("Screening sensitivity by form")
    if (is.character(x$by_form)) cat(x$by_form, "\n") else print(x$by_form)
  }
  if (!is.null(x$invariance)) {
    rule("Invariance across forms")
    if (is.character(x$invariance)) cat(x$invariance, "\n") else print(x$invariance)
  }
  rule("CHECKS 1-4  Retained records")
  if (is.null(x$checks)) cat("The model could not be estimated.\n") else print(x$checks)
  cat("Identification check:\n"); print(x$identification)
  if (!is.null(x$hierarchical)) {
    h <- x$hierarchical
    cat("\nomega_h =", round(h$omega_h, 3), "| ECV =", round(h$ecv, 3),
        "| omega_hs:", paste(names(h$omega_hs), round(h$omega_hs, 3), collapse = ", "), "\n")
  }
  if (!is.null(x$boot)) { rule("Sample-size behaviour"); print(x$boot) }
  if (!is.null(x$added_value)) {
    rule("Score reporting"); print(x$added_value)
    cat("\nScore profiles:\n"); print(x$profiles, row.names = FALSE)
  }
  invisible(x)
}

model_from_parcels <- function(sec_idx, par_idx) {
  owner <- vapply(par_idx, function(p) {
    hit <- names(sec_idx)[vapply(sec_idx, function(s) all(p %in% s), logical(1))]
    if (length(hit) != 1) stop("Each parcel must lie within one section.", call. = FALSE)
    hit
  }, character(1))
  blocks <- lapply(names(sec_idx), function(s) names(par_idx)[owner == s])
  names(blocks) <- names(sec_idx)
  if (any(lengths(blocks) == 0)) stop("Every section needs at least one parcel.", call. = FALSE)
  second_order_syntax(blocks)
}
