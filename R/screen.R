#' Check 0: flag records before any model is fitted
#'
#' Flags records in which every item is scored 0 (probable non-attempts) and
#' records in which every item of at least one section is scored 0 while the
#' rest of the record is not (possible partial non-attempts). The function
#' also returns the low-score distribution of each section, which shows
#' whether zeros are more common than scores of 1 or 2 (the "spike at zero").
#'
#' Flagged records are potentially anomalous, not proven invalid. When their
#' status cannot be verified from attendance or administration records, report
#' the main results with and without them (see [screening_sensitivity()]).
#'
#' @param items A data frame or matrix of dichotomous item scores (1/0), one
#'   row per record. Missing values are not allowed.
#' @param sections A named list giving the item columns (numbers or names) of
#'   each section, for example `list(Listening = 1:50, Grammar = 51:90)`.
#' @param form Optional vector with the test form of each record. Counts are
#'   then reported by form.
#' @param n_options Number of response options per item, used to judge how
#'   unusual a zero section score is under random guessing. Default `4`.
#' @param max_p_zero A whole-section zero is flagged only for sections in which
#'   the probability of a zero score under random guessing
#'   ([zero_probability()]) is below this value. Default `.01`. In short
#'   sections (for example 4 items with 4 options, probability .32) a zero
#'   score is common among genuine low scorers, so it is not flagged. Set to
#'   `1` to flag whole-section zeros in every section.
#' @return An object of class `fc_screen`: a list with
#'   \describe{
#'     \item{flags}{a data frame with one row per record: `form`, `all_zero`,
#'       `section_zero`, `keep` (neither flag), and one logical column per
#'       section (`zero_<section>`).}
#'     \item{by_form}{counts of flagged records by form.}
#'     \item{low_scores}{for each section, the number of non-all-zero records
#'       scoring 0, 1, and 2, the spike ratio (0 against 1 or 2), the
#'       probability of a zero section score under random guessing, and
#'       whether whole-section zeros in that section are flagged.}
#'   }
#' @references
#' Rashid, S., Retnawati, H., & Hadi, S. (2026). When all-zero records stay in:
#' Data screening, measurement evidence, and fit-index behavior in an
#' operational English proficiency test. Unpublished manuscript.
#' @seealso [screening_sensitivity()], [zero_probability()]
#' @examples
#' data(sim_test)
#' s <- screen_records(sim_test[, sim_test_design$items],
#'                     sim_test_design$sections, form = sim_test$form)
#' s
#' table(s$flags$keep)
#' @export
screen_records <- function(items, sections, form = NULL, n_options = 4, max_p_zero = .01) {
  X <- as_item_matrix(items)
  sections <- resolve_sections(sections, colnames(X))
  if (is.null(form)) form <- rep("all", nrow(X))
  if (length(form) != nrow(X)) stop("`form` must have one value per record.", call. = FALSE)
  sec <- vapply(sections, function(i) rowSums(X[, i, drop = FALSE]), numeric(nrow(X)))
  if (is.null(dim(sec))) sec <- matrix(sec, nrow = 1, dimnames = list(NULL, names(sections)))
  all_zero <- rowSums(X) == 0
  zero_sec <- sec == 0 & !all_zero
  eligible <- as.vector(zero_probability(lengths(sections), n_options)) < max_p_zero
  sec_zero <- rowSums(zero_sec[, eligible, drop = FALSE]) > 0
  flags <- data.frame(form = form, all_zero = all_zero, section_zero = sec_zero,
                      keep = !(all_zero | sec_zero), zero_sec, check.names = FALSE)
  names(flags)[-(1:4)] <- paste0("zero_", names(sections))
  by_form <- do.call(rbind, lapply(split(flags, factor(flags$form, levels = sort(unique(flags$form)))), function(d)
    data.frame(form = as.character(d$form[1]), records = nrow(d), all_zero = sum(d$all_zero),
               section_zero = sum(d$section_zero),
               flagged_pct = round(100 * mean(!d$keep), 1))))
  rownames(by_form) <- NULL
  low <- t(vapply(names(sections), function(k) {
    s <- sec[!all_zero, k]
    c(score_0 = sum(s == 0), score_1 = sum(s == 1), score_2 = sum(s == 2))
  }, numeric(3)))
  low <- data.frame(section = rownames(low), low, row.names = NULL)
  low$spike_ratio <- round(low$score_0 / pmax(low$score_1 + low$score_2, 1), 2)
  low$p_zero_guessing <- signif(as.vector(zero_probability(lengths(sections), n_options)), 3)
  low$flagged <- eligible
  structure(list(flags = flags, by_form = by_form, low_scores = low,
                 sections = sections, n = nrow(X)), class = "fc_screen")
}

#' @export
print.fc_screen <- function(x, ...) {
  cat("Check 0: screening of", x$n, "records\n\n")
  print(x$by_form, row.names = FALSE)
  cat("\nLow scores by section (records without an all-zero total):\n")
  print(x$low_scores, row.names = FALSE)
  short <- x$low_scores$section[!x$low_scores$flagged]
  if (length(short))
    cat("\nWhole-section zeros are not flagged in:", paste(short, collapse = ", "),
        "\n(a zero score there is not unusual under random guessing).\n")
  spikes <- x$low_scores$section[x$low_scores$spike_ratio > 1]
  if (length(spikes))
    cat("\nZeros outnumber scores of 1 or 2 in:", paste(spikes, collapse = ", "),
        "\n-> Check these records against administration records before analysis.\n")
  cat("\nFlagged records are potentially anomalous, not proven invalid.",
      "\nReport results with and without them (screening_sensitivity()).\n")
  invisible(x)
}

#' Probability of a zero score under random guessing
#'
#' The probability that a candidate who chooses among `n_options` options at
#' random answers none of `n_items` independent items correctly. It only
#' illustrates how unusual a zero score is; real candidates do not guess
#' independently and at random.
#'
#' The function also returns, in its `detail` attribute, the threshold used in
#' the local-independence argument: if responses are independent given
#' ability, a score of 0 is more likely than a score of 1 only when the
#' average probability of a correct answer is below `1 / n_items`.
#'
#' @param n_items Number of items (a vector is allowed).
#' @param n_options Number of response options per item.
#' @return A numeric vector of probabilities.
#' @examples
#' zero_probability(c(50, 40))   # about 6e-07 and 1e-05
#' @export
zero_probability <- function(n_items, n_options = 4) {
  p <- (1 - 1 / n_options)^n_items
  attr(p, "detail") <- data.frame(n_items = n_items,
                                  max_mean_p_for_zero_more_likely = 1 / n_items)
  p
}

#' Sum item scores into parcels
#'
#' @param items A data frame or matrix of item scores.
#' @param parcels A named list of item columns per parcel.
#' @return A data frame with one column per parcel.
#' @examples
#' data(sim_test)
#' head(parcel_scores(sim_test[, sim_test_design$items], sim_test_design$parcels))
#' @export
parcel_scores <- function(items, parcels) {
  X <- as_item_matrix(items)
  parcels <- resolve_sections(parcels, colnames(X))
  as.data.frame(lapply(parcels, function(i) rowSums(X[, i, drop = FALSE])))
}

# ---- internal helpers ----------------------------------------------------
as_item_matrix <- function(items) {
  X <- as.matrix(items)
  if (!is.numeric(X)) stop("`items` must be numeric.", call. = FALSE)
  if (anyNA(X)) stop("`items` contains missing values; handle them before screening.", call. = FALSE)
  if (!all(X %in% c(0, 1)))
    stop("`items` must contain dichotomous scores (0 or 1) only.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("i", seq_len(ncol(X)))
  X
}

resolve_sections <- function(sections, cn) {
  if (!is.list(sections) || is.null(names(sections)) || any(names(sections) == ""))
    stop("`sections`/`parcels` must be a named list.", call. = FALSE)
  lapply(sections, function(s) {
    idx <- if (is.character(s)) match(s, cn) else as.integer(s)
    if (anyNA(idx) || any(idx < 1 | idx > length(cn))) stop("Unknown item column in `sections`/`parcels`.", call. = FALSE)
    idx
  })
}
