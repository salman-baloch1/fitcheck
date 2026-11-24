#' Simulated item scores from a three-section operational test
#'
#' Simulated data (no real examinee records) modelled on a three-section
#' university English test: 2,000 records, three forms, and 80 four-option
#' items scored 1/0 in Listening (30), Grammar (20), and Reading (30).
#' Responses follow a three-parameter logistic model with a second-order
#' factor structure. As in many operational files, some records contain no
#' real responses: all-zero records (concentrated in form B) and records with
#' every Listening item scored 0. Use [screen_records()] to find them.
#'
#' @format A data frame with 2,000 rows and 82 columns: `id`, `form` (A, B,
#'   C), and items `i1` to `i80`.
#' @source Simulated; see `data-raw/make_sim_test.R` in the package sources.
#' @seealso [sim_test_design]
"sim_test"

#' Design of the simulated test
#'
#' A list describing [sim_test]: `items` (item column names), `sections`
#' (items per section), `parcels` (items per parcel: three Listening, two
#' Grammar, five Reading), and `model` (lavaan syntax for the second-order
#' model of the ten parcels).
#'
#' @format A list of four elements.
#' @source Simulated; see `data-raw/make_sim_test.R`.
"sim_test_design"

#' fitcheck: Check the Data, Then the Solution, Then the Fit
#'
#' Tools for judging second-order confirmatory factor models fitted to
#' operational test data in a fixed order: screen the records (Check 0),
#' confirm that the solution is admissible (Check 1), read the fit indices in
#' context with the unbiased SRMR and its interval (Checks 2 to 4), see how the
#' model would behave in smaller samples, and decide whether section scores
#' deserve to be reported. Built on \pkg{lavaan}.
#'
#' Start with [fitcheck()] for the whole workflow, or use the steps
#' separately: [screen_records()], [screening_sensitivity()],
#' [four_checks()], [identification_check()], [hierarchical_indices()],
#' [boot_by_n()], [mini_mc()], [added_value()], and [profile_check()].
#'
#' @keywords internal
"_PACKAGE"
