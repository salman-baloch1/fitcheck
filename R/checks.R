#' Checks 1 to 4 for one fitted second-order model
#'
#' Reads a fitted lavaan model in a fixed order: estimation (is the solution
#' admissible?), context (sample size and loading magnitude), estimates (fit
#' indices, with the unbiased SRMR and its interval), and decision aids (the
#' largest residual correlation). Fit indices are reported for an inadmissible
#' solution but flagged, because fit computed from an improper solution should
#' not be interpreted.
#'
#' The unbiased SRMR (uSRMR), its 90% interval, and the test of exact fit come
#' from [lavaan::lavResiduals()] (Maydeu-Olivares, 2017). A uSRMR of exactly 0
#' means the corrected value was truncated at the boundary; it is not evidence
#' of exact fit.
#'
#' @param fit A fitted [lavaan::lavaan-class] object (typically from
#'   [lavaan::cfa()] with `estimator = "MLR"`).
#' @param verbose Print the checks? Default `TRUE`.
#' @return An object of class `fc_checks` (a list), returned invisibly when
#'   printed.
#' @references
#' Maydeu-Olivares, A. (2017). Assessing the size of model misfit in
#' structural equation models. *Psychometrika, 82*(3), 533-558.
#'
#' Rashid, S. (2026). *The impact of sample size and factor loadings on
#' goodness-of-fit indices in second-order confirmatory factor analysis: A
#' bootstrapping study of the ProTEFL test data* (Doctoral dissertation).
#' Universitas Negeri Yogyakarta.
#' @seealso [identification_check()], [hierarchical_indices()]
#' @examples
#' S <- make_population(c(3, 2, 5), lambda = .70, gamma = c(.95, .90, .85))
#' set.seed(1)
#' d <- as.data.frame(MASS::mvrnorm(500, rep(0, 10), S))
#' fit <- lavaan::cfa(second_order_syntax(blocks_from_sizes(c(3, 2, 5))),
#'                    data = d, estimator = "MLR")
#' four_checks(fit)
#' @export
four_checks <- function(fit, verbose = TRUE) {
  if (!inherits(fit, "lavaan")) stop("`fit` must be a fitted lavaan object.", call. = FALSE)
  out <- list()
  out$converged <- isTRUE(lavaan::lavInspect(fit, "converged"))
  out$admissible <- out$converged && isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
  pe <- lavaan::parameterEstimates(fit)
  lv <- lavaan::lavNames(fit, "lv")
  var_rows <- pe$op == "~~" & pe$lhs == pe$rhs
  out$neg_disturbance <- unique(pe$lhs[var_rows & pe$lhs %in% lv & pe$est < 0])
  out$neg_residual <- unique(pe$lhs[var_rows & !(pe$lhs %in% lv) & pe$est < 0])
  std <- lavaan::standardizedSolution(fit)
  hi <- higher_order_factors(pe, lv)
  load <- std$op == "=~"
  out$N <- lavaan::lavInspect(fit, "nobs")
  out$mean_first_order <- mean(abs(std$est.std[load & !(std$lhs %in% hi)]))
  out$max_second_order <- if (length(hi)) max(abs(std$est.std[load & std$lhs %in% hi])) else NA
  out$min_first_order <- min(abs(std$est.std[load & !(std$lhs %in% hi)]))
  out$out_of_range <- any(abs(std$est.std[load]) > 1)
  if (out$converged) {
    fm <- lavaan::fitMeasures(fit, c("cfi.robust", "tli.robust", "rmsea.robust", "srmr"))
    rs <- usrmr_summary(fit)
    out$cfi <- unname(fm["cfi.robust"]); out$tli <- unname(fm["tli.robust"])
    out$rmsea <- unname(fm["rmsea.robust"]); out$srmr <- unname(fm["srmr"])
    out$usrmr <- rs$usrmr; out$usrmr_ci <- rs$ci; out$exact_p <- rs$exact_p
    rc <- tryCatch(lavaan::lavResiduals(fit, type = "cor.bentler")$cov, error = function(e) NULL)
    if (!is.null(rc)) {
      rc[upper.tri(rc, diag = TRUE)] <- NA
      w <- which(abs(rc) == max(abs(rc), na.rm = TRUE), arr.ind = TRUE)[1, ]
      out$max_resid <- rc[w[1], w[2]]
      out$max_resid_pair <- paste(rownames(rc)[w[1]], colnames(rc)[w[2]], sep = "-")
    }
  }
  class(out) <- "fc_checks"
  if (verbose) print(out)
  invisible(out)
}

#' @export
print.fc_checks <- function(x, ...) {
  f3 <- function(v) if (is.null(v) || is.na(v)) "NA" else formatC(v, digits = 3, format = "f")
  cat("\nCHECK 1  Estimation\n")
  cat("  Converged:", x$converged, "  Admissible:", x$admissible, "\n")
  if (length(x$neg_disturbance)) cat("  Negative disturbance variance:", x$neg_disturbance, "\n")
  if (length(x$neg_residual)) cat("  Negative residual variance:", x$neg_residual, "\n")
  if (isTRUE(x$out_of_range)) cat("  Standardized loading outside [-1, 1]\n")
  if (!x$admissible) cat("  -> Stop. Diagnose the solution before reading fit (see identification_check()).\n")
  cat("CHECK 2  Context\n")
  cat("  N =", x$N, "  mean first-order loading =", f3(x$mean_first_order),
      "  largest second-order loading =", f3(x$max_second_order), "\n")
  if (!is.na(x$max_second_order) && x$max_second_order > .95)
    cat("  Second-order loading near 1: improper solutions are likely in small samples.\n")
  if (!is.null(x$srmr)) {
    cat("CHECK 3  Estimates", if (!x$admissible) "(from an improper solution: do not interpret)", "\n")
    cat("  SRMR =", f3(x$srmr), "  uSRMR =", f3(x$usrmr),
        " 90% CI [", f3(x$usrmr_ci[1]), ",", f3(x$usrmr_ci[2]), "]\n")
    cat("  Exact-fit test (SRMR) p =", f3(x$exact_p), "\n")
    cat("  robust CFI =", f3(x$cfi), "  robust TLI =", f3(x$tli), "  robust RMSEA =", f3(x$rmsea), "\n")
    if (isTRUE(x$usrmr == 0)) cat("  uSRMR is zero: truncated at the boundary; not evidence of exact fit.\n")
  }
  if (!is.null(x$max_resid)) {
    cat("CHECK 4  Decision aids\n")
    cat("  Largest residual correlation:", f3(x$max_resid), "(", x$max_resid_pair, ")\n")
  }
  cat("\n")
  invisible(x)
}

#' Admissibility under two identifications
#'
#' Fits the same model with latent variances fixed to 1 (`std.lv = TRUE`) and
#' with marker variables (`std.lv = FALSE`), and reports convergence,
#' admissibility, and negative variances for each. The two are the same
#' model, but they can fail differently: with latent variances fixed to 1, a
#' negative first-order disturbance variance often appears as nonconvergence
#' or a boundary solution, whereas marker identification shows it as a
#' Heywood case. Report the identification you used.
#'
#' @param model lavaan model syntax.
#' @param data A data frame.
#' @param estimator Estimator passed to [lavaan::cfa()]. Default `"MLR"`.
#' @param ... Further arguments for [lavaan::cfa()].
#' @return A data frame with one row per identification.
#' @references
#' Rashid, S. (2026). *The impact of sample size and factor loadings on
#' goodness-of-fit indices in second-order confirmatory factor analysis: A
#' bootstrapping study of the ProTEFL test data* (Doctoral dissertation).
#' Universitas Negeri Yogyakarta.
#' @examples
#' S <- make_population(c(3, 2, 5), lambda = .50, gamma = c(.98, .95, .95))
#' set.seed(3)
#' d <- as.data.frame(MASS::mvrnorm(150, rep(0, 10), S))
#' identification_check(second_order_syntax(blocks_from_sizes(c(3, 2, 5))), d)
#' @export
identification_check <- function(model, data, estimator = "MLR", ...) {
  one <- function(std_lv) {
    fit <- tryCatch(suppressWarnings(lavaan::cfa(model, data = data, estimator = estimator,
                                                 std.lv = std_lv, ...)), error = function(e) NULL)
    if (is.null(fit)) return(data.frame(converged = FALSE, admissible = FALSE,
                                        negative_disturbance = NA, negative_residual = NA))
    conv <- isTRUE(lavaan::lavInspect(fit, "converged"))
    ok <- conv && isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
    pe <- lavaan::parameterEstimates(fit); lv <- lavaan::lavNames(fit, "lv")
    v <- pe[pe$op == "~~" & pe$lhs == pe$rhs, ]
    data.frame(converged = conv, admissible = ok,
               negative_disturbance = none_if_empty(v$lhs[v$lhs %in% lv & v$est < 0]),
               negative_residual = none_if_empty(v$lhs[!(v$lhs %in% lv) & v$est < 0]))
  }
  out <- rbind(one(TRUE), one(FALSE))
  rownames(out) <- c("latent variances fixed (std.lv = TRUE)", "marker variables (std.lv = FALSE)")
  out
}

#' Omega, omega-hierarchical, and explained common variance for a second-order model
#'
#' Uses the Schmid-Leiman transformation of the standardized second-order
#' solution: the general loading of indicator i on first-order factor f is
#' lambda_if * gamma_f, and its group loading is lambda_if * sqrt(1 - gamma_f^2).
#' Assumes simple structure (each indicator loads on one first-order factor)
#' and an admissible solution.
#'
#' @param fit A fitted lavaan second-order model.
#' @return A list with `omega` (per first-order factor and total),
#'   `omega_h` (general factor), `omega_hs` (per first-order factor), and
#'   `ecv` (explained common variance of the general factor).
#' @references
#' Rodriguez, A., Reise, S. P., & Haviland, M. G. (2016). Evaluating bifactor
#' models: Calculating and interpreting statistical indices.
#' *Psychological Methods, 21*(2), 137-150.
#' @examples
#' S <- make_population(c(3, 2, 5), lambda = .70, gamma = c(.95, .90, .85))
#' set.seed(1)
#' d <- as.data.frame(MASS::mvrnorm(800, rep(0, 10), S))
#' fit <- lavaan::cfa(second_order_syntax(blocks_from_sizes(c(3, 2, 5))),
#'                    data = d, estimator = "MLR")
#' hierarchical_indices(fit)
#' @export
hierarchical_indices <- function(fit) {
  pe <- lavaan::parameterEstimates(fit); lv <- lavaan::lavNames(fit, "lv")
  hi <- higher_order_factors(pe, lv)
  if (length(hi) != 1) stop("Expected exactly one second-order factor.", call. = FALSE)
  std <- lavaan::standardizedSolution(fit)
  g <- std[std$op == "=~" & std$lhs == hi, c("rhs", "est.std")]
  fo <- std[std$op == "=~" & std$lhs %in% g$rhs, c("lhs", "rhs", "est.std")]
  gam <- g$est.std[match(fo$lhs, g$rhs)]
  if (any(abs(gam) >= 1)) warning("A second-order loading is 1 or more: indices are not meaningful.", call. = FALSE)
  gl <- fo$est.std * gam
  sl <- fo$est.std * sqrt(pmax(0, 1 - gam^2))
  R <- stats::cov2cor(lavaan::lavInspect(fit, "implied")$cov)
  R <- R[fo$rhs, fo$rhs]
  tot <- sum(R)
  by_f <- split(seq_len(nrow(fo)), fo$lhs)
  om <- vapply(by_f, function(i) (sum(gl[i])^2 + sum(sl[i])^2) / sum(R[i, i]), numeric(1))
  omhs <- vapply(by_f, function(i) sum(sl[i])^2 / sum(R[i, i]), numeric(1))
  list(omega = c(om, total = (sum(gl)^2 + sum(vapply(by_f, function(i) sum(sl[i])^2, numeric(1)))) / tot),
       omega_h = sum(gl)^2 / tot, omega_hs = omhs,
       ecv = sum(gl^2) / (sum(gl^2) + sum(sl^2)))
}

# ---- internal helpers ----------------------------------------------------
higher_order_factors <- function(pe, lv) unique(pe$lhs[pe$op == "=~" & pe$rhs %in% lv])

usrmr_summary <- function(fit) {
  rs <- tryCatch(lavaan::lavResiduals(fit, type = "cor.bentler", summary = TRUE)$summary,
                 error = function(e) NULL)
  if (is.null(rs)) return(list(usrmr = NA, ci = c(NA, NA), exact_p = NA))
  g <- function(n) if (n %in% rownames(rs)) rs[n, 1] else NA
  list(usrmr = g("usrmr"), ci = c(max(0, g("usrmr.ci.lower")), g("usrmr.ci.upper")),
       exact_p = g("srmr.exactfit.pvalue"))
}

# "none" when no variable is listed, so blank cells are not misread.
none_if_empty <- function(x) if (length(x)) paste(x, collapse = " ") else "none"
