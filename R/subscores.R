#' Do section scores add value over the total score?
#'
#' For each section, compares the proportional reduction in mean squared error
#' (PRMSE) of the section score with that of the total score as estimators of
#' the true section score (Haberman, 2008), and the gain of the augmented
#' score (the best linear combination of the section score and the rest of the
#' test). A percentile bootstrap interval for PRMSE_s - PRMSE_x gives the
#' decision (in the spirit of Sinharay, 2019): "adds value" if the interval
#' lies above 0, "no added value" if it lies below 0, and "inconclusive"
#' otherwise. All quantities are recomputed in every resample.
#'
#' PRMSE_s equals the reliability of the section score. With
#' `reliability = "kr20"` it is KR-20 (item level). With
#' `reliability = "part"` it is the Feldt-Raju coefficient over the parts
#' given in `parts`, which treats variance specific to a part (for example a
#' reading passage) as error; use it as a sensitivity analysis when items share
#' stimuli. The formulas assume measurement errors that are uncorrelated
#' across sections.
#'
#' @param items Dichotomous item scores (data frame or matrix).
#' @param sections Named list of item columns per section.
#' @param reliability `"kr20"` (default) or `"part"`.
#' @param parts For `reliability = "part"`: a named list (one element per
#'   section) of lists of item columns per part.
#' @param B Number of bootstrap resamples. Use at least 1,000 for reporting.
#' @param seed Random seed.
#' @return A data frame of class `fc_added_value`: one row per section with
#'   PRMSE_s, PRMSE_x, PRMSE of the augmented score, the value-added ratio
#'   (VAR), the difference with its 95% interval and decision, and the
#'   augmented gain over the better of the two with its interval. When a PRMSE
#'   falls outside (0, 1], usually because a section is too short or
#'   unreliable, the decision is "not estimable" and a warning is given.
#' @references
#' Haberman, S. J. (2008). When can subscores have value? *Journal of
#' Educational and Behavioral Statistics, 33*(2), 204-229.
#'
#' Sinharay, S. (2019). Added value of subscores and hypothesis testing.
#' *Journal of Educational and Behavioral Statistics, 44*(1), 25-44.
#'
#' Feldt, L. S., & Charter, R. A. (2003). Estimating the reliability of a test
#' split into two parts of equal or unequal length. *Psychological Methods,
#' 8*(1), 102-109.
#' @seealso [profile_check()]
#' @examples
#' data(sim_test)
#' X <- sim_test[, sim_test_design$items]
#' keep <- screen_records(X, sim_test_design$sections)$flags$keep
#' added_value(X[keep, ], sim_test_design$sections, B = 200)
#' @export
added_value <- function(items, sections, reliability = c("kr20", "part"), parts = NULL,
                        B = 1000, seed = 2025) {
  reliability <- match.arg(reliability)
  X <- as_item_matrix(items)
  sections <- resolve_sections(sections, colnames(X))
  K <- names(sections)
  if (reliability == "part") {
    if (is.null(parts) || !all(K %in% names(parts))) stop("`parts` must give the parts of every section.", call. = FALSE)
    parts <- lapply(parts[K], resolve_sections_unnamed, cn = colnames(X))
  }
  # score matrix: one column per section, plus part scores when needed
  S <- vapply(sections, function(i) rowSums(X[, i, drop = FALSE]), numeric(nrow(X)))
  P <- if (reliability == "part") do.call(cbind, lapply(K, function(k)
    vapply(parts[[k]], function(i) rowSums(X[, i, drop = FALSE]), numeric(nrow(X))))) else NULL
  pcol <- if (reliability == "part") split(seq_len(ncol(P)), rep(K, lengths(parts[K])))[K] else NULL
  n <- nrow(X)
  stat <- function(w) {
    sw <- sum(w)
    wcov <- function(M) { m <- colSums(M * w) / sw; (crossprod(M * w, M) / sw - tcrossprod(m)) * sw / (sw - 1) }
    C <- wcov(S)
    p <- colSums(X * w) / sw; iv <- p * (1 - p) * sw / (sw - 1)
    CP <- if (!is.null(P)) wcov(P) else NULL
    t(vapply(seq_along(K), function(j) {
      vs <- C[j, j]
      rel <- if (reliability == "kr20") {
        kk <- length(sections[[j]]); kk / (kk - 1) * (1 - sum(iv[sections[[j]]]) / vs)
      } else {
        Cp <- CP[pcol[[j]], pcol[[j]], drop = FALSE]; lam <- rowSums(Cp) / sum(Cp)
        (sum(Cp) - sum(diag(Cp))) / (sum(Cp) * (1 - sum(lam^2)))
      }
      vo <- sum(C[-j, -j]); cso <- sum(C[j, -j]); vx <- vs + vo + 2 * cso
      vt <- rel * vs; prx <- (vt + cso)^2 / (vt * vx)
      Sig <- matrix(c(vs, cso, cso, vo), 2); cv <- c(vt, cso)
      praug <- as.numeric(crossprod(cv, solve(Sig, cv))) / vt
      c(prmse_s = rel, prmse_x = prx, prmse_aug = praug)
    }, numeric(3)))
  }
  est <- stat(rep(1, n))
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  set.seed(seed)
  bt <- replicate(B, {
    s <- stat(tabulate(sample.int(n, n, replace = TRUE), n))
    cbind(diff = s[, 1] - s[, 2], gain = s[, 3] - pmax(s[, 1], s[, 2]))
  }, simplify = "array")
  q <- function(m) t(apply(m, 1, stats::quantile, c(.025, .975), na.rm = TRUE))
  qd <- q(matrix(bt[, 1, ], nrow = length(K))); qg <- q(matrix(bt[, 2, ], nrow = length(K)))
  out <- data.frame(section = K, prmse_s = est[, 1], prmse_x = est[, 2], prmse_aug = est[, 3],
                    VAR = est[, 1] / est[, 2], diff = est[, 1] - est[, 2],
                    diff_low = qd[, 1], diff_high = qd[, 2],
                    decision = ifelse(qd[, 1] > 0, "adds value", ifelse(qd[, 2] < 0, "no added value", "inconclusive")),
                    aug_gain = est[, 3] - pmax(est[, 1], est[, 2]), gain_low = qg[, 1], gain_high = qg[, 2],
                    row.names = NULL)
  bad <- !(out$prmse_s > 0 & out$prmse_s <= 1 & out$prmse_x > 0 & out$prmse_x <= 1)
  if (any(bad)) {
    out$decision[bad] <- "not estimable"
    warning("PRMSE outside (0, 1] for ", paste(out$section[bad], collapse = ", "),
            " (section reliability too low or unstable); the decision is reported as 'not estimable'.",
            call. = FALSE)
  }
  attr(out, "reliability") <- reliability; attr(out, "B") <- B; attr(out, "n") <- n
  class(out) <- c("fc_added_value", "data.frame")
  out
}

#' @export
print.fc_added_value <- function(x, digits = 3, ...) {
  if (is.null(attr(x, "reliability")) || !all(c("section", "decision", "aug_gain") %in% names(x)))
    return(print(as.data.frame(unclass(x)), digits = digits))
  y <- as.data.frame(unclass(x)); num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], function(v) round(v, digits))
  cat("Added value of section scores (Haberman, 2008)\n")
  cat("Reliability:", if (attr(x, "reliability") == "kr20") "KR-20 (item level)" else "Feldt-Raju (part level)",
      "| N =", attr(x, "n"), "| bootstrap resamples =", attr(x, "B"), "\n\n")
  print(y[, c("section", "prmse_s", "prmse_x", "VAR", "diff", "diff_low", "diff_high", "decision")], row.names = FALSE)
  cat("\nAugmented score gain over the better of the two:\n")
  print(y[, c("section", "prmse_aug", "aug_gain", "gain_low", "gain_high")], row.names = FALSE)
  if (attr(x, "B") < 1000) cat("\nUse B >= 1000 for reported intervals.\n")
  invisible(x)
}

#' Can two section scores be told apart? Reliability of score differences
#'
#' For each pair of sections, the reliability of the difference between the
#' standardized section scores, rho_D = ((rho_a + rho_b)/2 - r_ab) / (1 - r_ab),
#' its standard error, SE_D = sqrt(2 - rho_a - rho_b), and the percentage of
#' records whose absolute difference exceeds 1.96 SE_D. Under a normal-error
#' approximation and with no true difference, about 5% would exceed it by
#' chance for one pair. Comparisons are not adjusted for multiplicity.
#'
#' @inheritParams added_value
#' @return A data frame with one row per pair of sections.
#' @examples
#' data(sim_test)
#' X <- sim_test[, sim_test_design$items]
#' keep <- screen_records(X, sim_test_design$sections)$flags$keep
#' profile_check(X[keep, ], sim_test_design$sections)
#' @export
profile_check <- function(items, sections) {
  X <- as_item_matrix(items)
  sections <- resolve_sections(sections, colnames(X))
  S <- vapply(sections, function(i) rowSums(X[, i, drop = FALSE]), numeric(nrow(X)))
  rel <- vapply(names(sections), function(k) {
    it <- sections[[k]]; kk <- length(it); p <- colMeans(X[, it, drop = FALSE])
    kk / (kk - 1) * (1 - sum(p * (1 - p)) * nrow(X) / (nrow(X) - 1) / stats::var(S[, k]))
  }, numeric(1))
  Z <- scale(S)
  pr <- utils::combn(names(sections), 2)
  out <- t(apply(pr, 2, function(ab) {
    a <- ab[1]; b <- ab[2]; r <- stats::cor(S[, a], S[, b])
    se <- sqrt(2 - rel[[a]] - rel[[b]])
    c(rel_a = rel[[a]], rel_b = rel[[b]], r = r,
      reliability_diff = ((rel[[a]] + rel[[b]]) / 2 - r) / (1 - r), se_diff = se,
      pct_beyond = 100 * mean(abs(Z[, a] - Z[, b]) > 1.96 * se))
  }))
  data.frame(pair = paste(pr[1, ], pr[2, ], sep = " - "), round(out, 3), row.names = NULL)
}

resolve_sections_unnamed <- function(lst, cn) {
  if (is.null(names(lst))) names(lst) <- paste0("part", seq_along(lst))
  resolve_sections(lst, cn)
}
