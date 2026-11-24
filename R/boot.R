#' How would the model behave in smaller samples of your own data?
#'
#' Draws `reps` bootstrap samples (with replacement) at each sample size,
#' refits the model, and summarizes how often estimation succeeds and how the
#' fit indices behave. Every replication has its own pre-drawn seed, so results
#' are identical with one or several cores, and lavaan cannot disturb the
#' random-number stream.
#'
#' Sizes above the number of records reuse records, so treat them as
#' approximations.
#'
#' @param data A data frame of indicators.
#' @param model lavaan model syntax.
#' @param ns Sample sizes.
#' @param reps Replications per size.
#' @param seed Random seed.
#' @param cores Number of CPU cores. Values above 1 use a PSOCK cluster from
#'   the \pkg{parallel} package.
#' @param std.lv Fix latent variances to 1? Default `FALSE` (marker
#'   variables), which shows Heywood cases as improper solutions rather than as
#'   nonconvergence.
#' @param estimator Estimator for [lavaan::cfa()]. Default `"MLR"`.
#' @return An object of class `fc_runs`: a list with `raw` (one row per
#'   replication) and `summary` (one row per size: admissible percentage with
#'   a 95% Wilson interval, counts of improper and nonconverged solutions, and
#'   mean fit indices over admissible solutions; robust CFI and RMSEA can be
#'   unavailable in some small resamples, so their means use the solutions in
#'   which they could be computed, and `cfi_missing` counts the others).
#' @seealso [mini_mc()]
#' @examples
#' \donttest{
#' data(sim_test)
#' keep <- screen_records(sim_test[, sim_test_design$items], sim_test_design$sections)$flags$keep
#' d <- parcel_scores(sim_test[keep, sim_test_design$items], sim_test_design$parcels)
#' boot_by_n(d, sim_test_design$model, ns = c(100, 300), reps = 20)
#' }
#' @export
boot_by_n <- function(data, model, ns = c(100, 200, 300, 500, 1000), reps = 50,
                      seed = 2025, cores = 1, std.lv = FALSE, estimator = "MLR") {
  data <- as.data.frame(data)
  force(model); force(std.lv); force(estimator)   # evaluate before workers receive the closure
  tasks <- make_tasks(ns, reps, seed)
  worker <- function(t) {
    set.seed(t$seed)
    d <- data[sample.int(nrow(data), t$N, replace = TRUE), , drop = FALSE]
    c(N = t$N, rep = t$rep, fit_once(model, d, std.lv = std.lv, estimator = estimator))
  }
  finish_runs(run_tasks(tasks, worker, cores), "bootstrap", nrow(data))
}

#' A small Monte Carlo check from a population correlation matrix
#'
#' Simulates multivariate normal data from `Sigma` at each sample size, fits
#' the model, and summarizes estimation success and fit. `Sigma` can come from
#' [make_population()] or from a fitted model:
#' `cov2cor(lavaan::lavInspect(fit, "implied")$cov)`.
#'
#' @inheritParams boot_by_n
#' @param Sigma A population correlation (or covariance) matrix with
#'   dimnames matching the indicator names in `model`.
#' @return An object of class `fc_runs` (see [boot_by_n()]).
#' @examples
#' \donttest{
#' S <- make_population(c(3, 2, 5), lambda = .50, gamma = c(.98, .95, .95))
#' mini_mc(S, second_order_syntax(blocks_from_sizes(c(3, 2, 5))),
#'         ns = c(150, 300), reps = 20)
#' }
#' @export
mini_mc <- function(Sigma, model, ns = c(100, 200, 300, 500, 1000), reps = 100,
                    seed = 2025, cores = 1, std.lv = FALSE, estimator = "MLR") {
  if (is.null(colnames(Sigma))) stop("`Sigma` needs column names that match the model.", call. = FALSE)
  force(model); force(std.lv); force(estimator)
  tasks <- make_tasks(ns, reps, seed)
  worker <- function(t) {
    set.seed(t$seed)
    d <- as.data.frame(MASS::mvrnorm(t$N, rep(0, ncol(Sigma)), Sigma))
    c(N = t$N, rep = t$rep, fit_once(model, d, std.lv = std.lv, estimator = estimator))
  }
  finish_runs(run_tasks(tasks, worker, cores), "Monte Carlo", NA)
}

#' @export
print.fc_runs <- function(x, ...) {
  cat(x$type, "runs:", nrow(x$raw), "replications\n")
  if (!is.na(x$n_data) && any(x$summary$N > x$n_data))
    cat("Note: sizes above", x$n_data, "reuse records.\n")
  s <- x$summary; num <- vapply(s, is.numeric, logical(1))
  s[num] <- lapply(s[num], function(v) round(v, 4))
  print(s, row.names = FALSE)
  invisible(x)
}

# ---- internal helpers ----------------------------------------------------
fit_once <- function(model, data, std.lv = FALSE, estimator = "MLR") {
  fit <- tryCatch(suppressWarnings(lavaan::cfa(model, data = data, estimator = estimator, std.lv = std.lv)),
                  error = function(e) NULL)
  res <- c(status = "error", srmr = NA, usrmr = NA, cfi = NA, rmsea = NA)
  if (is.null(fit)) return(res)
  if (!isTRUE(lavaan::lavInspect(fit, "converged"))) { res["status"] <- "nonconverged"; return(res) }
  ok <- isTRUE(suppressWarnings(lavaan::lavInspect(fit, "post.check")))
  res["status"] <- if (ok) "admissible" else "improper"
  fm <- lavaan::fitMeasures(fit, c("srmr", "cfi.robust", "rmsea.robust"))
  res[c("srmr", "cfi", "rmsea")] <- fm[c("srmr", "cfi.robust", "rmsea.robust")]
  res["usrmr"] <- usrmr_summary(fit)$usrmr
  res
}

make_tasks <- function(ns, reps, seed) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  set.seed(seed)
  g <- expand.grid(rep = seq_len(reps), N = ns)
  g$seed <- sample.int(.Machine$integer.max, nrow(g))
  lapply(seq_len(nrow(g)), function(i) as.list(g[i, ]))
}

run_tasks <- function(tasks, worker, cores) {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    old <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit(assign(".Random.seed", old, envir = globalenv()))
  }
  if (cores > 1) {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, { requireNamespace("lavaan"); requireNamespace("fitcheck") })
    res <- parallel::parLapplyLB(cl, tasks, worker)
  } else res <- lapply(tasks, worker)
  do.call(rbind, res)
}

finish_runs <- function(raw, type, n_data) {
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)
  for (v in c("N", "rep", "srmr", "usrmr", "cfi", "rmsea")) raw[[v]] <- as.numeric(raw[[v]])
  Ns <- sort(unique(raw$N))
  summ <- do.call(rbind, lapply(Ns, function(n) {
    d <- raw[raw$N == n, ]; a <- d[d$status == "admissible", ]
    k <- nrow(a); m <- nrow(d); ci <- wilson(k, m)
    data.frame(N = n, reps = m, admissible_pct = 100 * k / m,
               ci_low = 100 * ci[1], ci_high = 100 * ci[2],
               improper = sum(d$status == "improper"),
               nonconverged = sum(d$status %in% c("nonconverged", "error")),
               srmr = mean(a$srmr, na.rm = TRUE), usrmr = mean(a$usrmr, na.rm = TRUE),
               usrmr_zero_pct = 100 * mean(a$usrmr <= 1e-12, na.rm = TRUE),
               cfi = mean(a$cfi, na.rm = TRUE), rmsea = mean(a$rmsea, na.rm = TRUE),
               cfi_missing = sum(is.na(a$cfi)))
  }))
  num <- vapply(summ, is.numeric, logical(1))
  summ[num] <- lapply(summ[num], function(v) replace(v, is.nan(v), NA))
  none <- summ$N[summ$admissible_pct == 0]
  if (length(none))
    warning("No admissible solution at N = ", paste(none, collapse = ", "),
            "; fit summaries for these sizes are NA.", call. = FALSE)
  structure(list(raw = raw, summary = summ, type = type, n_data = n_data), class = "fc_runs")
}

wilson <- function(k, n, z = stats::qnorm(.975)) {
  if (n == 0) return(c(NA, NA))
  p <- k / n; den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(max(0, ctr - hw), min(1, ctr + hw))
}
