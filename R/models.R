#' lavaan syntax for a second-order factor model
#'
#' @param blocks A named or unnamed list of indicator names, one element per
#'   first-order factor. Names become the factor names (default `F1`, `F2`, ...).
#' @param general Name of the second-order factor.
#' @return A character string of lavaan model syntax. The first loading of each
#'   first-order factor is fixed to 1 (marker), and the second-order factor
#'   variance is fixed to 1 with all its loadings free. Fit it with
#'   `std.lv = TRUE` to fix all latent variances instead.
#' @examples
#' cat(second_order_syntax(list(L = c("l1", "l2", "l3"), R = c("r1", "r2", "r3"),
#'                              G = c("g1", "g2"))))
#' @export
second_order_syntax <- function(blocks, general = "G") {
  fo <- names(blocks)
  if (is.null(fo) || any(fo == "")) fo <- paste0("F", seq_along(blocks))
  if (length(blocks) < 3) stop("A second-order model needs at least three first-order factors.", call. = FALSE)
  lines <- mapply(function(f, v) paste(f, "=~", paste(v, collapse = " + ")), fo, blocks)
  top <- paste0(general, " =~ NA*", fo[1], " + ", paste(fo[-1], collapse = " + "))
  paste(c(lines, top, paste0(general, " ~~ 1*", general)), collapse = "\n")
}

#' Population correlation matrix of a second-order model
#'
#' @param sizes Number of indicators per first-order factor, e.g. `c(3, 2, 5)`.
#' @param lambda First-order standardized loadings (one value, or one per indicator).
#' @param gamma Second-order standardized loadings (one value, or one per factor).
#' @param cross Optional omitted cross-loading:
#'   `list(item = 3, factor = 3, value = .30)`.
#' @return A correlation matrix with indicators named `x1`, `x2`, ...
#' @examples
#' S <- make_population(c(3, 2, 5), lambda = .70, gamma = c(.95, .90, .85))
#' round(S[1:4, 1:4], 3)
#' @export
make_population <- function(sizes, lambda, gamma, cross = NULL) {
  p <- sum(sizes); k <- length(sizes); fac <- rep(seq_len(k), sizes)
  lambda <- rep_len(lambda, p); gamma <- rep_len(gamma, k)
  if (any(abs(gamma) >= 1)) stop("Second-order loadings must be below 1 in absolute value.", call. = FALSE)
  L <- matrix(0, p, k)
  L[cbind(seq_len(p), fac)] <- lambda
  if (!is.null(cross)) L[cross$item, cross$factor] <- cross$value
  Phi <- tcrossprod(gamma); diag(Phi) <- 1
  S <- L %*% Phi %*% t(L)
  if (any(diag(S) >= 1)) stop("Communality of 1 or more: lower the loadings.", call. = FALSE)
  diag(S) <- 1
  nm <- paste0("x", seq_len(p)); dimnames(S) <- list(nm, nm)
  S
}

#' Indicator blocks from block sizes
#'
#' @param sizes Number of indicators per factor.
#' @return A list of indicator names (`x1`, `x2`, ...) per factor, for
#'   [second_order_syntax()] and [make_population()].
#' @examples
#' blocks_from_sizes(c(3, 2, 5))
#' @export
blocks_from_sizes <- function(sizes) {
  idx <- split(seq_len(sum(sizes)), rep(seq_along(sizes), sizes))
  names(idx) <- paste0("F", seq_along(sizes))
  lapply(idx, function(i) paste0("x", i))
}
