# Demonstration of the fitcheck workflow on the simulated data shipped with
# the package (Psychoco 2027). Runs in about one minute on two cores.
library(fitcheck)
data(sim_test); data(sim_test_design)
items <- sim_test[, sim_test_design$items]

# 1. The whole workflow: Check 0, screening sensitivity (overall, by form, and
#    invariance across forms), Checks 1-4, score reporting
report <- fitcheck(items, sim_test_design$sections, sim_test_design$parcels,
                   form = sim_test$form, B = 1000)
print(report)

# 2. Sample-size behaviour of the retained data (identical on one or two cores)
#    40 replications per sample size keep this demonstration fast; they illustrate the
#    output and are not enough for precise estimates of admissibility rates.
keep <- report$screen$flags$keep
d <- parcel_scores(items[keep, ], sim_test_design$parcels)
print(boot_by_n(d, sim_test_design$model, ns = c(100, 300, 1000), reps = 40, cores = 2, seed = 11))

# 3. Does screening change the reporting decision? Keep the whole-section zeros
noall <- !report$screen$flags$all_zero
cat("\nAlternative screening: remove only all-zero records (keep whole-section-zero records).\n",
    "N =", sum(noall), "instead of", sum(report$screen$flags$keep), "\n")
print(added_value(items[noall, ], sim_test_design$sections, B = 1000))

# 4. Example 2: flagged records can change the comparison of forms.
# Three forms with the same 27 items (three sections of 9 items, three parcels
# per section, easy, medium and hard), the same loadings in every form, and 100
# extra all-zero records added to form B only.
make_form <- function(n, seed) {
  set.seed(seed)
  p_correct <- rep(rep(c(.85, .60, .35), each = 3), 3)
  R <- matrix(.75, 3, 3); diag(R) <- 1
  F <- matrix(stats::rnorm(n * 3), n) %*% chol(R)
  Y <- sapply(1:27, function(j) {
    s <- (j - 1) %/% 9 + 1
    as.integer(.65 * F[, s] + sqrt(1 - .65^2) * stats::rnorm(n) > stats::qnorm(1 - p_correct[j]))
  })
  colnames(Y) <- paste0(rep(c("L", "G", "R"), each = 9), rep(1:9, 3))
  Y
}
X2 <- rbind(make_form(700, 1), make_form(700, 2), make_form(700, 3),
            matrix(0L, 100, 27))
form2 <- c(rep(c("A", "B", "C"), each = 700), rep("B", 100))
sections2 <- list(Listening = paste0("L", 1:9), Grammar = paste0("G", 1:9),
                  Reading = paste0("R", 1:9))
parcels2 <- list()
for (s in c("L", "G", "R")) for (p in 1:3) parcels2[[paste0(s, p)]] <- paste0(s, (p - 1) * 3 + 1:3)
model2 <- second_order_syntax(list(Listening = c("L1", "L2", "L3"),
                                   Grammar = c("G1", "G2", "G3"),
                                   Reading = c("R1", "R2", "R3")))
s2 <- screen_records(X2, sections2, form = form2)
print(screening_invariance(X2, s2$flags$keep, form2, parcels2, model2))
