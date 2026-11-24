# Tests for fitcheck (base R only, run by R CMD check)
library(fitcheck)
ok <- function(cond, msg) if (!isTRUE(cond)) stop("FAILED: ", msg, call. = FALSE) else cat("ok:", msg, "\n")
near <- function(a, b, tol = 1e-8) isTRUE(all.equal(as.vector(a), as.vector(b), tolerance = tol))

# ---- Check 0 ------------------------------------------------------------
X <- matrix(1L, 6, 6, dimnames = list(NULL, paste0("i", 1:6)))
X[1, ] <- 0L            # all-zero record
X[2, 1:3] <- 0L         # whole first section zero
X[3, 1] <- 0L           # ordinary wrong answer
secs <- list(A = 1:3, B = 4:6)
s <- screen_records(X, secs, form = c("f1", "f1", "f1", "f2", "f2", "f2"), max_p_zero = 1)
ok(sum(s$flags$all_zero) == 1, "one all-zero record")
ok(sum(s$flags$section_zero) == 1 && s$flags$section_zero[2], "one whole-section zero record")
ok(identical(s$flags$keep, c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)), "keep flags")
ok(s$by_form$flagged_pct[s$by_form$form == "f1"] == round(100 * 2 / 3, 1), "percentage by form")
ok(s$low_scores$score_0[s$low_scores$section == "A"] == 1, "low-score table counts section zeros")
ok(near(zero_probability(50), .75^50), "guessing probability")
s2 <- screen_records(as.data.frame(X), list(A = c("i1", "i2", "i3"), B = c("i4", "i5", "i6")), max_p_zero = 1)
ok(identical(s$flags$keep, s2$flags$keep), "sections given by name")
err <- tryCatch(screen_records(X, list(1:3)), error = function(e) "error")
ok(identical(err, "error"), "unnamed sections are refused")

# ---- models -------------------------------------------------------------
S <- make_population(c(3, 2, 5), lambda = .70, gamma = c(.95, .90, .85))
ok(all(diag(S) == 1) && isSymmetric(S), "population matrix is a correlation matrix")
ok(near(S["x1", "x4"], .70 * .95 * .90 * .70), "population correlation across factors")
ok(identical(tryCatch(make_population(c(3, 3, 3), .7, 1), error = function(e) "e"), "e"), "gamma of 1 refused")
syn <- second_order_syntax(blocks_from_sizes(c(3, 2, 5)))
ok(grepl("G =~ NA\\*F1 \\+ F2 \\+ F3", syn), "second-order syntax")

# ---- Checks 1 to 4: reproduces handbook Example 1 (improper solution) ------
S1 <- make_population(c(3, 2, 5), lambda = .50, gamma = .95)
set.seed(101); d1 <- as.data.frame(MASS::mvrnorm(150, rep(0, 10), S1))
fit1 <- lavaan::cfa(syn, data = d1, estimator = "MLR")
c1 <- four_checks(fit1, verbose = FALSE)
ok(c1$converged && !c1$admissible, "Example 1 is converged but not admissible")
ok("F3" %in% c1$neg_disturbance, "Example 1 has a negative disturbance variance for F3")
# Exact values depend on the linear-algebra library (handbook run: CFI .959, SRMR .050),
# so the test checks the behaviour: an improper solution with good conventional fit.
ok(c1$cfi > .93 && c1$srmr < .08, "Example 1 is improper but has good conventional fit")
idc <- identification_check(syn, d1)
ok(nrow(idc) == 2 && !idc$admissible[2], "identification check runs and flags the marker solution")

set.seed(7); d <- as.data.frame(MASS::mvrnorm(1500, rep(0, 10), S))
fit <- lavaan::cfa(syn, data = d, estimator = "MLR")
c2 <- four_checks(fit, verbose = FALSE)
ok(c2$admissible, "a well-estimated model is admissible")
h <- hierarchical_indices(fit)
ok(h$omega_h < h$omega["total"] && h$ecv > 0 && h$ecv < 1, "omega_h below omega total; ECV in (0, 1)")
ok(all(h$omega_hs > 0 & h$omega_hs < h$omega[names(h$omega_hs)]), "omega_hs below subscale omega")

# ---- added value: compare with direct formulas -----------------------------
data(sim_test, package = "fitcheck"); data(sim_test_design, package = "fitcheck")
Xs <- as.matrix(sim_test[, sim_test_design$items])
keep <- screen_records(Xs, sim_test_design$sections)$flags$keep
Xk <- Xs[keep, ]
sec <- lapply(sim_test_design$sections, function(v) match(v, colnames(Xk)))
set.seed(99); before <- runif(1)
set.seed(99)
av <- added_value(Xk, sec, B = 50, seed = 1)
ok(near(runif(1), before), "added_value() leaves the user's random-number stream unchanged")
Ssum <- sapply(sec, function(i) rowSums(Xk[, i]))
s1 <- Ssum[, 1]; o1 <- rowSums(Ssum[, -1]); k1 <- length(sec[[1]]); p1 <- colMeans(Xk[, sec[[1]]])
n <- nrow(Xk)
kr <- k1 / (k1 - 1) * (1 - sum(p1 * (1 - p1)) * n / (n - 1) / var(s1))
vt <- kr * var(s1); prx <- (vt + cov(s1, o1))^2 / (vt * var(s1 + o1))
ok(near(av$prmse_s[1], kr, 1e-10), "PRMSE_s equals KR-20")
ok(near(av$prmse_x[1], prx, 1e-10), "PRMSE_x matches the Haberman formula")
ok(all(av$prmse_aug >= pmax(av$prmse_s, av$prmse_x) - 1e-10), "augmented score is never worse than the better of the two")
av2 <- added_value(Xk, sec, B = 50, seed = 1)
ok(identical(av$diff_low, av2$diff_low), "bootstrap intervals are reproducible with a seed")
ok(all(av$decision %in% c("adds value", "no added value", "inconclusive")), "decisions are valid")
parts <- list(Listening = list(1:10, 11:20, 21:30), Grammar = list(31:40, 41:50),
              Reading = list(51:56, 57:62, 63:68, 69:74, 75:80))
avp <- added_value(Xk, sec, reliability = "part", parts = parts, B = 20)
ok(all(avp$prmse_s <= av$prmse_s + .02), "part-level reliability is not above KR-20 by more than .02 here")

pc <- profile_check(Xk, sec)
r12 <- cor(Ssum[, 1], Ssum[, 2])
k2 <- length(sec[[2]]); p2 <- colMeans(Xk[, sec[[2]]])
kr2 <- k2 / (k2 - 1) * (1 - sum(p2 * (1 - p2)) * n / (n - 1) / var(Ssum[, 2]))
ok(abs(pc$reliability_diff[1] - ((kr + kr2) / 2 - r12) / (1 - r12)) < .0006, "reliability of differences")
ok(abs(pc$se_diff[1] - sqrt(2 - kr - kr2)) < .0006, "standard error of differences")

# ---- screening sensitivity and the workflow -------------------------------
ss <- screening_sensitivity(Xs, keep, sim_test_design$parcels, sim_test_design$model)
ok(all(c("CFI", "SRMR", "omega_h", "ECV") %in% ss$index), "sensitivity table has fit and measurement rows")
rep1 <- fitcheck(Xs, sim_test_design$sections, sim_test_design$parcels, form = sim_test$form, B = 20)
ok(inherits(rep1, "fc_report") && rep1$checks$admissible, "workflow runs and the retained solution is admissible")
auto <- fitcheck(Xs, sim_test_design$sections, sim_test_design$parcels, subscores = FALSE)
ok(grepl("Listening =~ L1 \\+ L2 \\+ L3", auto$model), "model is built from sections and parcels")

# ---- sample-size runs are reproducible -------------------------------------
dk <- parcel_scores(Xk, sim_test_design$parcels)
b1 <- boot_by_n(dk, sim_test_design$model, ns = c(150, 300), reps = 3, seed = 5)
b2 <- boot_by_n(dk, sim_test_design$model, ns = c(150, 300), reps = 3, seed = 5)
ok(identical(b1$raw, b2$raw), "boot_by_n() is reproducible")
ok(all(b1$summary$ci_low <= b1$summary$admissible_pct & b1$summary$admissible_pct <= b1$summary$ci_high),
   "Wilson interval contains the estimate")
m1 <- mini_mc(S1, syn, ns = 150, reps = 4, seed = 3)
ok(nrow(m1$raw) == 4 && sum(m1$summary[, c("improper", "nonconverged")]) + m1$summary$admissible_pct * 4 / 100 == 4,
   "mini_mc() classifies every replication")
b3 <- boot_by_n(dk, sim_test_design$model, ns = c(150, 300), reps = 3, seed = 5, cores = 2)
ok(identical(b1$raw, b3$raw), "boot_by_n() gives identical results on two cores")
# ---- screening by form and invariance across forms -------------------------
fo <- fitcheck:::first_order_model(sim_test_design$model)
ok(!grepl("G =~", fo) && grepl("Reading =~ R1 \\+ R2 \\+ R3 \\+ R4 \\+ R5", fo),
   "second-order model is turned into the correlated first-order model")
s_all <- screen_records(Xs, sim_test_design$sections, form = sim_test$form)
inv <- screening_invariance(Xs, s_all$flags$keep, sim_test$form, sim_test_design$parcels, sim_test_design$model)
ok(inherits(inv, "fc_invariance") && nrow(inv$table) == 6, "invariance table has three levels for two data versions")
ok(all(inv$table$converged), "all invariance models converge")
ok(all(inv$verdict %in% c("configural", "metric", "scalar")), "verdict is an invariance level")
ok(identical(inv$same_conclusion, unname(inv$verdict[1]) == unname(inv$verdict[2])), "same_conclusion matches the verdicts")
t1 <- inv$table[inv$table$version == "retained", ]
ok(near(t1$delta_CFI[2], t1$CFI[2] - t1$CFI[1]), "delta CFI is the change from the previous level")
fake <- inv$table; fake$CFI <- c(1, .98, .97, 1, .995, .99); fake$RMSEA <- c(.02, .04, .06, .02, .03, .04)
ok(identical(fitcheck:::invariance_decision(fake$CFI[1:3], fake$RMSEA[1:3], .01, .015),
             c("baseline", "not supported", "not supported")) &&
   identical(fitcheck:::invariance_decision(fake$CFI[4:6], fake$RMSEA[4:6], .01, .015),
             c("baseline", "supported", "supported")), "Chen criteria: both changes must exceed the cutoffs")
err2 <- tryCatch(screening_invariance(Xs, s_all$flags$keep, rep("A", nrow(Xs)), sim_test_design$parcels,
                                      sim_test_design$model), error = function(e) "error")
ok(identical(err2, "error"), "one form only is refused")
bf <- screening_by_form(Xs, s_all$flags$keep, sim_test$form, sim_test_design$parcels, sim_test_design$model)
ml <- bf[bf$index == "mean_first_order_loading", ]
ok(setequal(unique(bf$form), c("A", "B", "C")), "screening_by_form() reports every form")
ok(ml$form[which.max(abs(ml$difference))] == "B", "the form with most flagged records changes most")
ok(inherits(rep1$invariance, "fc_invariance") && inherits(rep1$by_form, "fc_by_form"),
   "workflow includes the by-form and invariance checks when forms are given")
ok(is.null(auto$invariance), "no invariance check without forms")
# ---- review fixes ---------------------------------------------------------------
bad <- Xs[1:5, ]; bad[1, 1] <- 2
err3 <- tryCatch(screen_records(bad, sim_test_design$sections), error = function(e) conditionMessage(e))
ok(grepl("dichotomous", err3), "non-0/1 scores are refused")
idc <- identification_check(sim_test_design$model, parcel_scores(Xs[s_all$flags$keep, ], sim_test_design$parcels))
ok(all(idc$negative_disturbance == "none") && all(idc$negative_residual == "none"),
   "identification table prints none instead of blanks")
ok("admissible" %in% names(inv$table) && all(inv$table$admissible), "invariance table reports admissibility")
ok(attr(fo, "n_first") == 3 && attr(fo, "n_higher") == 1, "first-order model records the factor counts")
m2 <- "L =~ L1 + L2 + L3\nR =~ R1 + R2 + R3 + R4 + R5\nG =~ L + R"
w <- NULL
withCallingHandlers(tryCatch(screening_invariance(Xs, s_all$flags$keep, sim_test$form, sim_test_design$parcels, m2),
                             error = function(e) NULL),
                    warning = function(cw) { w <<- conditionMessage(cw); invokeRestart("muffleWarning") })
ok(!is.null(w) && grepl("three first-order factors", w), "warns when the model does not have three first-order factors")
sens2 <- screening_sensitivity(Xs, s_all$flags$keep, sim_test_design$parcels, sim_test_design$model, std.lv = FALSE)
adm2 <- sens2[sens2$index == "admissible", ]
ok(all(c(adm2$all_records, adm2$retained) == 1) && "ECV" %in% sens2$index,
   "omega and ECV are reported for admissible solutions")
# ---- example 2: screening changes the invariance conclusion ------------------
# Example 2: flagged records can change the comparison of forms.
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
inv2 <- screening_invariance(X2, s2$flags$keep, form2, parcels2, model2)
ok(!inv2$same_conclusion && inv2$verdict[["all records"]] == "configural" && inv2$verdict[["retained"]] == "scalar",
   "all-zero records in one form can create apparent non-invariance")
# ---- short sections: whole-section zeros are not flagged -------------------
Xshort <- matrix(c(0,0,0,0, 1,1,0,1,  1,0,1,1, 0,1,1,0,  1,1,1,1, 1,1,1,1), nrow = 3, byrow = TRUE,
                 dimnames = list(NULL, paste0("i", 1:8)))
ss <- screen_records(Xshort, list(A = 1:4, B = 5:8))
ok(!any(ss$flags$section_zero) && !any(ss$low_scores$flagged), "4-item sections: zeros are not flagged by default")
ss1 <- screen_records(Xshort, list(A = 1:4, B = 5:8), max_p_zero = 1)
ok(ss1$flags$section_zero[1] && sum(ss1$flags$section_zero) == 1, "max_p_zero = 1 flags whole-section zeros in every section")
ok(all(s_all$low_scores$flagged), "long sections are flagged as before")
cat("All tests passed.\n")
