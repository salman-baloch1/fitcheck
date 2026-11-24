# fitcheck 0.2.0

* New `screening_by_form()`: repeats the screening sensitivity analysis within
  each test form, because flagged records are often concentrated in a few forms.
* New `screening_invariance()`: tests configural, metric, and scalar invariance
  across forms with all records and with the retained records, and reports
  whether screening changes the conclusion (Chen, 2007, criteria).
* `fitcheck()` runs both checks when records come from two or more forms
  (new argument `invariance`, default `TRUE`, added at the end of the
  argument list so existing calls are unchanged).
* `screening_sensitivity()` computes omega, omega-hierarchical and ECV only
  from admissible solutions, and warns when either solution is improper.
* `screening_invariance()` reports admissibility for every model, prints its
  decision rule, and warns when the model does not have three first-order
  factors (the case in which the correlated and second-order models are
  equivalent).
* `identification_check()` prints "none" instead of blank cells when no
  negative variance is found.
* `screen_records()` flags whole-section zeros only in sections where a zero
  score is unusual under random guessing (new arguments `n_options = 4` and
  `max_p_zero = .01`). In short sections a zero score is common among genuine
  low scorers. Results for long sections, such as those of the example data,
  do not change. `fitcheck()` passes both arguments on.
* `added_value()` reports "not estimable" with a warning when a PRMSE falls
  outside (0, 1], for example for a very short or unreliable section.
* Item scores must be 0 or 1; other values now stop with a clear message.
* `boot_by_n()` and `mini_mc()` average robust CFI and RMSEA over the
  admissible solutions in which they could be computed (they can be
  unavailable in small resamples), and report the number missing in
  `cfi_missing`; previously one missing value made the mean `NA`.
* `boot_by_n()` and `mini_mc()` warn when a sample size has no admissible
  solution, and report its fit summaries as `NA`.
* Vignette, demo and tests include a second simulated example in which
  all-zero records concentrated in one form make metric invariance fail,
  while the retained records support scalar invariance.
* New script `inst/examples/public_data_ability.R` runs the workflow on the
  public ICAR `ability` data of the psychTools package (not a dependency;
  the script installs it if needed).
* README and help pages: the evidence now comes from the doctoral study behind
  the package (Rashid, 2026), with numbers that can be reproduced from its
  analysis script.

# fitcheck 0.1.0

* First release. Functions for Check 0 (screen_records(), zero_probability(),
  screening_sensitivity()), Checks 1 to 4 (four_checks(),
  identification_check(), hierarchical_indices()), sample-size behaviour
  (boot_by_n(), mini_mc()), score reporting (added_value(), profile_check()),
  and the whole workflow (fitcheck()). Simulated example data: sim_test.
* Vignette "Getting started with fitcheck".
* Automatic R CMD check on Windows, macOS and Linux (GitHub Actions).
