# fitcheck

[![R-CMD-check](https://github.com/salman-baloch1/fitcheck/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/salman-baloch1/fitcheck/actions/workflows/R-CMD-check.yaml)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](LICENSE.md)
[![DOI](https://zenodo.org/badge/1470291412.svg)](https://doi.org/10.5281/zenodo.24539810)

**Check the data, then the solution, then the fit, for second-order factor models in operational testing.**

Second-order confirmatory factor models are routinely fitted to operational test data and judged by
fit indices against fixed cutoffs. fitcheck fixes the order of judgment:

| Step | Question | Functions |
|---|---|---|
| Check 0 | Do the records contain real responses? | `screen_records()`, `zero_probability()`, `screening_sensitivity()` |
| Forms | Does screening change the evidence within each form, and the comparison of forms? | `screening_by_form()`, `screening_invariance()` |
| Check 1 | Is the solution admissible, under both identifications? | `four_checks()`, `identification_check()` |
| Checks 2–4 | What do the indices say in context, with the unbiased SRMR and its interval? | `four_checks()`, `hierarchical_indices()` |
| Sample size | How would estimation and fit behave in smaller samples of *my* data? | `boot_by_n()`, `mini_mc()` |
| Reporting | Do section scores add value over the total score? Can two section scores be told apart? | `added_value()`, `profile_check()` |

`fitcheck()` runs all steps and returns one report. The package is built on
[lavaan](https://lavaan.ugent.be).

## Why this order

The evidence behind each step comes from a doctoral study of an operational university English test
(eight forms, 49,166 records; Rashid, 2026). Every number below can be reproduced from the analysis
script of that study.

* **All-zero records change the measurement evidence but not the fit.** In one form (3,881 records),
  14.8% of records had every item scored 0. Removing them hardly changed robust CFI (.998 vs .997)
  or RMSEA (.022 vs .021) of the second-order model, but the explained common variance (ECV) of a
  separately fitted bifactor model fell from .92 to .86, and Reading's average variance extracted in
  the correlated-factor model fell from .69 to .48, below the usual .50 criterion. Each pair was
  obtained with the same model and estimator (MLR).
* **Screening must be checked across forms too.** Across the eight forms, metric invariance was not
  supported with all records (ΔCFI = −.031, ΔRMSEA = +.043) or without the 3,060 all-zero records
  (ΔCFI = −.039, ΔRMSEA = +.039), so the invariance conclusion was unchanged after screening. But
  all-zero records were concentrated in some forms (0.7% to 14.8% per form), so this has to be
  shown, not assumed.
* **Second-order solutions often fail in small samples.** In 500 bootstrap samples of 100 records
  from one form, only 59.6% yielded admissible second-order solutions, and every failure was
  nonconvergence, with second-order loadings close to 1. In a simulated example in the
  package tests, an improper solution (negative disturbance variance) still had good
  conventional fit (robust CFI .959, SRMR .050 in the handbook run; the exact values vary slightly
  between computers). Fit computed from an improper solution should
  not be read.
* **Small samples displace fit indices.** Mean SRMR was .029 at n = 100 against .008 at n = 20,000
  (the large-sample reference in that study), a relative increase of 258%. The unbiased SRMR (Maydeu-Olivares, 2017)
  removes this upward shift, but it was truncated at zero in 68% of admissible solutions, so a value
  of 0 is a truncation, not proof of exact fit.
* **Screening can change reporting decisions.** Across eight forms (49,166 records; 45,176 after
  removing 3,060 all-zero (6.2%) and 930 whole-section-zero (1.9%) records), no section score showed statistically
  supported added value over the total score; in one form, keeping records with a whole section
  scored 0 made two section scores appear to add value.
* **The package reproduces the study.** Applied to the operational archive, fitcheck reproduced the
  screening counts (3,060 and 930 records; the two groups do not overlap) and the invariance results,
  with fit indices matching to three decimal places; admissibility at n = 100 was 65.4% against 59.6% in the study, with
  different random seeds, and every failure was again nonconvergence.
* **The same problems appear in public data.** On the ICAR reasoning items in the `ability` data of
  the psychTools package (16 items, four item types, 1,248 complete records), the full-sample
  solution was admissible (robust CFI .983), but only 67.5% of 200 bootstrap samples of 100 records
  gave admissible solutions, all failures being improper solutions, and the unbiased SRMR was zero in
  70% of these. One section (three-dimensional rotation) did add value over the total score, so
  the score-reporting check does not always favour the total score. The script is
  `inst/examples/public_data_ability.R`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("salman-baloch1/fitcheck", build_vignettes = TRUE)

# or from a local copy of the source package
install.packages("fitcheck_0.2.0.tar.gz", repos = NULL, type = "source")
```

Requires R (>= 4.1) and lavaan (>= 0.6-17). A step-by-step introduction is in the vignette:
`vignette("fitcheck", package = "fitcheck")`.

## Example (simulated data shipped with the package)

```r
library(fitcheck)
data(sim_test); data(sim_test_design)
items <- sim_test[, sim_test_design$items]

report <- fitcheck(items, sim_test_design$sections, sim_test_design$parcels,
                   form = sim_test$form, B = 1000)
report
```

Selected output:

```
==== CHECK 0  Data ======================
Records: 2000 | all-zero: 117 | whole-section zero: 24 | retained: 1859
 form records all_zero section_zero flagged_pct
    A     671        8            2         1.5
    B     620       97           21        19.0
    C     709       12            1         1.8

==== Screening sensitivity ==============
                    index all_records retained difference
                      CFI       0.999    1.000      0.001
                    RMSEA       0.015    0.007     -0.008
 mean_first_order_loading       0.834    0.748     -0.086
            omega_Grammar       0.830    0.738     -0.093
                  omega_h       0.924    0.885     -0.039

==== Screening sensitivity by form ======
 form                    index all_records retained difference
    A mean_first_order_loading       0.766    0.741     -0.025
    B mean_first_order_loading       0.898    0.757     -0.141
    C mean_first_order_loading       0.780    0.749     -0.031

==== Invariance across forms ============
Rule (Chen, 2007): the metric or scalar step is not supported when the decrease in robust CFI
  exceeds 0.01 AND the increase in robust RMSEA exceeds 0.015
     version      level    N admissible   CFI RMSEA   dCFI dRMSEA  decision
 all records configural 2000       TRUE 1.000 0.000     NA     NA  baseline
 all records     metric 2000       TRUE 0.999 0.015 -0.001  0.015 supported
 all records     scalar 2000       TRUE 0.999 0.013  0.000 -0.003 supported
    retained configural 1859       TRUE 1.000 0.000     NA     NA  baseline
    retained     metric 1859       TRUE 1.000 0.000  0.000  0.000 supported
    retained     scalar 1859       TRUE 1.000 0.000  0.000  0.000 supported
Highest level supported: all records = scalar | retained = scalar
Screening does not change the invariance conclusion.

==== Score reporting ====================
   section prmse_s prmse_x   VAR   diff diff_low diff_high       decision
 Listening   0.834   0.883 0.945 -0.049   -0.061    -0.035 no added value
   Grammar   0.721   0.851 0.848 -0.130   -0.158    -0.104 no added value
   Reading   0.839   0.867 0.968 -0.027   -0.040    -0.014 no added value
```

The flagged records sit mostly in form B, so screening changes the loadings of form B far more than those
of forms A and C. The invariance check then shows whether this changes the comparison of forms. In these
simulated data it does not; in operational data it can, so it should be reported. The decision is a
threshold-based assessment (Chen, 2007), not proof of invariance, and a second-order model is tested
through its correlated first-order version, which is equivalent only with three first-order factors.
The vignette and the demo script include a second simulated example in which 100 all-zero records
added to one form make metric invariance fail (ΔCFI = −.049, ΔRMSEA = +.030), while the retained
records support scalar invariance. Note that screening does not always improve fit: in form B, SRMR rises from .008 to .017. Screening
changes the evidence; it is not a way to obtain better fit.

How would the model behave with fewer candidates? Every replication has its own pre-drawn seed, so the
results are identical on one or several cores:

```r
d <- parcel_scores(items[report$screen$flags$keep, ], sim_test_design$parcels)
boot_by_n(d, sim_test_design$model, ns = c(100, 300, 1000), reps = 40, cores = 2)
#     N admissible_pct ci_low ci_high improper nonconverged   srmr usrmr_zero_pct
#   100           72.5   57.2    83.9       11            0 0.0378           72.4
#   300           80.0   65.2    89.5        8            0 0.0237           31.3
#  1000          100.0   91.2   100.0        0            0 0.0144            7.5
```

## Design choices

* **Admissibility before fit.** Fit indices from improper solutions are flagged, and every summary of
  repeated fits reports admissible, improper, and nonconverged solutions separately, with a Wilson interval
  for the admissible rate.
* **Two identifications.** With latent variances fixed to 1, a negative disturbance variance often shows up as
  nonconvergence; with marker variables it shows up as a Heywood case. `identification_check()` shows both.
* **Intervals, not cutoffs.** Added value is decided from a bootstrap interval for PRMSE_s − PRMSE_x
  (adds value / no added value / inconclusive), with all quantities recomputed in every resample.
* **Reproducibility.** Pre-drawn seeds per replication, no change to the user's random-number stream,
  identical results in serial and parallel runs.
* **Flag, do not delete.** Screening returns flags and a sensitivity analysis; whether flagged records are
  excluded remains the analyst's decision.

## Limits

* **Parcels.** The workflow fits the model to parcel scores with robust maximum likelihood, because
  operational forms often have 80 or more items. Parcelling can hide item-level misfit; item-level
  estimation of ordinal data (for example WLSMV) is a planned extension.
* **Screening is a judgment.** An all-zero record is flagged, not proven to be empty.
  `zero_probability()` shows how plausible a real all-zero score is, and `screening_sensitivity()`
  shows what the decision changes.
* **Model family.** The checks are written for second-order models. `four_checks()` works for any
  lavaan fit; `hierarchical_indices()` assumes one second-order factor.

## Related tools

fitcheck builds on known methods and does not replace the packages that implement them.
[lavaan](https://lavaan.ugent.be) provides estimation, the admissibility check, and `lavResiduals()`
for the unbiased SRMR. Omegas and related indices are also available in semTools and psych,
Haberman's added-value statistics in subscore, and simulation-based fit cutoffs in FCO. What
fitcheck adds is the order: records are screened, the effect of screening is shown, admissibility
is established under two identifications, and only then is fit read, all in one report.

## Citation

Rashid, S., Retnawati, H., & Hadi, S. (2026). *fitcheck: Check the data, then the solution, then the fit,
for second-order factor models* (Version 0.2.0) [Computer software]. Zenodo.
https://doi.org/10.5281/zenodo.24539812

To cite all versions of fitcheck, use https://doi.org/10.5281/zenodo.24539810.

In R: `citation("fitcheck")`.

## References

Chen, F. F. (2007). Sensitivity of goodness of fit indexes to lack of measurement invariance.
*Structural Equation Modeling, 14*(3), 464–504. https://doi.org/10.1080/10705510701301834

Haberman, S. J. (2008). When can subscores have value? *Journal of Educational and Behavioral Statistics,
33*(2), 204–229. https://doi.org/10.3102/1076998607302636

Maydeu-Olivares, A. (2017). Assessing the size of model misfit in structural equation models.
*Psychometrika, 82*(3), 533–558. https://doi.org/10.1007/s11336-016-9552-7

Rashid, S. (2026). *The impact of sample size and factor loadings on goodness-of-fit indices in
second-order confirmatory factor analysis: A bootstrapping study of the ProTEFL test data* (Doctoral
dissertation). Universitas Negeri Yogyakarta.

## License

GPL-3. See [LICENSE.md](LICENSE.md).
