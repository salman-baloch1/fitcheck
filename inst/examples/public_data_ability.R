# fitcheck on a public dataset: ability (ICAR sample items, psychTools package)
# ---------------------------------------------------------------------------
# It installs psychTools if needed, uses at most 2 cores, and writes everything
# to "public_data_ability_output.txt" in the current working directory.
#
# The dataset is public (it ships with the psychTools package; older versions
# of psych also included it), so no confidential
# examinee data are used here.
# ---------------------------------------------------------------------------

if (!requireNamespace("psychTools", quietly = TRUE)) install.packages("psychTools")
library(fitcheck)

out_file <- "public_data_ability_output.txt"
sink(out_file, split = TRUE)

cat("fitcheck on the ability data (psychTools)\n")
cat("Run on:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")

data("ability", package = "psychTools")
A <- as.data.frame(ability)

# Item types are identified from the column-name prefixes.
types <- c(Verbal = "^reason", Letter = "^letter", Matrix = "^matrix", Rotation = "^rotate")
sections <- lapply(types, function(p) grep(p, names(A), value = TRUE))
items <- unlist(sections, use.names = FALSE)

cat("Records:", nrow(A), "| items:", length(items), "\n")
print(sapply(sections, length))
cat("Item values:", paste(sort(unique(unlist(A[items]))), collapse = ", "), "\n")

# fitcheck needs complete 0/1 data, so records with any missing response are
# set aside and counted.
complete <- stats::complete.cases(A[items])
cat("Records with at least one missing response:", sum(!complete),
    sprintf("(%.1f%%)", 100 * mean(!complete)), "\n")
cat("Complete records used:", sum(complete), "\n\n")
X <- A[complete, items]

# Each item is its own indicator: four items per item type, one second-order
# factor. Items are treated as continuous indicators (MLR), as in the
# parcel-level analyses; this is a demonstration, not a full item-level model.
parcels <- stats::setNames(as.list(items), items)

report <- fitcheck(X, sections, parcels,
                   boot_ns = c(100, 300, 1000), reps = 200,
                   cores = 2, seed = 2027, B = 1000)
print(report)

cat("\nNote: the ability data have no test forms, so the by-form and invariance\n")
cat("checks are not run. The model has four first-order factors.\n\n")

print(citation("psychTools"))
cat("\n")
print(sessionInfo())
sink()
cat("\nDone. Output written to", normalizePath(out_file), "\n")
