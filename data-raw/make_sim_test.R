# Creates the simulated example data `sim_test` and `sim_test_design`.
# Simulated only: no real examinee records are used.
#
# Design: 2,000 records, three forms (A, B, C), 80 four-option items in three
# sections (Listening 30, Grammar 20, Reading 30) and ten parcels, a structure
# modelled on a three-section university English test. Responses follow a
# three-parameter logistic model (guessing .25) with a second-order factor
# structure. Like operational files, the data contain records with no real
# responses: all-zero records (concentrated in form B) and records with every
# Listening item scored 0.
set.seed(20270204)
N <- 2000
form <- sample(c("A", "B", "C"), N, replace = TRUE, prob = c(.35, .30, .35))
g <- rnorm(N)
gam <- c(Listening = .97, Grammar = .96, Reading = .94)
f <- sapply(gam, function(gk) gk * g + sqrt(1 - gk^2) * rnorm(N))
sec_of <- rep(1:3, c(30, 20, 30))
a <- runif(80, .5, 1.2); b <- rnorm(80, 0, .9); cc <- .25
X <- sapply(1:80, function(j) {
  p <- cc + (1 - cc) * plogis(1.7 * a[j] * (f[, sec_of[j]] - b[j]))
  as.integer(runif(N) < p)
})
# records without real responses
p_all <- ifelse(form == "B", .14, .02)
all_zero <- runif(N) < p_all
X[all_zero, ] <- 0L
lis_zero <- !all_zero & runif(N) < ifelse(form == "B", .05, .005)
X[lis_zero, 1:30] <- 0L
colnames(X) <- paste0("i", 1:80)
sim_test <- data.frame(id = sprintf("R%04d", 1:N), form = form, X, stringsAsFactors = FALSE)

items <- paste0("i", 1:80)
sections <- list(Listening = items[1:30], Grammar = items[31:50], Reading = items[51:80])
parcels <- list(L1 = items[1:10], L2 = items[11:20], L3 = items[21:30],
                G1 = items[31:40], G2 = items[41:50],
                R1 = items[51:56], R2 = items[57:62], R3 = items[63:68],
                R4 = items[69:74], R5 = items[75:80])
model <- paste("Listening =~ L1 + L2 + L3",
               "Grammar =~ G1 + G2",
               "Reading =~ R1 + R2 + R3 + R4 + R5",
               "G =~ NA*Listening + Grammar + Reading",
               "G ~~ 1*G", sep = "\n")
sim_test_design <- list(items = items, sections = sections, parcels = parcels, model = model)
save(sim_test, file = "data/sim_test.rda", compress = "xz")
save(sim_test_design, file = "data/sim_test_design.rda", compress = "xz")
