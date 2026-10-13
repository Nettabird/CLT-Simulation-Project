# helpers.R: shared settings and functions for the CLT Monte Carlo project.
# Every analysis file starts with source("R/helpers.R").

library(tidyverse)   # dplyr, purrr, ggplot2, tibble, ...
library(knitr)       # kable()
source("R/text.R")   # repeated text

# ---- Settings ---------------------------------------------------------------
SEED <- 3337
REPS <- as.integer(Sys.getenv("CLT_REPS", 50000))   # replications: simulated sample means per n
NS   <- c(1, 2, 5, 10, 20, 30, 50, 100, 200, 400, 700, 1000)   # default 12 sample sizes

# "Normal enough" criteria. Defended in summary file. Monte Carlo noise with 50,000 reps:
# mean SE ~ 0.004, standard deviation SE ~ 0.003, skewness SE ~ 0.011,
# tail proportion SE ~ 0.0007, and KS distance for an exactly normal sample averages 0.004.
# **Note:** Changing criteria may require changing the values of n in the list of sample sizes.
CRITERIA <- list(
  mean_z   = 0.015,              # |mean of z| below this
  sd_ratio = c(0.99, 1.01),       # observed SD / theoretical SE
  skew       = 0.05,             # |skewness| below this
  ex_kurt    = 0.05,             # |excess kurtosis| below this
  ks         = 0.03,             # KS distance to N(0, 1) below this
  tail_total = c(0.04, 0.06),    # P(|Z| > 1.96), nominal 0.05
  tail_each  = c(0.015, 0.035)  # each tail, nominal 0.025
)

# Show the current criteria as a table.
criteria_table <- function(crit = CRITERIA) {
  fmt_range <- \(x) paste(x[1], "to", x[2])
  tibble(
    Diagnostic  = c("Absolute mean of z", "Observed SD / theoretical SE",
                    "Absolute skewness", "Absolute excess kurtosis", "KS distance",
                    "Total tail share", "Each tail share"),
    Requirement = c(paste("below", crit$mean_z), fmt_range(crit$sd_ratio),
                    paste("below", crit$skew), paste("below", crit$ex_kurt),
                    paste("below", crit$ks), fmt_range(crit$tail_total),
                    fmt_range(crit$tail_each))
  ) |> kable()
}

# Where the mean of iid draws crosses the threshold, given one draw's
# skewness (gamma) or excess kurtosis (kappa): gamma / sqrt(n) and kappa / n.
n_skew <- \(gamma, crit = CRITERIA) (abs(gamma) / crit$skew)^2
n_kurt <- \(kappa, crit = CRITERIA) abs(kappa) / crit$ex_kurt

# Smallest whole n at which the theoretical skewness and excess kurtosis of the
# mean of iid draws are strictly inside the thresholds, rounded up to a coarser
# step for larger n: 5 below 100, 10 from 100 to 499, 20 from 500 up.
n_theory <- function(gamma, kappa, crit = CRITERIA) {
  need  <- max(n_skew(gamma, crit), n_kurt(kappa, crit))
  n_min <- floor(round(need, 8)) + 1       # strict inequality, as in apply_criteria()
  step  <- if (n_min < 100) 5 else if (n_min < 500) 10 else 20
  ceiling(n_min / step) * step
}

# Smallest whole n strictly beyond a crossover value (apply_criteria() uses strict inequality).
n_pass <- \(x) floor(round(x, 8)) + 1

# First grid value at or above n (NA if n is beyond the grid).
first_grid <- \(n, ns = NS) ns[ns >= n][1]

# Rough KS distance to N(0, 1) for the mean of n iid draws.
# sigma only matters for discrete populations (lattice term).
ks_pred <- function(n, gamma, kappa, sigma = NA, floor = 0.004) {
  edge    <- 0.0665 * abs(gamma) / sqrt(n) + 0.023 * abs(kappa) / n
  lattice <- if (is.na(sigma)) 0 else 0.2 / (sigma * sqrt(n))
  floor + edge + lattice
}

# Smallest whole n at which the predicted KS distance is below the threshold.
n_ks <- function(gamma, kappa, sigma = NA, crit = CRITERIA) {
  n <- 1
  while (ks_pred(n, gamma, kappa, sigma) >= crit$ks) n <- n + 1
  n
}

n_ks_lattice <- \(sigma, crit = CRITERIA, floor = 0.004) n_pass((0.2 / (sigma * (crit$ks - floor)))^2)

# ---- Populations ------------------------------------------------------------
# A population is a list with:
#   mean_fun(n, reps): returns `reps` simulated sample means of size n
#   mu(n), se(n):      true mean and standard error of the sample mean (NA if undefined)
#   lattice(n):        spacing of possible standardized means (discrete populations), else NULL

rdsix <- function(n) sample.int(6, n, replace = TRUE)   # fair six-sided die

pop_iid <- function(name, rng, mu = NA_real_, sigma = NA_real_, discrete = FALSE) {
  list(
    name     = name,
    mean_fun = function(n, reps) replicate(reps, mean(rng(n))),
    mu       = function(n) mu,
    se       = function(n) sigma / sqrt(n),
    # means of integer-valued data sit on a grid with spacing 1/n; divide by the SE
    lattice  = if (discrete) function(n) 1 / (sigma * sqrt(n)) else NULL
  )
}

# ---- Simulation -------------------------------------------------------------
# set.seed(SEED) is called here, so every population uses the same seed.
simulate_means <- function(pop, ns = NS, reps = REPS) {
  set.seed(SEED)
  map(ns, \(n) tibble(n = n, xbar = pop$mean_fun(n, reps))) |> list_rbind()
}

# Standardize each sample mean: z = (xbar - mu) / SE.
# If the SE is undefined (infinite variance), use a robust standardization
# (median and IQR / 1.349), which compares the SHAPE of the distribution to normal.
add_z <- function(sim, pop, se = pop$se) {
  sim |>
    group_by(n) |>
    mutate(
      se_n = se(n[1]),
      mu_n = pop$mu(n[1]),
      z = if (is.na(se_n[1])) (xbar - median(xbar)) / (IQR(xbar) / 1.349)
          else (xbar - mu_n[1]) / se_n[1]
    ) |>
    ungroup()
}

# ---- Diagnostics ------------------------------------------------------------
skew    <- function(x) mean((x - mean(x))^3) / mean((x - mean(x))^2)^1.5
ex_kurt <- function(x) mean((x - mean(x))^4) / mean((x - mean(x))^2)^2 - 3
ks_dist <- function(z) unname(suppressWarnings(ks.test(z, "pnorm")$statistic))

# Shapiro-Wilk on many small batches: a normal should be rejected about 5% of the time.
sw_reject_rate <- function(z, size = 50) {
  z |> split(ceiling(seq_along(z) / size)) |>
    map_lgl(\(b) tryCatch(shapiro.test(b)$p.value < 0.05, error = \(e) NA)) |>
    mean(na.rm = TRUE)
}

diagnose <- function(sim) {
  sim |>
    group_by(n) |>
    summarise(
      mean_xbar = mean(xbar),
      mu_n      = mu_n[1],
      mean_z    = if (is.na(se_n[1])) NA_real_ else mean(z),   # bias in SE units
      sd_xbar   = sd(xbar),
      se_n      = se_n[1],
      skew      = skew(z),
      ex_kurt   = ex_kurt(z),
      sw_p      = shapiro.test(head(z, 5000))$p.value,   # one test on up to 5,000 means
      sw_reject = sw_reject_rate(z),                     # rejection rate over batches of 50
      ks        = ks_dist(z),
      tail_lo   = mean(z < -1.96),
      tail_hi   = mean(z > 1.96),
      sd_ratio  = sd(xbar) / se_n[1],   # observed SD of the means / theoretical SE (want 1)
      iqr_xbar  = IQR(xbar),            # spread of the raw means (constant if the mean doesn't concentrate)
      .groups = "drop"
    ) |>
    mutate(tail = tail_lo + tail_hi)
}

# apply_criteria <- function(d, crit = CRITERIA) {
#   d |> mutate(pass = abs(skew) < crit$skew & abs(ex_kurt) < crit$ex_kurt & ks < crit$ks &
#                 (is.na(mean_z)   | abs(mean_z) < crit$mean_z) &
#                 (is.na(sd_ratio) | between(sd_ratio, crit$sd_ratio[1], crit$sd_ratio[2])) &
#                 between(tail, crit$tail_total[1], crit$tail_total[2]) &
#                 between(tail_lo, crit$tail_each[1], crit$tail_each[2]) &
#                 between(tail_hi, crit$tail_each[1], crit$tail_each[2]))
# }

# Pass/fail for each criterion separately (same row order as d).
criterion_pass <- function(d, crit = CRITERIA) {
  d |>
    mutate(
      .keep = "none",
      n,
      `mean of z`   = is.na(mean_z) | abs(mean_z) < crit$mean_z,
      `SD ratio`    = is.na(sd_ratio) | between(sd_ratio, crit$sd_ratio[1], crit$sd_ratio[2]),
      skewness      = abs(skew) < crit$skew,
      kurtosis      = abs(ex_kurt) < crit$ex_kurt,
      `KS distance` = ks < crit$ks,
      `total tail`  = between(tail, crit$tail_total[1], crit$tail_total[2]),
      `lower tail`  = between(tail_lo, crit$tail_each[1], crit$tail_each[2]),
      `upper tail`  = between(tail_hi, crit$tail_each[1], crit$tail_each[2])
    )
}

apply_criteria <- function(d, crit = CRITERIA) {
  cp <- criterion_pass(d, crit)
  m  <- as.matrix(select(cp, -n))
  d |> left_join(
    tibble(n       = cp$n,
           pass    = rowSums(m) == ncol(m),
           failing = apply(m, 1, \(r) paste(colnames(m)[!r], collapse = ", "))),
    by = "n"
  )
}

# The last criterion (or criteria) to meet its threshold and stay there.
last_to_pass <- function(d, crit = CRITERIA) {
  ok <- criterion_pass(d, crit) |> arrange(n)
  settle <- map_dbl(select(ok, -n), \(p) {
    stay <- rev(cumall(rev(p)))
    if (any(stay)) min(ok$n[stay]) else Inf
  })
  if (any(is.infinite(settle))) {
    return(paste0("Not met by n = ", max(ok$n), ": ",
                  paste(names(settle)[is.infinite(settle)], collapse = ", ")))
  }
  last <- max(settle)
  if (last == min(ok$n)) return(paste0("None (all met at n = ", last, ")"))
  paste0(paste(names(settle)[settle == last], collapse = ", "), " (n = ", last, ")")
}

# Smallest n where the criteria hold at that n and at every larger n (NA if none).
suggest_n <- function(d, crit = CRITERIA) {
  d  <- apply_criteria(d, crit) |> arrange(n)
  ok <- rev(cumall(rev(d$pass)))
  if (any(ok)) min(d$n[ok]) else NA_real_
}

# ---- Tables -----------------------------------------------------------------

COLUMN_LABELS <- c(
  n                 = "Sample size",
  pass              = "Meets criteria",
  mean_xbar         = "Mean of sample means",
  mean_z            = paste0("Mean of standardized means (below ", CRITERIA$mean_z, ")"),
  mu_n              = "Theoretical mean",
  sd_xbar           = "SD of sample means",
  sd_ratio          = paste0("Observed SD / theoretical SE (", CRITERIA$sd_ratio[1], " to ", CRITERIA$sd_ratio[2], ")"),
  skew_theory       = "Theoretical skewness",
  skew              = paste0("Skewness (below ", CRITERIA$skew, ")"),
  ex_kurt           = paste0("Excess kurtosis (below ", CRITERIA$ex_kurt, ")"),
  sw_p              = "Shapiro-Wilk p-value",
  sw_reject         = "Shapiro-Wilk rejection rate",
  ks                = paste0("KS distance (below ", CRITERIA$ks, ")"),
  tail_lo           = paste0("Lower tail share (", CRITERIA$tail_each[1], " to ", CRITERIA$tail_each[2], ")"),
  tail_hi           = paste0("Upper tail share (", CRITERIA$tail_each[1], " to ", CRITERIA$tail_each[2], ")"),
  tail              = paste0("Total tail share (", CRITERIA$tail_total[1], " to ", CRITERIA$tail_total[2], ")"),
  iqr_xbar          = "IQR of sample means",
  se_n              = "Theoretical SE",
  standardized_with = "Standard error used",
  failing           = "Criteria not met"
)

show_table <- function(d, digits = 3) {
  kable(d, digits = digits, col.names = coalesce(unname(COLUMN_LABELS[names(d)]), names(d)))
}

# ---- Plots ------------------------------------------------------------------
# Histogram bins: Freedman-Diaconis width; for discrete data, a whole number of
# lattice steps per bin, centered on lattice points (avoids empty/doubled bins).
# (Values beyond +/- lim are dropped for display only, so heavy-tailed densities look slightly inflated.)
hist_data <- function(z, step = NULL, lim = 4) {
  z  <- z[abs(z) <= lim]
  bw <- 2 * IQR(z) / length(z)^(1/3)
  if (!is.null(step)) bw <- step * max(1, round(bw / step))
  k  <- seq(floor((min(z) - z[1]) / bw), ceiling((max(z) - z[1]) / bw) + 1)
  h  <- hist(z, breaks = z[1] + (k - 0.5) * bw, plot = FALSE)
  tibble(mid = h$mids, density = h$density, width = bw)
}

plot_hist <- function(sim, pop, ns_show = sort(unique(sim$n))) {
  sim |>
    filter(n %in% ns_show) |>
    group_by(n) |>
    reframe(hist_data(z, step = if (is.null(pop$lattice)) NULL else pop$lattice(n[1]))) |>
    ggplot(aes(xmin = mid - width / 2, xmax = mid + width / 2, ymin = 0, ymax = density)) +
    geom_rect(fill = "steelblue", alpha = 0.7) +
    stat_function(fun = dnorm, xlim = c(-4, 4), n = 201, inherit.aes = FALSE, color = "red") +
    facet_wrap(~n, ncol = 4, labeller = label_both) +
    coord_cartesian(xlim = c(-4, 4)) +
    labs(x = "Standardized sample mean", y = "Density",
         subtitle = "Red curve: standard normal density")
}

# Set ylim = NULL to see the full range (e.g., Cauchy).
plot_qq <- function(sim, ns_show = sort(unique(sim$n)), ylim = c(-5, 5)) {
  sim |>
    filter(n %in% ns_show) |>
    ggplot(aes(sample = z)) +
    stat_qq(size = 0.3, alpha = 0.4) +
    geom_abline(color = "red") +
    facet_wrap(~n, ncol = 4, labeller = label_both) +
    coord_cartesian(ylim = ylim) +
    labs(x = "Theoretical normal quantile", y = "Sample quantile")
}

plot_convergence <- function(d, crit = CRITERIA) {
  stats <- c("mean_z", "sd_ratio", "skew", "ex_kurt", "ks")
  
  # dashed red criteria lines
  thr <- tibble(
    name  = c("mean_z", "mean_z", "sd_ratio", "sd_ratio",
              "skew", "skew", "ex_kurt", "ex_kurt", "ks"),
    value = c(crit$mean_z, -crit$mean_z, crit$sd_ratio[1], crit$sd_ratio[2],
              crit$skew, -crit$skew, crit$ex_kurt, -crit$ex_kurt, crit$ks)
  )
  
  # solid gray reference lines at the theoretical value
  ref <- tibble(name = c("mean_z", "sd_ratio", "skew", "ex_kurt"),
                value = c(0, 1, 0, 0))

  long <- d |>
    select(n, all_of(stats)) |>
    pivot_longer(-n, values_drop_na = TRUE) |>
    mutate(name = factor(name, levels = stats))
  
  # drop panels with no data (mean_z and sd_ratio are NA when the SE is undefined)
  keep <- levels(long$name)[levels(long$name) %in% long$name]
  thr  <- thr |> filter(name %in% keep) |> mutate(name = factor(name, levels = stats))
  ref  <- ref |> filter(name %in% keep) |> mutate(name = factor(name, levels = stats))
  
  ggplot(long, aes(n, value)) +
    geom_hline(data = ref, aes(yintercept = value), color = "gray50") +
    geom_hline(data = thr, aes(yintercept = value), linetype = "dashed", color = "red") +
    geom_line() + geom_point() +
    facet_wrap(~name, scales = "free_y", ncol = 3) +
    scale_x_log10() +
    labs(x = "Sample size n (log scale)", y = NULL,
         subtitle = "Dashed red: criteria; solid gray: theoretical value")
}

plot_skew_kurt <- function(d, crit = CRITERIA, ylim = NULL, label = NULL) {
  thr <- tibble(statistic = rep(c("Skewness", "Excess kurtosis"), each = 2),
                value     = c(crit$skew, -crit$skew, crit$ex_kurt, -crit$ex_kurt))
  d |>
    select(n, Skewness = skew, `Excess kurtosis` = ex_kurt) |>
    pivot_longer(-n, names_to = "statistic") |>
    ggplot(aes(n, value, color = statistic)) +
    geom_hline(yintercept = 0, color = "gray50") +
    geom_hline(data = thr, aes(yintercept = value, color = statistic),
               linetype = "dashed", alpha = 0.6) +
    geom_line() + geom_point() +
    scale_x_log10() +
    coord_cartesian(ylim = ylim) +
    labs(title = label, x = "Sample size n (log scale)", y = "Value", color = NULL,
         subtitle = "Dashed lines: criteria; solid gray line: theoretical normal (0)")
}

# ---- Saving results for the summary ----------------------------------------
# selected_n: a value from NS, Inf if you conclude no reasonable n exists,
# or NA while you are undecided (the summary then shows the suggested n, marked).
save_results <- function(id, label, diag, selected_n = NA, observation = "TODO") {
  stopifnot(is.na(selected_n) || is.infinite(selected_n) || selected_n %in% diag$n)
  dir.create("results", showWarnings = FALSE)
  saveRDS(list(id = id, label = label, diag = diag, selected_n = selected_n,
               suggested_n = suggest_n(diag), observation = observation),
          file.path("results", paste0(id, ".rds")))
}
