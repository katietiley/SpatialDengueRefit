# =============================================================================
# Rio de Janeiro Wolbachia / CATI scenario analyses (SpatialDengueRefit)
# Run from the repo root after Rio_refit_workflow.Rmd has produced the inputs
# in data/ (radiation-model finalWeights, monthly calendar, seed fitdat).
# =============================================================================

#### Load local repo ####
devtools::load_all(".", reset = TRUE) # run from repo root

#### Required libraries ####
library(geosphere)    # required for distm()
library(mgcv)         # gam() for the calibration emulator
library(tidyverse)    # dplyr, tidyr, ggplot2, tibble, purrr
library(lubridate)
library(cowplot)
library(png)
library(grid)

#### Load inputs ####
load_data_object <- function(file) {
  e <- new.env()
  obj_name <- load(file, envir = e)
  e[[obj_name]]
}
timestep  <- load_data_object("data/timestep_rio.rda")   # c(1, 12): months; seed rows have Week == 1
fitdat     <- load_data_object("data/fitdat_rio.rda")      # day-1 seed infections (Nov 2011 x 7/30)
pastdat    <- load_data_object("data/pastdat_rio.rda")
riopop     <- load_data_object("data/riopop_rio.rda")
unipix     <- load_data_object("data/unipix_rio.rda")
pixdistmat <- load_data_object("data/pixdistmat_rio.rda")
stim_vec   <- load_data_object("data/stim_vec_rio.rda")
load("data/sim_calendar_rio.rda")                          # n_days, day_month, fit_days
data(wmel_rio)

data(finalWeights)                                         # Rio radiation posterior, Brady weights
stopifnot(nrow(finalWeights) == 1000, all(finalWeights$Mov_model_type == 3))
stopifnot(n_days == 366)
stopifnot(!is.unsorted(unipix$patchID))                   # model outputs are in unipix (patchID) order
stopifnot(max(wmel_rio$wmel_pctg, na.rm = TRUE) <= 1)      # Wolbachia prevalence must be a proportion for beta_map()

OUT <- "26-11-21_"   # prefix for all saved outputs from this run

# ------------------------------------------------------------
# Paired-run settings and helper functions
# ------------------------------------------------------------
N_RUNS    <- 100    # stochastic runs per scenario, all analyses
SEED_BASE <- 1000   # run i uses seed SEED_BASE + i in every scenario
patch_pop <- unipix$pop

# Run i uses seed SEED_BASE + i in every scenario, so it draws the same posterior
# parameter set, transmission landscape and initial infections; scenarios differ
# only through the intervention. Failed runs are dropped with a warning.

# Simulations run in parallel with parallel::mclapply() on 8 cores on Mac.
# mclapply() relies on forking, which is unavailable on Windows: there it runs on a
# single core (or set mc.cores = 1). Reduce the number of cores if memory is limited.

run_raw_scenario <- function(paramsList, nruns = N_RUNS, ncores = 8, seed_base = SEED_BASE) {
  out <- parallel::mclapply(seq_len(nruns), function(i) {
    set.seed(seed_base + i)
    res <- DEN.spatial(timestep = timestep, fitdat = fitdat, pastdat = pastdat,
                       sgpop = riopop, unipix = unipix, pixdistmat = pixdistmat,
                       steprun = n_days, paramsList = paramsList)
    if (!is.list(res)) return(NULL)
    newD_mat <- do.call(rbind, lapply(res$Patchtotals, function(x) x$newD))
    list(run = i,
         daily_cases = rowSums(newD_mat, na.rm = TRUE),
         patch_cases = colSums(newD_mat, na.rm = TRUE))
  }, mc.cores = ncores, mc.preschedule = FALSE)
  ok <- vapply(out, function(o) is.list(o) && !inherits(o, "try-error"), logical(1))
  if (any(!ok)) warning(sum(!ok), " of ", nruns, " runs failed and were dropped")
  out[ok]
}

# Per-run effectiveness against the paired baseline run
paired_stats <- function(raw_int, raw_base, patch_pop) {
  base_by_run <- setNames(raw_base, sapply(raw_base, `[[`, "run"))
  raw_int <- Filter(function(x) as.character(x$run) %in% names(base_by_run), raw_int)
  bind_rows(lapply(raw_int, function(x) {
    b  <- base_by_run[[as.character(x$run)]]$patch_cases
    pe <- ifelse(b > 0, 1 - x$patch_cases / b, NA_real_)
    pe <- pmin(pmax(pe, 0), 1)
    tibble(
      run = x$run,
      total_cases   = sum(x$patch_cases),
      cases_averted = sum(b) - sum(x$patch_cases),
      effectiveness_total        = 1 - sum(x$patch_cases) / sum(b),
      effectiveness_pop_weighted = weighted.mean(pe, patch_pop, na.rm = TRUE),
      effectiveness_unweighted   = mean(pe, na.rm = TRUE))
  }))
}

# Per-patch median of paired per-run effectiveness (for maps)
paired_patch_effectiveness <- function(raw_int, raw_base) {
  base_by_run <- setNames(raw_base, sapply(raw_base, `[[`, "run"))
  raw_int <- Filter(function(x) as.character(x$run) %in% names(base_by_run), raw_int)
  eff_mat <- sapply(raw_int, function(x) {
    b <- base_by_run[[as.character(x$run)]]$patch_cases
    pmin(pmax(ifelse(b > 0, 1 - x$patch_cases / b, NA_real_), 0), 1)
  })
  apply(eff_mat, 1, function(v) if (all(is.na(v))) NA_real_ else median(v, na.rm = TRUE))
}

# ------------------------------------------------------------
# Paired no-intervention baseline (shared by calibration and RQ1-RQ4)
# ------------------------------------------------------------
raw_default <- run_raw_scenario(list(stim = stim_vec))
stopifnot(length(raw_default) == N_RUNS)   # every run needs a baseline partner
saveRDS(raw_default, paste0(OUT, "raw_default.rds"))
summary(sapply(raw_default, function(x) sum(x$daily_cases)))



#### 0. SUPP Calibrate betaEnv_reduction for dengue case reduction #### ---------------------------------------------------------

# ------------------------------------------------------------
# Literature values
# ------------------------------------------------------------
lit_points <- data.frame(
  wolbachia_prevalence = c(0, 32, 40, 60, 55, 73, 80, 95, 95.8, 100),
  dengue_reduction = c(0, 38, 46, 69.4, 75.9, 73, 77, 88.8, 77.1, 76),
  dengue_reduction_lower = c(0, 32, 21, 54.4, 58, 49, 63, 77.9, 65.3, 60),
  dengue_reduction_upper = c(0, 44, 63, 79.4, 86, 86, 86, 94.4, 83.1, 86),
  source = c("Biological constraint", "Rio de Janeiro average", "Niteroi", "Niteroi", "Niteroi",
             "Yogyakarta non-RCT", "Niteroi", "Niteroi long-term", "Yogyakarta RCT", "Yogyakarta non-RCT")
)
lit_points$source <- factor(
  lit_points$source,
  levels = c(
    "Biological constraint",
    "Niteroi",
    "Niteroi long-term",
    "Rio de Janeiro average",
    "Yogyakarta RCT",
    "Yogyakarta non-RCT"
  )
)

### plot panel A - lit points
supp_A_plot_base <- ggplot(lit_points) +
  geom_errorbar(aes(x = wolbachia_prevalence / 100,
                    ymin = dengue_reduction_lower / 100,
                    ymax = dengue_reduction_upper / 100),
                width = 0.015, linewidth = 0.6, colour = "black") +
  geom_point(aes(x = wolbachia_prevalence / 100,
                 y = dengue_reduction / 100,
                 shape = source),
             size = 3, colour = "black") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_shape_manual(
    values = c(
      "Biological constraint" = 4,
      "Niteroi" = 16,
      "Niteroi long-term" = 1,
      "Rio de Janeiro average" = 17,
      "Yogyakarta RCT" = 15,
      "Yogyakarta non-RCT" = 0
    )
  ) +
  theme_bw()

supp_A_plot <- supp_A_plot_base +
  labs(x = "",
       y = "Dengue reduction",
       shape = "Literature source") +
  theme(
    legend.position = c(0.98, 0.02),
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.key.size = unit(0.6, "lines"),
    legend.key.height = unit(0.9, "lines")
  )

supp_A_plot

# ------------------------------------------------------------
# Fit linear model, line of best fit through points
# ------------------------------------------------------------
lit_fit <- lm(dengue_reduction ~ 0 + wolbachia_prevalence, data = lit_points) # forcing line of best to go through origin

# extrapolate through all wolbachia_prevalence values
fit_line <- data.frame(wolbachia_prevalence = seq(0, 100, by = 1))
fit_line$dengue_reduction <- predict(lit_fit, newdata = fit_line)
fit_line$dengue_reduction <- pmax(fit_line$dengue_reduction, 0) # force prediction line not to go below 0

### plot panel B - lit line of best fit
supp_B_plot <- supp_A_plot_base +
  geom_line(
    data = fit_line,
    aes(x = wolbachia_prevalence / 100,
        y = dengue_reduction / 100),
    colour = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  labs(x = "",
       y = "",
       shape = "Literature source") +
  theme(legend.position = "none")

supp_B_plot


# ------------------------------------------------------------
# Calibrate betaEnv_reduction against literature dengue reductions (paired runs)
# ------------------------------------------------------------
# target dengue reductions from literature linear fit
target_df <- fit_line %>%
  filter(wolbachia_prevalence %% 5 == 0,
         wolbachia_prevalence >= 5) %>%
  transmute(wolbachia_prevalence = wolbachia_prevalence / 100,
            target_reduction = dengue_reduction)

# betaEnv_reduction values to test
beta_reduction_grid <- seq(0.02, 0.98, by = 0.02)

# each candidate run i is paired with baseline run i (same posterior draw, landscape,
# initial conditions and seed); reduction = 1 - total cases / paired baseline total
run_reduction_candidate <- function(betaEnv_reduction_value) {
  message("Calibrating betaEnv_reduction = ", betaEnv_reduction_value)
  raw <- run_raw_scenario(
    list(stim = stim_vec, betaEnv_reduction = betaEnv_reduction_value),
    nruns = N_RUNS, seed_base = SEED_BASE
  )
  red <- paired_stats(raw, raw_default, patch_pop)$effectiveness_total
  tibble(
    betaEnv_reduction    = betaEnv_reduction_value,
    model_reduction      = 100 * median(red, na.rm = TRUE),
    model_reduction_q025 = 100 * quantile(red, 0.025, na.rm = TRUE),
    model_reduction_q975 = 100 * quantile(red, 0.975, na.rm = TRUE)
  )
}

beta_calib_results <- purrr::map_dfr(beta_reduction_grid, run_reduction_candidate)

# check: modelled reduction should rise smoothly with betaEnv_reduction
plot(beta_calib_results$betaEnv_reduction, beta_calib_results$model_reduction,
     type = "b", xlab = "betaEnv reduction", ylab = "Modelled dengue reduction (%)")

# for each target Wolbachia prevalence, select the candidate reduction whose modelled
# dengue reduction is closest to the literature target (nearest neighbour, absolute error)
beta_calibration_lookup <- target_df %>%
  tidyr::crossing(beta_calib_results) %>%
  dplyr::mutate(
    abs_error = abs(model_reduction - target_reduction)
  ) %>%
  dplyr::group_by(wolbachia_prevalence) %>%
  dplyr::slice_min(abs_error, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    wolbachia_prevalence,
    target_reduction,
    betaEnv_reduction,
    model_reduction,
    abs_error
  )

beta_calibration_lookup

### plot panel C - wolbachia_prevalence by betaEnv_reduction value points
supp_C_plot_base <- ggplot(beta_calibration_lookup, aes(x = wolbachia_prevalence, y = betaEnv_reduction)) +
  geom_line() +
  geom_point(size = 2.5) +
  scale_x_continuous(labels = scales::percent_format()) +
  theme_bw()

supp_C_plot <- supp_C_plot_base +
  labs(x = "Wolbachia prevalence",
       y = expression("Calibrated " * beta[mu] * " reduction"))

supp_C_plot
# ------------------------------------------------------------
# Fit non-linear model line of best fit through points
# ------------------------------------------------------------
beta_fit <- gam(betaEnv_reduction ~ s(wolbachia_prevalence, k = 4), data = beta_calibration_lookup)

# extrapolate through all wolbachia_prevalence values
beta_fit_line <- data.frame(wolbachia_prevalence = seq(0.05, 1, by = 0.01))
beta_fit_line$beta_pred <- predict(beta_fit, newdata = beta_fit_line)

### plot panel D - GAM through calibrated points
supp_D_plot <- supp_C_plot_base +
  geom_line(
    data = beta_fit_line,
    aes(y = beta_pred),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  labs(x = "Wolbachia prevalence",
       y = "")

supp_D_plot

### combined plot panel
supp_beta_panel <- cowplot::plot_grid(
  supp_A_plot,
  supp_B_plot,
  supp_C_plot,
  supp_D_plot,
  labels = c("A", "B", "C", "D"),
  label_fontface = "bold",
  label_size = 16,
  ncol = 2,
  align = "hv"
)

supp_beta_panel
ggsave(paste0(OUT, "FIGURE_supp_beta_panel.png"), supp_beta_panel, width = 12, height = 10, dpi = 300)

saveRDS(
  list(
    lit_fit = lit_fit,
    fit_line = fit_line,
    target_df = target_df,
    beta_calib_results = beta_calib_results,
    beta_calibration_lookup = beta_calibration_lookup,
    beta_fit = beta_fit,
    beta_fit_line = beta_fit_line
  ),
  file = paste0(OUT, "beta_calibration.rds")
)



# ------------------------------------------------------------
# Validate GAM emulator
# ------------------------------------------------------------

lookup <- beta_calibration_lookup
lookup$gam_prediction <- predict(beta_fit, newdata = lookup)

### compare observed vs predicted beta

lookup$error <- lookup$gam_prediction - lookup$betaEnv_reduction
lookup$abs_error <- abs(lookup$error)

summary(lookup$error)
cat("\nRMSE:", sqrt(mean(lookup$error^2)), "\n")
cat("Mean absolute error:", mean(lookup$abs_error), "\n")
cat("Maximum absolute error:", max(lookup$abs_error), "\n")
cat("R-squared:", cor(lookup$betaEnv_reduction, lookup$gam_prediction)^2, "\n")

### diagnostic plot

ggplot(lookup, aes(betaEnv_reduction, gam_prediction)) +
  geom_point(size=2) +
  geom_abline(slope=1, intercept=0,
              colour="red", linetype="dashed") +
  coord_equal() +
  labs(x="Calibration beta reduction",
       y="GAM prediction")

### residuals plot

ggplot(lookup, aes(wolbachia_prevalence, error)) +
  geom_hline(yintercept=0, linetype="dashed") +
  geom_point(size=2) +
  geom_line() +
  labs(x="Wolbachia prevalence",
       y="Prediction error")

### lookup vs GAM table

lookup %>%
  transmute(
    prevalence = scales::percent(wolbachia_prevalence),
    lookup_beta = round(betaEnv_reduction, 3),
    gam_beta = round(gam_prediction, 3),
    difference = round(error, 3)
  )


# ------------------------------------------------------------
# Calibration emulator: Wolbachia prevalence -> betaEnv reduction
# ------------------------------------------------------------
# to reuse a saved calibration instead of rerunning section 0:
# beta_fit <- readRDS(paste0(OUT, "beta_calibration.rds"))$beta_fit
beta_map <- function(w) {
  pred <- predict(beta_fit, newdata = data.frame(wolbachia_prevalence = w))
  pred[w == 0] <- 0      # enforce boundary condition explicitly
  pmin(pmax(pred, 0), 1) # ensure 0,1 boundedness
}


#### 1. RQ1 outbreaks under beta baselines #### ---------------------------------------------------------

# ------------------------------------------------------------
# 1. Reconstruct patch-level Wolbachia prevalence
# ------------------------------------------------------------
# mean of final 3 non-missing observations per deployment cell
wmel_patch_prev <- wmel_rio %>%
  filter(!is.na(wmel_pctg)) %>%
  arrange(cell.id, date) %>%
  group_by(cell.id) %>%
  slice_tail(n = 3) %>%
  summarise(wmel_prev = mean(wmel_pctg, na.rm = TRUE),
            .groups = "drop")

wmel_patch <- unipix %>%
  left_join(wmel_patch_prev,
            by = c("pixID" = "cell.id")) %>%
  mutate(wmel_prev = ifelse(is.na(wmel_prev), 0, wmel_prev))

# ------------------------------------------------------------
# 2. Map Wolbachia prevalence to beta_mu reduction
# ------------------------------------------------------------
wmel_patch <- wmel_patch %>%
  mutate(betaEnv_reduction = beta_map(wmel_prev))

betaEnv_reduction_vec <- wmel_patch %>%
  arrange(patchID) %>%
  pull(betaEnv_reduction)

stopifnot(length(betaEnv_reduction_vec) == nrow(unipix))

# population-weighted summary stats
wmel_effectiveness_summary <- wmel_patch %>%
  summarise(
    # population-weighted means
    pop_weighted_wmel_prev = weighted.mean(wmel_prev, pop, na.rm = TRUE),
    pop_weighted_betaEnv_reduction = weighted.mean(betaEnv_reduction, pop, na.rm = TRUE),

    # median
    median_wmel_prev = median(wmel_prev, na.rm = TRUE),
    median_betaEnv_reduction = median(betaEnv_reduction, na.rm = TRUE),

    # IQR
    wmel_prev_q25 = quantile(wmel_prev, 0.25, na.rm = TRUE),
    wmel_prev_q75 = quantile(wmel_prev, 0.75, na.rm = TRUE),
    beta_q25 = quantile(betaEnv_reduction, 0.25, na.rm = TRUE),
    beta_q75 = quantile(betaEnv_reduction, 0.75, na.rm = TRUE),

    # range
    min_wmel_prev = min(wmel_prev, na.rm = TRUE),
    max_wmel_prev = max(wmel_prev, na.rm = TRUE),
    min_betaEnv_reduction = min(betaEnv_reduction, na.rm = TRUE),
    max_betaEnv_reduction = max(betaEnv_reduction, na.rm = TRUE)
  )

wmel_effectiveness_summary
saveRDS(wmel_effectiveness_summary, paste0(OUT, "rq1_wmel_effectiveness_summary.rds"))

# ------------------------------------------------------------
# 3. Find homogeneous beta_mu reduction matching 89% dengue reduction
# ------------------------------------------------------------
beta_089 <- beta_calib_results %>%
  mutate(abs_error = abs(model_reduction - 89)) %>%
  slice_min(abs_error, n = 1, with_ties = FALSE) %>%
  pull(betaEnv_reduction)

beta_089


# ------------------------------------------------------------
# 4. Run model scenarios (paired with raw_default)
# ------------------------------------------------------------
raw_het <- run_raw_scenario(list(stim = stim_vec, betaEnv_reduction = betaEnv_reduction_vec))
raw_089 <- run_raw_scenario(list(stim = stim_vec, betaEnv_reduction = beta_089))
saveRDS(raw_het, paste0(OUT, "rq1_raw_het.rds"))
saveRDS(raw_089, paste0(OUT, "rq1_raw_089.rds"))


# ------------------------------------------------------------
# 5. Summarise outbreak metrics for table
# ------------------------------------------------------------

summarise_raw_outbreak <- function(raw_list, scenario) {

  bind_rows(lapply(raw_list, function(run) {

    daily_cases <- run$daily_cases

    # Peak day (first day of maximum incidence)
    peak_day <- which.max(daily_cases)

    # Peak number of daily cases
    peak_cases <- max(daily_cases, na.rm = TRUE)

    # Total outbreak size
    total_cases <- sum(daily_cases, na.rm = TRUE)

    # Duration = number of days from start of simulation until final reported case
    outbreak_duration_days <- max(which(daily_cases > 0))

    tibble(
      scenario = scenario,
      peak_day = peak_day,
      peak_cases = peak_cases,
      total_cases = total_cases,
      outbreak_duration_days = outbreak_duration_days
    )
  }))
}

scenario_runs <- bind_rows(
  summarise_raw_outbreak(raw_default, "No intervention"),
  summarise_raw_outbreak(raw_het, "Patch-based observed heterogeneous reduction"),
  summarise_raw_outbreak(raw_089, "Homogeneous 89% dengue reduction")
) %>%
  mutate(scenario = factor(scenario, levels = c("No intervention",
                                                "Patch-based observed heterogeneous reduction",
                                                "Homogeneous 89% dengue reduction")))

outbreak_table <- scenario_runs %>%
  group_by(scenario) %>%
  summarise(
    peak_day = paste0(
      median(peak_day, na.rm = TRUE), " (",
      quantile(peak_day, 0.025, na.rm = TRUE), "–",
      quantile(peak_day, 0.975, na.rm = TRUE), ")"
    ),
    peak_cases = paste0(
      round(median(peak_cases, na.rm = TRUE)), " (",
      round(quantile(peak_cases, 0.025, na.rm = TRUE)), "–",
      round(quantile(peak_cases, 0.975, na.rm = TRUE)), ")"
    ),
    total_cases = paste0(
      round(median(total_cases, na.rm = TRUE)), " (",
      round(quantile(total_cases, 0.025, na.rm = TRUE)), "–",
      round(quantile(total_cases, 0.975, na.rm = TRUE)), ")"
    ),
    outbreak_duration_days = paste0(
      median(outbreak_duration_days, na.rm = TRUE), " (",
      quantile(outbreak_duration_days, 0.025, na.rm = TRUE), "–",
      quantile(outbreak_duration_days, 0.975, na.rm = TRUE), ")"
    ),
    .groups = "drop"
  ) %>%
  arrange(scenario) %>%
  rename(
    `Wolbachia scenario` = scenario,
    `Day of peak cases` = peak_day,
    `Number of peak cases` = peak_cases,
    `Total cases across outbreak` = total_cases,
    `Outbreak duration (days)` = outbreak_duration_days
  )

outbreak_table
saveRDS(outbreak_table, paste0(OUT, "rq1_outbreak_table.rds"))

# ------------------------------------------------------------
# 6. Paired effectiveness and dataframes for mapping
# ------------------------------------------------------------
extract_patch_cases <- function(res) {
  if (!is.null(res$patch_cases)) return(res$patch_cases)
  if (!is.null(res$Patchtotals)) {
    newD_mat <- do.call(rbind, lapply(res$Patchtotals, function(x) x$newD))
    return(colSums(newD_mat, na.rm = TRUE))
  }
  stop("Run object contains neither `patch_cases` nor `Patchtotals`.")
}

# patch-level summary of cumulative cases across runs
summarise_patch_cases <- function(run_list, unipix) {
  patch_mat <- do.call(cbind, lapply(run_list, extract_patch_cases))
  stopifnot(nrow(patch_mat) == nrow(unipix))
  unipix %>%
    mutate(
      cum_cases_median = apply(patch_mat, 1, median, na.rm = TRUE),
      cum_cases_q025   = apply(patch_mat, 1, quantile, 0.025, na.rm = TRUE),
      cum_cases_q975   = apply(patch_mat, 1, quantile, 0.975, na.rm = TRUE),
      prob_affected    = rowMeans(patch_mat >= 1, na.rm = TRUE)
    )
}

patch_default <- summarise_patch_cases(raw_default, unipix) %>%
  rename(default_cases = cum_cases_median)

# median paired per-run effectiveness in each patch
patch_effectiveness <- patch_default %>%
  dplyr::select(patchID, x, y, pop, default_cases) %>%
  mutate(effectiveness_het = paired_patch_effectiveness(raw_het, raw_default),
         effectiveness_089 = paired_patch_effectiveness(raw_089, raw_default))

eff_limits <- range(c(patch_effectiveness$effectiveness_het, patch_effectiveness$effectiveness_089), na.rm = TRUE)

# scenario-level paired effectiveness (also used for RQ4 reference lines)
het_stats <- paired_stats(raw_het, raw_default, patch_pop)
hom_stats <- paired_stats(raw_089, raw_default, patch_pop)

effectiveness_summary <- bind_rows(
  "Patch-based observed heterogeneous reduction" = het_stats,
  "Homogeneous 89% dengue reduction"             = hom_stats,
  .id = "scenario") %>%
  group_by(scenario) %>%
  summarise(
    effectiveness_pop_weighted_median = median(effectiveness_pop_weighted, na.rm = TRUE),
    effectiveness_pop_weighted_q025   = quantile(effectiveness_pop_weighted, 0.025, na.rm = TRUE),
    effectiveness_pop_weighted_q975   = quantile(effectiveness_pop_weighted, 0.975, na.rm = TRUE),
    effectiveness_unweighted_median   = median(effectiveness_unweighted, na.rm = TRUE),
    effectiveness_unweighted_q025     = quantile(effectiveness_unweighted, 0.025, na.rm = TRUE),
    effectiveness_unweighted_q975     = quantile(effectiveness_unweighted, 0.975, na.rm = TRUE),
    effectiveness_total_median        = median(effectiveness_total, na.rm = TRUE),
    effectiveness_total_q025          = quantile(effectiveness_total, 0.025, na.rm = TRUE),
    effectiveness_total_q975          = quantile(effectiveness_total, 0.975, na.rm = TRUE),
    cases_averted_median              = median(cases_averted, na.rm = TRUE),
    cases_averted_q025                = quantile(cases_averted, 0.025, na.rm = TRUE),
    cases_averted_q975                = quantile(cases_averted, 0.975, na.rm = TRUE),
    .groups = "drop")

effectiveness_summary %>% as.data.frame()
saveRDS(effectiveness_summary, paste0(OUT, "rq1_effectiveness_summary.rds"))

# indirect protection in the observed landscape: patches with vs without established Wolbachia
rq1_untreated_summary <- patch_effectiveness %>%
  mutate(wmel_prev = wmel_patch$wmel_prev) %>%
  group_by(no_wolbachia = wmel_prev == 0) %>%
  summarise(n = n(),
            median_patch_eff = median(effectiveness_het, na.rm = TRUE),
            pop_weighted_eff = weighted.mean(effectiveness_het, pop, na.rm = TRUE),
            .groups = "drop")
rq1_untreated_summary
saveRDS(rq1_untreated_summary, paste0(OUT, "rq1_untreated_summary.rds"))

# ------------------------------------------------------------
# 7. Map no-intervention cases and modelled Wolbachia effectiveness
# ------------------------------------------------------------

map_A <- ggplot(patch_effectiveness %>%
                  mutate(default_cases_plot = ifelse(default_cases == 0, NA_real_, default_cases)),
                aes(x = x, y = y, fill = default_cases_plot)) +
  geom_tile(width = 0.005, height = 0.005) +
  coord_equal() +
  scale_fill_viridis_c(name = "Cumulative\ncases",
                       #option = "C",
                       begin = 0.15,
                       end = 0.95,
                       na.value = "grey90") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

map_B <- ggplot(patch_effectiveness, aes(x = x, y = y, fill = effectiveness_het)) +
  geom_tile(width = 0.005, height = 0.005) +
  coord_equal() +
  scale_fill_viridis_c(
    name = "Modelled\neffectiveness",
    begin = 0.15,
    end = 0.95,
    labels = scales::percent_format(),
    limits = eff_limits,
    na.value = "grey90"
  ) +
  labs(x = "Longitude", y = "") +
  theme_bw()

map_C <- ggplot(patch_effectiveness, aes(x = x, y = y, fill = effectiveness_089)) +
  geom_tile(width = 0.005, height = 0.005) +
  coord_equal() +
  scale_fill_viridis_c(
    name = "Modelled\neffectiveness",
    begin = 0.15,
    end = 0.95,
    labels = scales::percent_format(),
    limits = eff_limits,
    na.value = "grey90"
  ) +
  labs(x = "Longitude", y = "") +
  theme_bw()

map_A_final <- map_A +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.theme = element_text(margin = margin(r = 10)))
  ) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(margin = margin(b = 12)),
        legend.margin = margin(t = -5, b = -5),
        legend.box.margin = margin(t = 0),
        plot.margin = margin(2, 5, 2, 2))

map_B_final <- map_B +
  guides(fill = guide_colorbar(
    title.position = "left",
    title.theme = element_text(margin = margin(r = 10)))
  ) +
  theme(legend.position = "bottom",
        legend.direction = "horizontal",
        legend.title = element_text(margin = margin(b = 12)),
        legend.margin = margin(t = -5, b = -5),
        legend.box.margin = margin(l = 68),
        plot.margin = margin(2, 5, 2, 2))

map_C_final <- map_C +
  theme(legend.position = "none",
        plot.margin = margin(2, 2, 2, 2))

wmel_maps_panel <- cowplot::plot_grid(
  map_A_final,
  map_B_final,
  map_C_final,
  labels = c("A", "B", "C"),
  label_fontface = "bold",
  label_size = 16,
  nrow = 1,
  align = "hv"
)

wmel_maps_panel
ggsave(paste0(OUT, "FIGURE_wolbachia_effectiveness_maps.png"), wmel_maps_panel, width = 14, height = 3.25, dpi = 300)



#### 2. RQ2 levels of indirect coverage #### ---------------------------------------------------------
# ------------------------------------------------------------
# 1. Setup
# ------------------------------------------------------------
wol_coverage_grid   <- seq(0.05, 0.95, by = 0.05)
stopifnot(length(raw_default) >= N_RUNS)

# check observed number of zero Wolbachia patches for reporting
n_observed_zero_wol <- sum(wmel_patch$wmel_prev == 0, na.rm = TRUE)
n_observed_zero_wol
n_observed_zero_wol / nrow(unipix) * 100 # % of patches with 0% Wolbachia in observed data

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
# which patches are untreated: own seed, independent of coverage and of the model run,
# so each run uses the same untreated patches at every coverage level
untreated_patches <- function(n_zero, layout, run) {
  set.seed(1e6 + n_zero * 1e3 + run * 2 + (layout == "clustered"))
  if (layout == "random") return(sample(nrow(unipix), n_zero))
  centre <- sample(nrow(unipix), 1)
  order(pixdistmat[centre, ])[seq_len(n_zero)]     # n nearest patches = one contiguous gap
}

# clustering index: mean share of each untreated patch's neighbours (<= 750 m) also untreated
neighbour_share <- function(idx, radius_m = 750) {
  if (length(idx) < 2) return(NA_real_)
  vals <- sapply(idx, function(k) {
    nb <- which(pixdistmat[k, ] > 0 & pixdistmat[k, ] <= radius_m)
    if (length(nb) == 0) NA_real_ else mean(nb %in% idx)
  })
  mean(vals, na.rm = TRUE)
}

# one paired indirect-protection run
run_indirect_once <- function(wol_coverage, n_zero, layout, run) {
  zero_idx  <- untreated_patches(n_zero, layout, run)
  beta_full <- rep(beta_map(wol_coverage), nrow(unipix))
  beta_full[zero_idx] <- 0

  set.seed(SEED_BASE + run)                 # same model random stream as baseline run `run`
  res <- DEN.spatial(timestep = timestep, fitdat = fitdat, pastdat = pastdat,
                     sgpop = riopop, unipix = unipix, pixdistmat = pixdistmat,
                     steprun = n_days,
                     paramsList = list(stim = stim_vec, betaEnv_reduction = beta_full))
  if (!is.list(res)) return(NULL)
  newD_mat <- do.call(rbind, lapply(res$Patchtotals, function(x) x$newD))
  x <- colSums(newD_mat, na.rm = TRUE)

  base <- raw_default[[run]]                # paired no-intervention run
  stopifnot(base$run == run)
  b <- base$patch_cases

  pe <- ifelse(b[zero_idx] > 0, 1 - x[zero_idx] / b[zero_idx], NA_real_)
  pe <- pmin(pmax(pe, 0), 1)

  tibble(
    wolbachia_coverage = wol_coverage, n_zero_patches = n_zero, layout = layout, run = run,
    clustering = neighbour_share(zero_idx),
    indirect_effectiveness_pop_weighted   = weighted.mean(pe, patch_pop[zero_idx], na.rm = TRUE),
    indirect_effectiveness_unweighted     = mean(pe, na.rm = TRUE),
    indirect_effectiveness_case_aggregate = 1 - sum(x[zero_idx]) / sum(b[zero_idx]),
    untreated_cases          = sum(x[zero_idx]),
    untreated_baseline_cases = sum(b[zero_idx])
  )
}

# ------------------------------------------------------------
# 3. Clustering reference: observed zero-Wolbachia patches
# ------------------------------------------------------------
observed_zero_idx <- which(wmel_patch$wmel_prev == 0)
length(observed_zero_idx)
observed_clustering <- neighbour_share(observed_zero_idx)
observed_clustering

# ------------------------------------------------------------
# 4. Run RQ2 grid
# ------------------------------------------------------------
rq2_configs <- bind_rows(
  expand.grid(wolbachia_coverage = wol_coverage_grid, n_zero = 1,
              layout = "random", stringsAsFactors = FALSE),
  expand.grid(wolbachia_coverage = wol_coverage_grid, n_zero = c(10, 50, 100),
              layout = c("random", "clustered"), stringsAsFactors = FALSE)
)

rq2_results <- vector("list", nrow(rq2_configs))
for (k in seq_len(nrow(rq2_configs))) {
  cfg <- rq2_configs[k, ]
  message("RQ2 ", k, "/", nrow(rq2_configs), ": cov = ", cfg$wolbachia_coverage,
          ", n = ", cfg$n_zero, ", ", cfg$layout)
  out <- parallel::mclapply(
    seq_len(N_RUNS),
    function(r) run_indirect_once(cfg$wolbachia_coverage, cfg$n_zero, cfg$layout, r),
    mc.cores = 8, mc.preschedule = FALSE)
  ok <- vapply(out, function(o) is.data.frame(o), logical(1))
  if (any(!ok)) warning("RQ2 config ", k, ": ", sum(!ok), " runs failed and were dropped")
  rq2_results[[k]] <- bind_rows(out[ok])
  if (k %% 10 == 0) saveRDS(rq2_results, paste0(OUT, "rq2_indirect_protection_in_progress.rds"))
}
rq2_results <- bind_rows(rq2_results)
saveRDS(rq2_results, paste0(OUT, "rq2_indirect_protection_raw.rds"))

# ------------------------------------------------------------
# 5. Summarise across stochastic runs
# ------------------------------------------------------------
rq2_summary <- rq2_results %>%
  group_by(wolbachia_coverage, n_zero_patches, layout) %>%
  summarise(effectiveness_median = median(indirect_effectiveness_pop_weighted, na.rm = TRUE),
            effectiveness_q025 = quantile(indirect_effectiveness_pop_weighted, 0.025, na.rm = TRUE),
            effectiveness_q975 = quantile(indirect_effectiveness_pop_weighted, 0.975, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(percent_patches = 100 * n_zero_patches / nrow(unipix),
         n_zero_label = paste0(n_zero_patches,
                               ifelse(n_zero_patches == 1, " untreated patch (", " untreated patches ("),
                               round(percent_patches, 1), "%)"),
         n_zero_label = factor(n_zero_label, levels = unique(n_zero_label[order(n_zero_patches)])),
         layout = factor(layout, levels = c("random", "clustered"), labels = c("Random", "Clustered")))
rq2_summary
saveRDS(rq2_summary, paste0(OUT, "rq2_summary.rds"))

# clustering achieved by each placement (same untreated set at every coverage level per run)
rq2_clustering <- rq2_results %>%
  distinct(n_zero_patches, layout, run, clustering) %>%
  group_by(n_zero_patches, layout) %>%
  summarise(clustering_median = median(clustering, na.rm = TRUE),
            clustering_q025 = quantile(clustering, 0.025, na.rm = TRUE),
            clustering_q975 = quantile(clustering, 0.975, na.rm = TRUE),
            .groups = "drop")
rq2_clustering
observed_clustering   # where the real deployment gaps sit

rq2_text_values <- rq2_summary %>%
  group_by(n_zero_patches, n_zero_label, layout) %>%
  summarise(
    # overall range of medians across all Wolbachia coverage levels
    median_min = min(effectiveness_median, na.rm = TRUE),
    median_max = max(effectiveness_median, na.rm = TRUE),

    # low-coverage behaviour: 5-20% Wolbachia
    low_cov_median = median(effectiveness_median[wolbachia_coverage <= 0.20], na.rm = TRUE),
    low_cov_q025 = median(effectiveness_q025[wolbachia_coverage <= 0.20], na.rm = TRUE),
    low_cov_q975 = median(effectiveness_q975[wolbachia_coverage <= 0.20], na.rm = TRUE),

    # high-coverage behaviour: 75-95% Wolbachia
    high_cov_median = median(effectiveness_median[wolbachia_coverage >= 0.75], na.rm = TRUE),
    high_cov_q025 = median(effectiveness_q025[wolbachia_coverage >= 0.75], na.rm = TRUE),
    high_cov_q975 = median(effectiveness_q975[wolbachia_coverage >= 0.75], na.rm = TRUE),

    # change from lowest to highest Wolbachia coverage
    effect_at_5pct = effectiveness_median[wolbachia_coverage == 0.05],
    effect_at_95pct = effectiveness_median[wolbachia_coverage == 0.95],
    change_5_to_95 = effect_at_95pct - effect_at_5pct,

    # uncertainty across all scenarios
    min_q025 = min(effectiveness_q025, na.rm = TRUE),
    max_q975 = max(effectiveness_q975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    median_range_pct = paste0(round(100 * median_min, 1), "–", round(100 * median_max, 1), "%"),
    low_cov_median_pct = round(100 * low_cov_median, 1),
    high_cov_median_pct = round(100 * high_cov_median, 1),
    effect_at_5pct_pct = round(100 * effect_at_5pct, 1),
    effect_at_95pct_pct = round(100 * effect_at_95pct, 1),
    change_5_to_95_pct = round(100 * change_5_to_95, 1),
    uncertainty_range_pct = paste0(round(100 * min_q025, 1), "–", round(100 * max_q975, 1), "%")
  )

rq2_text_values %>% as.data.frame() # values used in manuscript results
saveRDS(rq2_text_values, paste0(OUT, "rq2_text_values.rds"))


# ------------------------------------------------------------
# 6. Plot
# ------------------------------------------------------------

# area_coverage is treated-patch coverage × the share of patches treated. It’s used as the x variable and in the axis label.
# hom_ref is the homogeneous calibration curve. It’s drawn in every panel as a dotted line, with its own legend entry (“Homogeneous coverage”)

# area-wide coverage: treated-patch coverage diluted by the untreated patches
rq2_summary <- rq2_summary %>%
  mutate(area_coverage = wolbachia_coverage * (1 - n_zero_patches / nrow(unipix)))

# reference: homogeneous coverage of all patches (calibration runs)
hom_ref <- beta_calibration_lookup %>%
  transmute(area_coverage = wolbachia_prevalence,
            effectiveness_median = model_reduction / 100) %>%
  filter(area_coverage <= 0.95)

dodge <- position_dodge(width = 0.02)

rq2_plot <- ggplot(
  rq2_summary,
  aes(x = area_coverage, y = effectiveness_median, colour = layout, shape = layout)) +
  geom_line(data = hom_ref,
            aes(x = area_coverage, y = effectiveness_median, linetype = "Homogeneous coverage"),
            inherit.aes = FALSE, colour = "black", linewidth = 0.6) +
  geom_errorbar(aes(ymin = effectiveness_q025, ymax = effectiveness_q975),
                width = 0.015, linewidth = 0.4, position = dodge) +
  geom_line(linewidth = 0.7, position = dodge) +
  geom_point(size = 2, position = dodge) +
  facet_wrap(~ n_zero_label, ncol = 2) +
  scale_colour_manual(values = c("Random" = "black", "Clustered" = "grey55")) +
  scale_shape_manual(values = c("Random" = 16, "Clustered" = 17)) +
  scale_linetype_manual(values = c("Homogeneous coverage" = "dotted"), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1)) +
  labs(x = expression("Area-wide " * italic("Wolbachia") * " coverage"),
       y = "Indirect protection effectiveness",
       colour = "Untreated patch placement", shape = "Untreated patch placement") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(colour = "white"),
    panel.spacing = unit(1.4, "lines"),
    legend.position = "bottom"
  )

rq2_plot
ggsave(paste0(OUT, "FIGURE_RQ2_indirect_protection.png"), rq2_plot, width = 10, height = 8.5, dpi = 300)



#### 3. RQ3 impact of spatial heterogeneity #### ---------------------------------------------------------
# ------------------------------------------------------------
# 1. Helper functions
# ------------------------------------------------------------
summarise_effectiveness_rq3 <- function(raw_list, scenario, target, achieved) {
  s <- paired_stats(raw_list, raw_default, patch_pop)
  tibble(
    scenario = scenario,
    target_coverage = target,
    wolbachia_coverage = achieved,   # matched achieved mean used for comparison and plotting
    total_cases_median = median(s$total_cases, na.rm = TRUE),
    total_cases_q025 = quantile(s$total_cases, 0.025, na.rm = TRUE),
    total_cases_q975 = quantile(s$total_cases, 0.975, na.rm = TRUE),
    effectiveness_pop_weighted_median = median(s$effectiveness_pop_weighted, na.rm = TRUE),
    effectiveness_pop_weighted_q025 = quantile(s$effectiveness_pop_weighted, 0.025, na.rm = TRUE),
    effectiveness_pop_weighted_q975 = quantile(s$effectiveness_pop_weighted, 0.975, na.rm = TRUE),
    effectiveness_unweighted_median = median(s$effectiveness_unweighted, na.rm = TRUE),
    effectiveness_unweighted_q025 = quantile(s$effectiveness_unweighted, 0.025, na.rm = TRUE),
    effectiveness_unweighted_q975 = quantile(s$effectiveness_unweighted, 0.975, na.rm = TRUE)
  )
}

# Rescale the observed spatial Wolbachia pattern towards a population-weighted target
# mean, preserving spatial order. Patches with 0% stay at 0% and values are capped at
# 100%, so at high targets the achieved mean falls short of the target; the matched
# homogeneous scenario therefore uses the achieved mean, not the target.
make_heterogeneous_wmel <- function(wmel_prev, target_mean, pop) {
  current_mean <- weighted.mean(wmel_prev, pop, na.rm = TRUE)
  pmin(pmax(wmel_prev * (target_mean / current_mean), 0), 1)
}

# ------------------------------------------------------------
# 2. Run heterogeneous and matched homogeneous scenarios (paired)
# ------------------------------------------------------------
wol_coverage_grid <- seq(0.05, 0.95, by = 0.05)
wmel_prev_vec <- wmel_patch$wmel_prev   # unipix (patchID) order

rq3_results <- vector("list", length(wol_coverage_grid) * 2)
counter <- 1

for (cov in wol_coverage_grid) {

  # Heterogeneous scenario
  wmel_heterogeneous <- make_heterogeneous_wmel(wmel_prev_vec, cov, patch_pop)
  achieved <- weighted.mean(wmel_heterogeneous, patch_pop)
  message("RQ3 target = ", cov, " | achieved population-weighted mean = ", round(achieved, 3))

  raw_heterogeneous <- run_raw_scenario(
    list(stim = stim_vec, betaEnv_reduction = beta_map(wmel_heterogeneous)))
  rq3_results[[counter]] <- summarise_effectiveness_rq3(raw_heterogeneous, "Heterogeneous", cov, achieved)
  counter <- counter + 1
  rm(raw_heterogeneous); gc()

  # Homogeneous scenario at the same achieved mean
  raw_homogeneous <- run_raw_scenario(
    list(stim = stim_vec, betaEnv_reduction = beta_map(rep(achieved, nrow(unipix)))))
  rq3_results[[counter]] <- summarise_effectiveness_rq3(raw_homogeneous, "Homogeneous", cov, achieved)
  counter <- counter + 1
  rm(raw_homogeneous); gc()

  saveRDS(rq3_results, paste0(OUT, "rq3_in_progress.rds"))
}

rq3_summary <- bind_rows(rq3_results)
rq3_summary %>% distinct(target_coverage, wolbachia_coverage)   # target vs achieved mean
saveRDS(rq3_summary, paste0(OUT, "rq3_heterogeneous_vs_homogeneous_summary.rds"))

# ------------------------------------------------------------
# 3. Heterogeneity summary results table
# ------------------------------------------------------------
# calculate the epidemiological cost of heterogeneity
rq3_gap_summary <- rq3_summary %>%
  dplyr::select(scenario, target_coverage, wolbachia_coverage, effectiveness_pop_weighted_median) %>%
  tidyr::pivot_wider(names_from = scenario,
                     values_from = effectiveness_pop_weighted_median) %>%
  mutate(absolute_gap_pop_weighted = Homogeneous - Heterogeneous,
         relative_loss_pop_weighted = (Homogeneous - Heterogeneous) / Homogeneous)
rq3_gap_summary %>% as.data.frame()
saveRDS(rq3_gap_summary, paste0(OUT, "rq3_gap_summary.rds"))

# summarise across all coverage levels
rq3_gap_stats <- rq3_gap_summary %>%
  summarise(median_absolute_gap_pop_weighted = median(absolute_gap_pop_weighted, na.rm = TRUE),
            q25_absolute_loss_pop_weighted = quantile(absolute_gap_pop_weighted, 0.25, na.rm = TRUE),
            q75_absolute_loss_pop_weighted = quantile(absolute_gap_pop_weighted, 0.75, na.rm = TRUE),
            min_absolute_gap_pop_weighted = min(absolute_gap_pop_weighted, na.rm = TRUE),
            max_absolute_gap_pop_weighted = max(absolute_gap_pop_weighted, na.rm = TRUE),
            median_relative_loss_pop_weighted = median(relative_loss_pop_weighted, na.rm = TRUE),
            q25_relative_loss_pop_weighted = quantile(relative_loss_pop_weighted, 0.25, na.rm = TRUE),
            q75_relative_loss_pop_weighted = quantile(relative_loss_pop_weighted, 0.75, na.rm = TRUE),
            min_relative_loss_pop_weighted = min(relative_loss_pop_weighted, na.rm = TRUE),
            max_relative_loss_pop_weighted = max(relative_loss_pop_weighted, na.rm = TRUE))

rq3_gap_stats %>% as.data.frame()
saveRDS(rq3_gap_stats, paste0(OUT, "rq3_gap_stats.rds"))


# ------------------------------------------------------------
# 4. Plot
# ------------------------------------------------------------
rq3_summary <- rq3_summary %>%
  mutate(scenario = factor(scenario, levels = c("Homogeneous", "Heterogeneous")))

rq3_plot <- ggplot(
  rq3_summary,
  aes(x = wolbachia_coverage, y = effectiveness_pop_weighted_median,
      linetype = scenario, group = scenario)) +
  geom_errorbar(aes(ymin = pmax(effectiveness_pop_weighted_q025, 0), ymax = pmin(effectiveness_pop_weighted_q975, 1)),
                width = 0.008, linewidth = 0.35, position = position_dodge(width = 0.02), color = "grey50") +
  geom_line(linewidth = 0.9, position = position_dodge(width = 0.02)) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.2),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_linetype_manual(values = c("Homogeneous" = "solid", "Heterogeneous" = "dashed"),
                        labels = c("Homogeneous" = "Homogeneous\ncoverage", "Heterogeneous" = "Heterogeneous\ncoverage")) +
  labs(x = expression("Population-weighted mean " * italic("Wolbachia") * " prevalence"), y = "Modelled effectiveness",
       colour = NULL, fill = NULL, linetype = NULL) +
  theme_bw() +
  theme(
    legend.position = c(0.99, 0.01),
    legend.justification = c(1, 0),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = alpha("white", 0.8), colour = "black", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

rq3_plot
ggsave(paste0(OUT, "FIGURE_RQ3_heterogeneous_vs_homogeneous_wolbachia.png"), rq3_plot, width = 10, height = 6, dpi = 300)



#### 4. RQ4 reactive CATI compensating wol gaps #### ---------------------------------------------------------
# ------------------------------------------------------------
# 1. Scenario grid
# ------------------------------------------------------------
radii <- seq(200, 500, by = 50)
eff_vals <- c(0.1, 0.3, 0.5)
delay_vals <- c(1, 2, 3)

cati_grid <- expand.grid(
  radius = radii,
  eff = eff_vals,
  delay = delay_vals
)

# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------
summarise_effectiveness_cati <- function(raw_list, scenario_info) {
  s <- paired_stats(raw_list, raw_default, patch_pop)
  tibble(
    radius = scenario_info$radius,
    eff = scenario_info$eff,
    delay = scenario_info$delay,

    effectiveness_unweighted_median = median(s$effectiveness_unweighted, na.rm = TRUE),
    effectiveness_unweighted_q025 = quantile(s$effectiveness_unweighted, 0.025, na.rm = TRUE),
    effectiveness_unweighted_q975 = quantile(s$effectiveness_unweighted, 0.975, na.rm = TRUE),

    effectiveness_pop_weighted_median = median(s$effectiveness_pop_weighted, na.rm = TRUE),
    effectiveness_pop_weighted_q025 = quantile(s$effectiveness_pop_weighted, 0.025, na.rm = TRUE),
    effectiveness_pop_weighted_q975 = quantile(s$effectiveness_pop_weighted, 0.975, na.rm = TRUE),

    total_cases_median = median(s$total_cases, na.rm = TRUE),
    total_cases_q025 = quantile(s$total_cases, 0.025, na.rm = TRUE),
    total_cases_q975 = quantile(s$total_cases, 0.975, na.rm = TRUE),

    # paired cases averted relative to no intervention
    cases_averted_median = median(s$cases_averted, na.rm = TRUE),
    cases_averted_q025 = quantile(s$cases_averted, 0.025, na.rm = TRUE),
    cases_averted_q975 = quantile(s$cases_averted, 0.975, na.rm = TRUE)
  )
}

# Area of a circle of radius r (centred at 0,0) falling inside rectangle [x1,x2] x [y1,y2]
circle_rect_area <- function(r, x1, x2, y1, y2) {
  f <- function(x) {
    h <- sqrt(pmax(r^2 - x^2, 0))                     # half-height of circle at x
    pmax(pmin(y2, h) - pmax(y1, -h), 0)
  }
  lo <- max(x1, -r); hi <- min(x2, r)
  if (lo >= hi) return(0)
  integrate(f, lo, hi)$value
}

make_vectreat_coverage <- function(eff, radius_m, patch_size_m = 500,
                                   duration = 1, delay = 0) {
  # Case assumed at the centre of its patch. Each patch's coverage is the share of
  # its area inside a circle of radius_m: index patch, 4 edge and 4 corner neighbours.
  # Cov_* columns require the modified DEN.spatial(); Duration is not used by DEN.spatial(),
  # so each CATI application lasts one day per detection event, after the delay.
  h <- patch_size_m / 2
  a <- patch_size_m^2
  cov_index  <- circle_rect_area(radius_m, -h, h, -h, h) / a
  cov_edge   <- circle_rect_area(radius_m,  h, 3 * h, -h, h) / a
  cov_corner <- circle_rect_area(radius_m,  h, 3 * h,  h, 3 * h) / a

  data.frame(EffCoverage = eff,          # efficacy only; scaled per patch inside DEN.spatial
             Duration    = duration,
             Radius      = 750,          # reach all 8 neighbours (edge ~500 m, corner ~707 m)
             Delay       = delay,
             Cov_index   = cov_index,
             Cov_edge    = cov_edge,
             Cov_corner  = cov_corner)
}

# check this reproduces the table of proportional neighbouring patch coverage we expect
t(sapply(seq(200, 500, 50), function(L)
  unlist(make_vectreat_coverage(0.3, L)[c("Cov_index", "Cov_edge", "Cov_corner")])))

run_cati_one <- function(radius, eff, delay) {
  message("Running CATI: radius = ", radius, ", eff = ", eff, ", delay = ", delay)
  vectreat <- make_vectreat_coverage(eff = eff, radius_m = radius, patch_size_m = 500,
                                     duration = 1, delay = delay)
  run_raw_scenario(list(stim = stim_vec,
                        betaEnv_reduction = betaEnv_reduction_vec,
                        vectreat = vectreat))
}

# ------------------------------------------------------------
# 3. Baseline reference scenarios (paired results from RQ1)
# ------------------------------------------------------------
ref_lines <- tibble(
  scenario = c(
    "Observed heterogeneous Wolbachia",
    "Homogeneous high-coverage scenario"
  ),
  effectiveness_pop_weighted = c(
    median(het_stats$effectiveness_pop_weighted, na.rm = TRUE),
    median(hom_stats$effectiveness_pop_weighted, na.rm = TRUE)
  ),
  ref_linetype = c("dashed", "dotdash")
)
ref_lines

# ------------------------------------------------------------
# 4. Run CATI grid
# ------------------------------------------------------------
cati_results <- vector("list", nrow(cati_grid))

for (i in seq_len(nrow(cati_grid))) {
  combo_row <- cati_grid[i, ]
  raw_cati <- run_cati_one(radius = combo_row$radius, eff = combo_row$eff, delay = combo_row$delay)
  cati_results[[i]] <- summarise_effectiveness_cati(raw_cati, combo_row)
  rm(raw_cati); gc()
  if (i %% 10 == 0) saveRDS(cati_results, paste0(OUT, "rq4_in_progress.rds"))
}

cati_summary <- bind_rows(cati_results)
saveRDS(cati_summary, paste0(OUT, "rq4_cati_effectiveness_summary.rds"))

# ------------------------------------------------------------
# 5. CATI summary results table
# ------------------------------------------------------------
cati_optimal_table <- cati_summary %>%
  mutate(cases_averted_n = cases_averted_median,
         eff_pct_median = 100 * effectiveness_pop_weighted_median,
         eff_pct_q025 = 100 * effectiveness_pop_weighted_q025,
         eff_pct_q975 = 100 * effectiveness_pop_weighted_q975) %>%
  group_by(eff, delay) %>%
  slice_max(effectiveness_pop_weighted_median, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(result = paste0(
    radius, " m: ",
    round(cases_averted_n), " (",
    round(eff_pct_median), "%, ",
    round(eff_pct_q025), "–",
    round(eff_pct_q975), "%)"
  )) %>%
  dplyr::select(eff, delay, result) %>%
  tidyr::pivot_wider(names_from = delay, values_from = result, names_prefix = "Delay_") %>%
  rename(
    `Adulticide efficacy` = eff,
    `1 day delay` = Delay_1,
    `2 day delay` = Delay_2,
    `3 day delay` = Delay_3
  )

cati_optimal_table
saveRDS(cati_optimal_table, paste0(OUT, "rq4_cati_optimal_table.rds"))


# ------------------------------------------------------------
# 6. Plot CATI effectiveness
# ------------------------------------------------------------
# plot df
cati_summary <- cati_summary %>%
  mutate(delay = factor(delay, levels = c(1, 2, 3)),
         eff_lab = factor(paste0("Adulticide efficacy = ", eff), levels = paste0("Adulticide efficacy = ", sort(unique(eff)))))

# plot settings to allow dodge
band_width <- 50
radius_vals <- sort(unique(cati_grid$radius))
band_df <- tibble(
  radius = radius_vals,
  xmin = radius_vals - band_width / 2,
  xmax = radius_vals + band_width / 2,
  ymin = -Inf,
  ymax = Inf,
  band = rep(c("grey95", "white"), length.out = length(radius_vals))
)

# plot
cati_effectiveness_plot <- ggplot(
  cati_summary, aes(x = radius, y = effectiveness_pop_weighted_median, colour = delay, group = delay)) +
  geom_rect(data = band_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = band),
            inherit.aes = FALSE, fill = NA, colour = NA, alpha = 1) +
  scale_fill_identity() +
  geom_hline(data = ref_lines,
             aes(yintercept = effectiveness_pop_weighted, linetype = ref_linetype),
             inherit.aes = FALSE, colour = "black",
             linewidth = 0.5, show.legend = FALSE) +
  scale_linetype_identity() +
  geom_errorbar(aes(ymin = effectiveness_pop_weighted_q025, ymax = effectiveness_pop_weighted_q975),
                width = 8, position = position_dodge(width = 20),
                linewidth = 0.4, alpha = 0.7) +
  geom_line(linewidth = 0.7, position = position_dodge(width = 20)) +
  geom_point(size = 2, position = position_dodge(width = 20)) +
  facet_wrap(~ eff_lab, nrow = 1) +
  scale_colour_manual(values = c("1" = "black",
                                 "2" = "grey40",
                                 "3" = "grey70")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 1, by = 0.25),
                     limits = c(0, 1),
                     expand = c(0, 0)) +
  scale_x_continuous(breaks = radii,
                     limits = c(180,520),
                     expand = c(0, 0)) +
  labs(x = "CATI radius (m)",
       y = "Modelled effectiveness",
       color = "Days delay") +
  theme_bw() +
  theme(
    legend.position = c(0.99, 0.02),
    legend.direction = "horizontal",
    legend.justification = c(1, 0),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.3),
    axis.ticks.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(fill = "black"),
    strip.text = element_text(colour = "white")
  )

cati_effectiveness_plot
ggsave(paste0(OUT, "FIGURE_cati_effectiveness_plot.png"), cati_effectiveness_plot, width = 10, height = 5, dpi = 300)




#### 5. Additional plots for manuscript #### ---------------------------------------------------------

#### a) Figure: model structure and immunity landscape, two panel

# ---- PANEL A: model structure figure ----

img_A <- png::readPNG("26-11-20_FIGURE_model_structure.png")
# Made outside R (PowerPoint), exported as PNG
panel_A <- ggdraw() +
  draw_image(img_A) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 8))

# ---- PANEL B: R immunity map plot ----

img_B <- png::readPNG("FIGURE_immunity_ladscape.png")
# Saved by Rio_refit_workflow.Rmd (Step 1.d, immunity landscape)
panel_B <- ggdraw() +
  draw_image(img_B) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 8))

# ---- COMBINE ----

methods_sero_plot <- plot_grid(
  panel_A,
  panel_B,
  labels = c("A", "B"),
  label_fontface = "bold",
  label_size = 16,
  label_x = 0.02,
  label_y = 0.98,
  hjust = 0,
  vjust = 1,
  ncol = 1,
  align = "v",
  rel_heights = c(1, 1)
)

print(methods_sero_plot)
ggsave(paste0(OUT, "FIGURE_model_structure_immunity_landscape_panel.png"), methods_sero_plot,
       width = 550 / 72, height = 520 / 72, units = "in", dpi = 300, bg = "white")


#### b) Figure: model fitting diagnostics, three panel

# Panels saved by Rio_refit_workflow.Rmd:
#   epidemic curve (chunk epidemic_curve_validation_manuscript),
#   lagged spatial spread CDF (Step 5.d), parameter histograms (chunk fig4_combined_panel)

img_A <- png::readPNG("FIGURE_epidemic_curve_validation.png")
panel_A <- ggdraw() +
  draw_image(img_A) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

img_B <- png::readPNG("FIGURE_lagged_spatial_spread_validation_cdf.png")
panel_B <- ggdraw() +
  draw_image(img_B) +
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0))

img_C <- png::readPNG("FIGURE_param_histograms.png")
panel_C <- ggdraw() +
  draw_image(img_C) +
  theme(plot.margin = margin(t = 18, r = 0, b = 0, l = 0))

top_row <- plot_grid(
  panel_A,
  panel_B,
  labels = c("A", "B"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 2,
  align = "h"
)

refit_params_plot <- plot_grid(
  top_row,
  panel_C,
  labels = c("", "C"),
  label_size = 18,
  label_fontface = "bold",
  ncol = 1,
  rel_heights = c(1.1, 0.9)
)

refit_params_plot
ggsave(paste0(OUT, "FIGURE_model_fit_panel.png"), refit_params_plot,
       width = 1400 / 72, height = 900 / 72, units = "in", dpi = 300, bg = "white")
