library(tidyverse);library(readxl);library(data.table);library(ggrepel)

source("~/PhD-thesis/Utils/Functions.R")
source("~/PhD-thesis/Utils/med.R")
source("~/PhD-thesis/Cefepime/cefepime_models.R")

set.seed(123)
# Loading dataset ---------------------------------------------------------

total_covariates <- lapply(cefepime_models, function(x) x$covariates) %>% unlist() %>% unique()
df <- read_excel("Cefepime/data/nantes.xlsx") |> 
  filter_out(is.na(ID)) |> 
  mutate(across(where(is.character), ~na_if(.x, "."))) |> 
  readr::type_convert() |> 
  mutate(
    DV = Y,
    DUR = INFUSION_DURATION,
    DIAL = DIALYSE,
    CVVHD = CRRT,
    CLCR_CG = CG,
    gender = if_else(SEX == 0, "female", "male", missing = NA_character_),
    WT     = WEIGHT, 
    height = HEIGHT, 
    TBW = WT,
    CREA_micromol_l = SCR,
    creatinine_mg_dl = CREA_micromol_l/88.4,
  ) |> 
  mutate(
    EVID = if_else(EVID == 3, 1, EVID)
  ) |> 
  rowwise() %>%
  mutate(
    IBW = calculate_ibw_safe(
      height_cm = height,
      gender    = gender
    ),
    LBW = calculate_lbw_safe(                  
      height_cm = height,
      weight_kg = WT,
      gender    = gender
    ),
    LBM = LBW,                            
    ABW = calculate_abw_safe(
      height_cm = height,
      weight    = WT,
      gender    = gender
    ),
    BSA = calculate_bsa_safe(
      weight_kg = WT,
      height_cm = height
    ),
    CREA_micromol_l = convert_creatinine_safe(
      serum_creatinine_mg_dl = creatinine_mg_dl
    ),
    CLCR_CG = calculate_crcl_safe(            
      age                    = AGE,
      weight_kg              = WT,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_TBW    = CLCR_CG,               
    CLCR        = CLCR_CG,               
    CLCR_CG_LBW = calculate_crcl_safe(
      age                    = AGE,
      weight_kg              = LBW,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_CG_IBW = calculate_crcl_safe(       
      age                    = AGE,
      weight_kg              = IBW,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_J = calculate_crcl_jelliffe_safe(
      age              = AGE,
      weight           = WT,
      serum_creatinine = creatinine_mg_dl,
      gender           = gender
    ),
    eGFR_CKD_EPI = calculate_egfr_safe(
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    eGFR_MDRD = calculate_egfr_mdrd_safe(     
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_70 = calculate_crcl_safe(       
      age                    = AGE,
      weight_kg              = 70,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
  ) %>%
  ungroup() |> 
  group_by(CENTER, ID) |> 
  mutate(ID = cur_group_id()) |> 
  ungroup()


# Fix last dose -----------------------------------------------------------

.n_before <- nrow(df)

df <- df %>%
  group_by(ID) %>%
  group_modify(function(g, key) {
    last_tdm <- suppressWarnings(max(g$TIME[g$EVID == 0 & !is.na(g$DV)]))
    if (!is.finite(last_tdm)) return(g)          # no TDM -> nothing to extend
    pad <- g |> 
      slice_max(TIME, n = 1, with_ties = FALSE) |> 
      mutate(TIME = last_tdm + 1, EVID = 1, AMT = 0, DUR = 0, DV = NA_real_)
    bind_rows(g, pad) |>  arrange(TIME)
  }) %>%
  ungroup()

message("Righe segnaposto aggiunte: ", nrow(df) - .n_before,
        " (una per paziente con almeno un TDM)")

df_tdm <- df |> 
  group_by(ID) |> 
  summarise(
    tdm = sum(EVID == 0 & !is.na(DV), na.rm = TRUE),
    time = max(TIME[EVID == 0], na.rm = TRUE),
    time2 = mean(TIME[EVID == 0], na.rm = TRUE),
    time3 = median(TIME[EVID == 0], na.rm = TRUE)
  )


# Mask events -------------------------------------------------------------

mask_tdm_after <- function(data, n_keep) {
  data |>
    arrange(ID, TIME) |>
    group_by(ID) |>
    mutate(
      .tdm_i = cumsum(EVID == 0 & !is.na(DV)),   # 1,2,3... over TDM rows only
      EVID   = replace(EVID, EVID == 0 & !is.na(DV) & .tdm_i > n_keep, 2)
    ) |>
    ungroup() 
}

# Create a dataframe where I mask the first TDM to predict the second
df_tdm1 <- mask_tdm_after(df, 1) |> 
  filter_out(.tdm_i == 3)

# Create a dataframe where I use the first 2 TDMs to predict the 3rd
df_tdm2 <-  mask_tdm_after(df, 2) |>  
  group_by(ID) |> 
  filter(any(.tdm_i == 3))


# Population evaluation ---------------------------------------------------

population_df <- df

# WHAT `Cc` IS, AND WHY RAISING pop_n_simul DOES NOT SILENTLY CHANGE IT
# Everything downstream -- rBias/rRMSE/MDPE/MDAPE/F20/F30, the metric panels,
# plot_ipred -- reads a single column called Cc and treats it as THE population
# prediction for that patient at that time. Two different quantities can fill
# that column and they are not interchangeable:
#
#   TYPICAL prediction, f(theta, eta = 0). What n_simul = 0 returns, and what
#     external-validation papers mean by PRED. rBias computed against it is the
#     number that is comparable with the published literature, and the +/-20%
#     acceptability band was defined for it.
#   MEAN or MEDIAN of the simulated population predictive. Because the model is
#     nonlinear in eta, neither equals f(eta = 0): with log-normal IIV on CL and
#     V the mean sits systematically ABOVE the typical profile (Jensen), so
#     swapping it in shifts rBias upward by an amount that has nothing to do
#     with how well the model fits.
#
# So pop_n_simul only ADDS columns (Cc_median, Cc_mean, Cc_p05, Cc_p95, and the
# AUC equivalents). Which of them becomes Cc is a separate, explicit choice made
# by pop_pred_stat, and it defaults to "typical" -- identical to the n_simul = 0
# behaviour this loop had before. Change it only deliberately, and say so in the
# write-up, because it redefines the estimand the metrics are scoring.
#
# The percentile columns are IIV-only: poso_simu_pop returns f, not Y, so no
# residual error is included. They are a parameter-uncertainty ribbon, NOT a
# VPC interval, and observations are expected to fall outside them more often
# than 10% of the time. Utils/med.R adds the residual error where a genuine
# predictive distribution is needed.
pop_n_simul   <- 1000       # 0 = typical profile only; > 0 also simulates
pop_pred_stat <- "typical"  # "typical" | "median" | "mean" -> which becomes Cc
pop_grid_n    <- 400        # cap on plotting-grid points when simulating

stopifnot(pop_pred_stat %in% c("typical", "median", "mean"))
if (pop_pred_stat != "typical" && pop_n_simul == 0) {
  stop("pop_pred_stat = '", pop_pred_stat, "' requires pop_n_simul > 0.")
}

results <- list()
m_out <- list()

for (i in sort(unique(population_df$ID))) {
  message("Processing item ", i, " of ", length(unique(population_df$ID)))

  df_i <- population_df |> filter(ID == i)
  max_time_i <- max(df_i$TIME, na.rm = TRUE)

  # The 0.1 h grid exists only so plot_cc can draw a smooth profile; the
  # observation times themselves come from as.et(df_i) and are exact whatever
  # the grid is. At n_simul = 0 that grid is free, but it is solved once PER
  # SIMULATION, and the longest follow-up here is 429 h = 4291 grid points, so
  # at 1000 simulations it would mean 4.3 million rows per patient per model
  # before any summarising. Capping the ADDED samples at pop_grid_n keeps the
  # curve smooth enough to plot without making the simulation infeasible; the
  # metrics are untouched either way because they only ever use observation
  # times.
  by_i <- if (pop_n_simul > 0) max(0.1, max_time_i / pop_grid_n) else 0.1
  et_i <- rxode2::as.et(df_i)
  et_i$add.sampling(seq(0, max_time_i, by = by_i))

  # Only EVID == 0 rows carry a DV in this dataset (checked: 89 of them, none
  # elsewhere), and no ID+TIME pair appears twice, so joining the observations
  # explicitly is equivalent to the previous join against all of df_i -- minus
  # the fill(ID) hack, since ID is simply set from the loop variable.
  obs_i <- df_i |> filter(EVID == 0, !is.na(DV)) |> select(TIME, DV)

  id_key <- as.character(i)
  results[[id_key]] <- list()

  for (model_name in names(cefepime_models)) {
    message("Testing model ", model_name, " for ID ", i)

    result_i <- tryCatch({
      model <- cefepime_models[[model_name]]

      solve_pop <- function(K) {
        s <- posologyr::poso_simu_pop(df_i, model, n_simul = K)
        s$model <- posologyr::poso_replace_et(
          target_model = s$model,
          prior_model  = model,
          event_table  = et_i
        )
        as.data.frame(s$model)
      }

      # eta = 0: one solve, no sim.id column, one row per time. The duplicate
      # times removed here are the event table's own (a dose record and a
      # sampling record landing on the same time), not repeated simulations.
      pred_i <- solve_pop(0) |>
        select(time, Cc_typical = Cc, AUC_typical = AUC) |>
        distinct(time, .keep_all = TRUE) |>
        arrange(time)

      if (pop_n_simul > 0) {
        # THE TRAP THIS REPLACES: with n_simul > 0 the solved object holds K
        # stacked profiles keyed by sim.id, so `distinct(time, .keep_all=TRUE)`
        # kept the FIRST simulation and discarded the other 999 -- silently, and
        # with no error, leaving the metrics scored against one arbitrary random
        # draw instead of the typical profile. Deduplication has to be by
        # (sim.id, time); aggregation is by time.
        sim_raw <- solve_pop(pop_n_simul)
        # rxode2 omits sim.id when it is handed exactly ONE parameter row, so
        # pop_n_simul = 1 comes back shaped like the eta = 0 solve and the
        # select() below would fail for every patient and every model. Nobody
        # would choose K = 1 deliberately, but it is reachable by typo and the
        # failure lands in the tryCatch, so it would look like the models broke.
        if (!"sim.id" %in% names(sim_raw)) sim_raw$sim.id <- 1L

        pred_i <- sim_raw |>
          select(sim.id, time, Cc, AUC) |>
          distinct(sim.id, time, .keep_all = TRUE) |>
          group_by(time) |>
          summarise(
            Cc_median  = median(Cc),
            Cc_mean    = mean(Cc),
            Cc_p05     = unname(quantile(Cc, 0.05)),
            Cc_p95     = unname(quantile(Cc, 0.95)),
            AUC_median = median(AUC),
            AUC_mean   = mean(AUC),
            .groups    = "drop"
          ) |>
          right_join(pred_i, by = "time") |>
          arrange(time)
      }

      pred_i |>
        mutate(
          Cc  = .data[[paste0("Cc_",  pop_pred_stat)]],
          AUC = .data[[paste0("AUC_", pop_pred_stat)]],
          ID  = i
        ) |>
        left_join(obs_i, by = c("time" = "TIME"))

    }, error = function(e) {
      message("  Model ", model_name, " failed for ID ", i, ": ", conditionMessage(e))
      NULL
    })

    # id_key, not i: the IDs in this cohort run 1-127 for 39 patients, so
    # results[[i]] indexed by POSITION and wrote patient 127 into slot 127 of a
    # 39-element list, padding it with NULLs and leaving the named slot created
    # two lines above permanently empty. rbind() dropped the NULLs so the final
    # table came out right by luck; this makes it right by construction.
    results[[id_key]][[model_name]] <- result_i
  }
}

combined_results <- list()
for (model_name in names(cefepime_models)) {
  combined_model_data <- do.call(rbind, lapply(results, `[[`, model_name))
  combined_results[[model_name]] <- combined_model_data
}

results_metrics <- lapply(combined_results, function(model_data) {
  model_data <- model_data %>% rename(predicted = Cc, observed = DV)
  model_data <- model_data %>% filter(!is.na(predicted) & !is.na(observed))
  calculate_metrics(model_data, n_bootstrap = 10000, conf_level = 0.95,
                    id_col = "ID")   # cluster bootstrap: resample patients, not samples
})

# The performance of the models was considered clinically acceptable if the rBias 
# was between −20% and 20%, with the 95% confidence intervals (CI) including zero.
# Additionally, the precision metric (rRMSE) should be as small as possible.

results_df <- bind_rows(
  lapply(names(results_metrics), function(model_name) {
    metrics <- results_metrics[[model_name]]
    
    data.frame(
      Model = model_name,
      rBias = metrics$rBias,
      rBias_Lower_CI = metrics[["rBias_CI"]][["lower"]],
      rBias_Upper_CI = metrics[["rBias_CI"]][["upper"]],
      rRMSE = metrics$rRMSE,
      rRMSE_Lower_CI = metrics[["rRMSE_CI"]][["lower"]],
      rRMSE_Upper_CI = metrics[["rRMSE_CI"]][["upper"]],
      MDPE = metrics$MDPE,
      MDPE_Lower_CI = metrics[["MDPE_CI"]][["lower"]],
      MDPE_Upper_CI = metrics[["MDPE_CI"]][["upper"]],
      MDAPE = metrics$MDAPE,
      MDAPE_Lower_CI = metrics[["MDAPE_CI"]][["lower"]],
      MDAPE_Upper_CI = metrics[["MDAPE_CI"]][["upper"]],
      F20 = metrics$F20,
      F20_Lower_CI = metrics[["F20_CI"]][["lower"]],
      F20_Upper_CI = metrics[["F20_CI"]][["upper"]],
      F30 = metrics$F30,
      F30_Lower_CI = metrics[["F30_CI"]][["lower"]],
      F30_Upper_CI = metrics[["F30_CI"]][["upper"]]
    )
  })
)

rownames(results_df) <- NULL

results_df <-results_df |>  
  mutate(
    across(
      rBias:F30_Upper_CI, ~ round(., 2)
    )
  )


# Population metric plots -------------------------------------------------

# One model order, reused by every panel, so a model sits in the same row in
# all of them and the three plots can be read side by side. Ordered by rRMSE
# (the precision metric), best at the top: with coord_flip() the FIRST factor
# level is drawn at the bottom, hence the descending sort.
model_levels <- results_df |>  arrange(desc(rRMSE)) |>  pull(Model)

results_df <- results_df |>  
  mutate(
    Model = factor(Model, levels = model_levels),
    # The acceptability rule stated above: rBias inside +/-20% AND a 95% CI
    # that covers zero. Both halves matter - a model can sit near zero with a
    # CI so wide it is uninformative, or be tightly estimated but off target.
    bias_ok = abs(rBias) <= 20 & rBias_Lower_CI <= 0 & rBias_Upper_CI >= 0
  )

theme_metrics <- theme_classic(base_size = 15) +
  theme(
    plot.title    = element_text(size = 17, face = "bold"),
    plot.subtitle = element_text(size = 12, colour = "grey30"),
    axis.text     = element_text(size = 13),
    axis.title    = element_text(size = 14),
    legend.position = "top",
    panel.grid.major.x = element_line(colour = "grey92")
  )

pass_cols <- c(`TRUE` = "#1b7837", `FALSE` = "#b2182b")

plot_rBias <- ggplot(results_df, aes(x = Model, y = rBias)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -20, ymax = 20,
           fill = "grey60", alpha = 0.18) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey20") +
  geom_errorbar(aes(ymin = rBias_Lower_CI, ymax = rBias_Upper_CI,
                    colour = bias_ok),
                width = 0.25, linewidth = 0.7) +
  geom_point(aes(colour = bias_ok), size = 3.2) +
  geom_text(aes(y = rBias_Upper_CI, label = sprintf("%+.1f", rBias)),
            hjust = -0.4, size = 4, colour = "grey20") +
  scale_colour_manual(
    values = pass_cols,
    breaks = c("TRUE", "FALSE"),
    labels = c("meets criterion", "fails criterion"),
    name   = NULL,
    drop   = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.20))) +
  coord_flip() +
  labs(
    title    = "Bias: relative prediction error",
    subtitle = paste(
      "Point = rBias, bar = 95% CI (cluster bootstrap).",
      "Shaded band = +/-20% acceptability;",
      "\na model passes when it lies inside the band AND its CI covers zero."
    ),
    x = NULL,
    y = "rBias (%)"
  ) +
  theme_metrics

plot_rRMSE <- ggplot(results_df, aes(x = Model, y = rRMSE)) +
  geom_col(fill = "#4575b4", width = 0.7, alpha = 0.8) +
  geom_point(aes(y = rRMSE)) +
  geom_errorbar(aes(ymin = rRMSE_Lower_CI, ymax = rRMSE_Upper_CI),
                width = 0.25, linewidth = 0.6, colour = "grey20") +
  geom_text(aes(y = rRMSE_Upper_CI, label = sprintf("%.1f", rRMSE)),
            hjust = -0.3, size = 4, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  coord_flip() +
  labs(
    title    = "Precision: relative RMSE",
    subtitle = "Bar = rRMSE, whisker = 95% CI (cluster bootstrap)",
    x = NULL,
    y = "rRMSE (%)"
  ) +
  theme_metrics

# F20 / F30: share of predictions falling within +/-20% and +/-30% of the
# observed concentration. Reshaped to one row per model per band so the two
# bands can be dodged side by side.
f_df <- results_df |>  
  transmute(
    Model,
    F20_est = F20, F20_lo = F20_Lower_CI, F20_hi = F20_Upper_CI,
    F30_est = F30, F30_lo = F30_Lower_CI, F30_hi = F30_Upper_CI
  ) %>% 
  pivot_longer(
    cols      = -Model,
    names_to  = c("Band", ".value"),
    names_sep = "_"
  )

plot_F <- ggplot(f_df, aes(x = Model, y = est, fill = Band)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_point(aes(y = est),  position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.6, colour = "grey20") +
  geom_text(aes(y = hi, label = sprintf("%.0f", est)),
            position = position_dodge(width = 0.75),
            hjust = -0.3, size = 3.5, colour = "grey20") +
  scale_fill_manual(
    values = c(F20 = "#a6bddb", F30 = "#2b8cbe"),
    labels = c(F20 = "within +/-20%", F30 = "within +/-30%"),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(
    title    = "Accuracy: samples predicted within +/-20% and +/-30%",
    subtitle = "Bar = percentage of samples, whisker = 95% CI (cluster bootstrap). Higher is better.",
    x = NULL,
    y = "Samples (%)"
  ) +
  theme_metrics


combined_results_tot <- bind_rows(combined_results, .id = "model")


plot_ipred <- ggplot(combined_results_tot |> filter_out(is.na(DV)), aes(x = DV, y = Cc)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 10) +  # Add labels with ID column
  labs(title = "Model Fit",
       x = "Observed",
       y = "Predicted") +
  theme_classic() +
  coord_fixed()+
  facet_wrap(~model)

obs_df <- df %>%
  filter(EVID %in% c(0), !is.na(DV)) %>%
  transmute(ID, TIME, OBSERVED = DV) %>%
  distinct()

for (i in sort(unique(combined_results_tot$ID))) {
ggplot(combined_results_tot |> filter(ID == i), aes(x = time, y = Cc, group = ID))+
  geom_line() +
  geom_point(data = obs_df |> filter(ID == i), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  theme_minimal() +
  facet_wrap(~ model, scales = "free") 
  
  ggsave(paste0("Cefepime/Concentration/POP/",i,".png"))

}

# MED ---------------------------------------------------------------------
#
# Multivariate Exact Discrepancy (Baklouti et al., Clin Pharmacokinet 2023;
# 62:1599-1609). The method and its caveats are documented at length at the top
# of Utils/med.R; this section only applies it. In one line: for each patient it
# returns the probability that the model produces a kinetic profile at least as
# atypical as the one observed, so
#
#   ACROSS PATIENTS  those probabilities must be Uniform(0,1) under a correct
#                    model -> a Kolmogorov-Smirnov test is a global
#                    goodness-of-fit test for the model on this cohort.
#   PER PATIENT      a very small value marks a whole profile the model cannot
#                    produce, i.e. an outlier kinetic, not an outlier sample.
#
# WHAT THIS ADDS TO THE METRICS ABOVE. rBias/rRMSE score one prediction against
# one observation and average; they never ask whether a patient's samples are
# JOINTLY plausible. A model can land every point inside a VPC interval and
# still be unable to produce the pair of points actually seen, because the slope
# they imply is impossible (the paper's Fig. 3 is exactly that case). MED is
# also the only one of these that comes with a null distribution, so it yields a
# p-value rather than a number to be eyeballed against a convention.
#
# THIS IS A PRIOR-PREDICTIVE CHECK. Nothing here is fitted: the etas are
# integrated out by simulation from the published Omega, so the observations are
# never used twice and there is no shrinkage to calibrate -- unlike the
# Mahalanobis d2 section further down, which asks a similar question in
# parameter space and needs a whole calibration block because of it. The two are
# complementary: d2 says WHICH parameter is unusual, MED says how unusual the
# patient is. Run on the full dataset (all TDMs), not on the masked frames --
# masking exists to keep the MAP fits honest, and MED has no fit to keep honest.

med_K      <- 1e5     # paper uses 1e5. See LL_alpha note below before lowering.
med_scale  <- "log"   # hold FIXED across models: MED is not scale-invariant
med_alpha  <- 0.05
med_allow_iov <- FALSE

# SET THIS TO THE ASSAY'S ACTUAL LLOQ. Simulated concentrations below it are
# censored to it rather than dropped, which is the correct treatment (a
# concentration under LLOQ is not observable) and avoids deleting exactly the
# low-density profiles the tail of MED depends on.
#
# READ max_censored AS A RESULT, NOT AS A TUNING DIAGNOSTIC. On this cohort only
# 2 of 39 patients exceed 1% censoring under Jonckheere_2016, and the worst
# (ID 96, 31%) is censored because the model's MEDIAN simulated concentration at
# the second sampling time is 0.0089 mg/L against 16.8 mg/L actually measured.
# That is the model putting a trough three orders of magnitude too low, not a
# badly chosen floor -- raising med_lloq would only move the pile-up. What the
# censoring DOES cost is precision in the MED value itself: profiles collapsed
# onto the floor share one coordinate, which inflates the simulated spread in
# that direction and makes the patient look less extreme than they are. So a
# heavily censored patient is reliably flagged as atypical but their MED is a
# lower bound on how atypical. censor_heavy below marks them.
med_lloq <- 0.1

# ks gives the paper's kernel-density estimator; without it, fall back to the
# elliptical (Mahalanobis) ranking of the same simulated sample, which agrees
# with MED whenever the log-concentration density has ellipsoidal contours and
# is far quicker. Both are valid uniform-under-the-null statistics.
med_method <-  "kde"# "mahalanobis"

# WHICH MODELS CAN BE SCORED AT ALL, STATED OUT LOUD. Two things disqualify a
# model, and both would otherwise show up only as a facet that never appears --
# the same silent-dropout problem called out in the ETA section below.
#   missing covariate  the cohort does not carry a column the model needs.
#   IOV                poso_simu_pop fixes every KAPPA at 0, so the simulated
#                      spread is too narrow and MED is biased toward 0, i.e.
#                      against the model. med_allow_iov accepts that knowingly.
med_eligibility <- tibble(Model = names(cefepime_models)) |>
  mutate(
    missing_cov = map_chr(cefepime_models, ~ paste(
      setdiff(.x$covariates, names(population_df)), collapse = ", ")),
    has_iov     = map_lgl(cefepime_models, ~ !is.null(.x$pi_matrix)),
    reason      = case_when(
      missing_cov != ""            ~ paste("missing covariate:", missing_cov),
      has_iov & !med_allow_iov     ~ "IOV model: KAPPAs not simulated",
      TRUE                         ~ NA_character_
    ),
    eligible = is.na(reason)
  )

print(med_eligibility |> select(Model, eligible, reason), n = Inf)
if (!any(med_eligibility$eligible)) {
  stop("No model is eligible for MED; see med_eligibility for the reasons.")
}

# Runtime: one solve of K profiles per patient per model. On this cohort that is
# roughly a few minutes per model at K = 1e5. Lower med_K to ~2e4 if only the
# KS test is wanted -- but see LL_alpha below, the outlier rule needs the tail.
med_raw <- list()
for (mn in med_eligibility$Model[med_eligibility$eligible]) {
  message("MED: ", mn)
  med_raw[[mn]] <- med_cohort(population_df, cefepime_models[[mn]],
                              K = med_K, method = med_method,
                              scale = med_scale, floor = med_lloq,
                              allow_iov = med_allow_iov, verbose = FALSE) |>
    med_flag_outliers(alpha = med_alpha)
}

med_tbl <- bind_rows(med_raw, .id = "Model") |>
  mutate(
    Model = factor(Model, levels = model_levels),
    # 1% of simulated concentrations landing under LLOQ is already an order of
    # magnitude more than the ~0.13% a 33% CV proportional error produces on its
    # own, so it marks a patient whose MED is a lower bound (see the med_lloq
    # note above) rather than a point estimate.
    censor_heavy = !is.na(frac_censored) & frac_censored > 0.01
  )

med_summary <- med_tbl |>
  group_by(Model) |>
  summarise(
    n            = sum(!is.na(med)),
    n_failed     = sum(failed),
    mean_MED     = mean(med, na.rm = TRUE),
    KS_D         = med_uniformity(med)$D,
    KS_p         = med_uniformity(med)$p.value,
    n_outlier    = sum(outlier, na.rm = TRUE),
    min_MED      = min(med, na.rm = TRUE),
    max_censored = max(frac_censored, na.rm = TRUE),
    n_censor_heavy = sum(censor_heavy),
    .groups      = "drop"
  ) |>
  # Same shape of rule as bias_ok above: a model whose MEDs are not credibly
  # uniform is invalidated on this cohort. Note the asymmetry -- failing is
  # evidence against the model, but passing with N = 39 is weak evidence for it,
  # because the KS test has little power at this sample size.
  mutate(uniform_ok = KS_p >= med_alpha) |>
  arrange(desc(KS_p))

print(med_summary, n = Inf)

# The per-patient threshold is NOT alpha. The smallest of N uniforms is small by
# construction, so the paper controls the probability that the MINIMUM falls
# below the bound: LL_alpha = 1 - (1 - alpha)^(1/N). It is also the resolution
# med_K has to deliver -- a MED can only be seen to be below LL_alpha if enough
# simulated profiles fit underneath it.
# med_flag_outliers already applied the threshold PER MODEL, using that model's
# own count of scored patients, and those values are in med_tbl$LL_alpha. The
# scalar below is only the reference line for the plot: the most conservative
# (largest) of them, which is also the common value whenever every model scored
# every patient.
med_N  <- n_distinct(med_tbl$ID[!is.na(med_tbl$med)])
med_LL <- max(med_tbl$LL_alpha, na.rm = TRUE)
message(sprintf(
  "MED: N = %d, LL_alpha = %.4g at alpha = %.2f; K = %s resolves down to %.2g",
  med_N, med_LL, med_alpha, format(med_K, scientific = FALSE), 1 / med_K))

# --- is the distribution uniform? ----------------------------------------
# The direct picture of the KS test: the empirical CDF of the MEDs against the
# diagonal it should follow. Bowing ABOVE the diagonal means the MEDs bunch low,
# i.e. the observed profiles are more extreme than the model routinely produces
# (too narrow, or the wrong structure); bowing below means the model is
# over-dispersed relative to the data.
plot_med_ecdf <- ggplot(med_tbl |> filter(!is.na(med)), aes(x = med)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  stat_ecdf(geom = "step", linewidth = 0.8, colour = "#2166ac") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(
    title    = "MED against the uniform it should follow",
    subtitle = sprintf("%s estimator, %s scale, K = %s; dashed = U(0,1)",
                       med_method, med_scale, format(med_K, scientific = FALSE)),
    x = "MED", y = "empirical CDF"
  ) +
  facet_wrap(~ Model) +
  theme_metrics

# --- which patients does each model fail to produce? ---------------------
# log10 because everything interesting happens in the bottom decade: LL_alpha is
# 0.0013 here, which is invisible on a linear axis. Points are floored at 1/K
# because a MED of exactly 0 has no position on a log axis -- it means "no
# simulated profile out of K was less typical", so 1/K is the honest upper bound
# on it.
#
# Labelling only pays off while the outliers are few. A model that fails
# wholesale piles a dozen patients onto the 1/K floor at the same coordinate,
# where repelled labels cannot be placed and say nothing anyway -- the model is
# already dead by the KS test. Above the cut the count is annotated instead.
med_label_max <- 8
med_outlier_n <- med_tbl |> filter(outlier) |> count(Model, name = "n_out")

plot_med_patients <- ggplot(
  med_tbl |> filter(!is.na(med)),
  aes(x = Model, y = pmax(med, 1 / med_K), colour = outlier)
) +
  geom_hline(yintercept = med_LL, linetype = "dashed", colour = "#b2182b") +
  geom_point(position = position_jitter(width = 0.3, height = 0), size = 2,
             alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = ~ .x |> filter(outlier) |> semi_join(
      filter(med_outlier_n, n_out <= med_label_max), by = "Model"),
    aes(label = ID), size = 3, box.padding = 0.4, max.overlaps = Inf,
    show.legend = FALSE) +
  geom_text(
    data = med_outlier_n |> filter(n_out > med_label_max) |>
      mutate(med = 1 / med_K, outlier = TRUE),
    aes(x = Model, y = med, label = paste(n_out, "outliers")),
    hjust = -0.15, vjust = -1.1, size = 3.2, show.legend = FALSE) +
  scale_y_log10() +
  scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "grey35"),
                      name = "outlier kinetic") +
  coord_flip() +
  labs(
    title    = "Per-patient MED",
    subtitle = sprintf(
      "dashed = LL_alpha = %.4g (alpha = %.2f, N = %d); points floored at 1/K",
      med_LL, med_alpha, med_N),
    x = NULL, y = "MED (log scale)"
  ) +
  theme_metrics


# MAP estimation ----------------------------------------------------------
map_data <- df_tdm1   ; map_label <- "1st TDM only"

results2 <- list()   # raw posologyr MAP fits, [[ID]][[model]]
m_out   <- list()    # predicted profiles, flat so bind_rows() still works
m_score <- list()    # one row per observation -- this is what metrics use
m_sir   <- list()    # SIR posterior eta DRAWS, one long tibble per (ID, model)

# --- SIR settings --------------------------------------------------------
# poso_estim_sir costs sir_n_sample ODE solves per fit where poso_estim_map
# costs a handful, so it is on its own switch: sir_run <- FALSE restores the
# MAP loop at its former price. MEASURED, on this data and this machine: one
# fit at n_sample = 10000 takes ~80 s, so the full 39 patients x 9 models is of
# the order of EIGHT HOURS.
#
# DO NOT CUT sir_n_sample TO BUY THAT TIME BACK. Measured on 4 patients x 2
# models at n_sample = 1500: two of the sixteen patient-model-eta rows came
# back with fewer than FIFTY distinct draws out of 1000 -- the importance
# weights had collapsed onto a handful of proposals, and the posterior SD such
# a row reports is an artefact of which few proposals happened to land near the
# likelihood. The cost is not linear in quality: 10000 is what keeps the worst
# patient usable, not what makes the median patient nicer, and the median
# patient is not the one the conclusions rest on. If the run has to be cheaper,
# cut the number of MODELS or patients -- both leave the surviving rows valid.
# The n_unique column in sir_post is what checks this after the fact.
#
# sir_n_resample is the size of the posterior sample that survives, and
# therefore the ceiling on the resolution of every quantile taken from it
# downstream -- 1000 draws puts a useful standard error on a 2.5th percentile,
# 100 would not. It is sampling WITH replacement from the n_sample proposals,
# so raising it past a fraction of n_sample buys duplicate rows, not
# information: the effective sample size is set by the weights, not by this.
sir_run        <- TRUE
sir_n_sample   <- 1e4   # etas proposed from the prior
sir_n_resample <- 1e3    # etas kept, resampled by importance weight

# EVERY true observation, flagged by whether the fit was allowed to see it.
# In a masked dataset the hidden samples keep their DV but carry EVID = 2, so
# posologyr leaves them out of the likelihood while the true value stays here
# for scoring. This flag is what makes a fair comparison possible downstream:
#   used_in_fit = TRUE   the fit saw it -> in-sample, optimistically biased
#   used_in_fit = FALSE  it was hidden  -> out-of-sample, the honest number
obs_df <- map_data |> 
  filter(EVID %in% c(0, 2), !is.na(DV)) |> 
  transmute(ID, TIME, OBSERVED = DV, used_in_fit = EVID == 0) |> 
  distinct()
#sort(unique(map_data$ID))

models_to_use <- setdiff(names(cefepime_models),"Bilal_2023")

for (i in sort(unique(map_data$ID))) {
  patient_data_for_map <- map_data |> 
    filter(ID == i) |> 
    arrange(TIME)
  
  obs_i <- obs_df |> 
    filter(ID == i) |> 
    select(TIME, OBSERVED, used_in_fit) |> 
    distinct()
  
  results2[[as.character(i)]] <- list()
  
  for (model_name in models_to_use) {
    message("Testing model ", model_name, " for ID ", i)
    
    result2_i <- tryCatch({
      model <- cefepime_models[[model_name]]
      
      sim_results <- posologyr::poso_estim_map(patient_data_for_map, model,
                                               return_ofv   = TRUE,
                                               return_model = TRUE,
                                               nocb = FALSE)
      
      # Full predicted profile, for the model-averaging plots.
      profile <- sim_results$model |> 
        as_tibble() |> 
        select(time, Cc, AUC) |> 
        mutate(ID = i)
      
      # One prediction per observation, by INTERPOLATION onto the solver grid.
      # Do NOT left_join on time: observation times need not land on a grid
      # point (ID 13's samples at 107.25, 179.25, 251.58, 353.92 fall between
      # the 0.1 h steps), and such rows silently became OBSERVED = NA and were
      # dropped from every metric -- losing the whole patient without a warning.
      # Duplicate grid times (dose + sample at the same instant) are collapsed to
      # the first, i.e. the pre-dose value, which is what a trough sample means.
      grid <- profile |> 
        distinct(time, .keep_all = TRUE) |> 
        arrange(time)
      
      scored_i <- obs_i |> 
        mutate(
          predicted = approx(grid$time, grid$Cc, xout = TIME, rule = 1)$y,
          ID = i, 
          model_name = model_name, 
          dataset = map_label,
          n_seen = sum(obs_i$used_in_fit)
        ) %>%
        rename(observed = OBSERVED)
      
      if (anyNA(scored_i$predicted))
        message("  ", model_name, " ID ", i, ": ",
                sum(is.na(scored_i$predicted)),
                " osservazione/i fuori dalla griglia temporale")
      
      # SSE describes the FIT, so it uses only the points the fit actually saw.
      # Including hidden points would make the model-averaging weights depend on
      # data the fit never used.
      sse <- scored_i |> 
        filter(!is.na(predicted), used_in_fit) |> 
        summarise(sse = sum((predicted - observed)^2)) |> 
        pull(sse)
      
      # Computing the AIC based on the KAPPA matrix and number of occasions 
      n_occ_i <- dplyr::n_distinct(patient_data_for_map$OCC)
      k_i     <- nrow(model$omega) +
        if (is.null(model$pi_matrix)) 0L else nrow(model$pi_matrix) * n_occ_i
      
      # --- Poso estim SIR: the whole posterior, not just its mode ----------
      sir_i <- NULL
      if (sir_run) {
        sir_i <- tryCatch({
         

          sir_fit <- posologyr::poso_estim_sir(patient_data_for_map, model,
                                               n_sample     = sir_n_sample,
                                               n_resample   = sir_n_resample,
                                               return_model = FALSE,
                                               nocb = FALSE)

          # poso_estim_sir$eta is a data.frame, one row per resampled draw and
          # one column per entry of the model's FULL Omega -- including entries
          # the fit does not estimate, which come back as structural zeros, not
          # as estimates. So restrict to the etas the MAP fit reports, BY NAME,
          # for the reason .omega_of gives below: name-indexing guarantees a
          # draw is never divided by another parameter's prior SD, and errors
          # rather than doing it quietly.
          eta_names_i <- names(sim_results$eta)
          stopifnot(all(eta_names_i %in% names(sir_fit$eta)))

          om_i    <- model$omega[eta_names_i, eta_names_i, drop = FALSE]
          om_sd_i <- setNames(sqrt(diag(om_i)), eta_names_i)

          # Long, and carrying the same columns as eta_tbl below (eta_name,
          # eta, omega_sd, eta_z, n_seen) so the posterior sample and the MAP
          # point estimate can be put on one plot without a reshape. eta_map
          # rides along per draw: it is constant within (ID, model) and is what
          # the draws get compared against.
          sir_fit$eta[, eta_names_i, drop = FALSE] |>
            as_tibble() |>
            mutate(draw = row_number(), .before = 1) |>
            pivot_longer(all_of(eta_names_i),
                         names_to = "eta_name", values_to = "eta") |>
            mutate(
              omega_sd   = unname(om_sd_i[eta_name]),
              eta_z      = eta / omega_sd,
              eta_map    = as.numeric(sim_results$eta[eta_name]),
              ID         = i,
              model_name = model_name,
              dataset    = map_label,
              n_seen     = sum(obs_i$used_in_fit)
            )

        }, error = function(e) {
          message("  SIR failed for ", model_name, " ID ", i, ": ",
                  conditionMessage(e))
          NULL
        })
      }

      list(
        fit    = sim_results,
        sir    = sir_i,
        scored = scored_i,
        tidy   = profile |> 
          mutate(
            ofv        = sim_results$ofv,
            LL         = exp(-0.5 * ofv),
            AIC        = exp(-0.5 * ofv - k_i),
            SSE        = exp(-0.5 * sse),   # squared prediction errors
            # log-scale versions of the three above: exp() underflows to 0 for
            # OFV/SSE above ~1420, which would make every weight NaN
            ll_exp     = -0.5 * ofv,
            aic_exp    = -0.5 * ofv - k_i,
            sse_exp    = -0.5 * sse,
            model_name = model_name,
            dataset    = map_label,          # which design produced this row
            n_seen     = sum(obs_i$used_in_fit)   # TDM available to this fit
          )
      )
      
    }, error = function(e) {
      message("  Model ", model_name, " failed for ID ", i, ": ", conditionMessage(e))
      NULL
    })
    
    results2[[as.character(i)]][[model_name]] <- result2_i$fit
    m_out[[paste(i, model_name, sep = "_")]]   <- result2_i$tidy
    m_score[[paste(i, model_name, sep = "_")]] <- result2_i$scored
    m_sir[[paste(i, model_name, sep = "_")]]   <- result2_i$sir
  }
}

# One long table of posterior draws: ID x model_name x eta_name x draw. Empty
# tibble, not an error, when sir_run was FALSE.
sir_eta_tbl <- purrr::list_rbind(purrr::compact(m_sir))

# ETAs distribution -------------------------------------------------------

# WHAT THIS SECTION ASKS
# Every model ships a prior for its random effects: eta_i ~ MVN(0, Omega),
# Omega being the covariance matrix estimated in the original publication. MAP
# estimation then produces one eta vector per patient per model. Plotting the
# estimates against the prior they are supposed to be drawn from asks one
# question: does this model's picture of the population describe THIS cohort?
#
# eta_z = eta / sqrt(Omega_jj) rescales an eta by the prior SD for that
# parameter, so under the model it is a draw from N(0, 1). |eta_z| > 1.96 puts
# the patient outside the central 95% of the population the model was built on,
# for that one parameter. That is the colouring in the figure below.
#
# Only the ETAs are covered. posologyr also estimates the KAPPAs of an IOV
# model, but it does not return them in fit$eta -- they come back as columns of
# fit$model -- so they sit outside this comparison entirely.
#
# THE SAMPLING DESIGN SITS UNDER EVERY PANEL
# With map_data = df_tdm1 each fit saw exactly ONE observation per patient. One
# sample against 2-4 etas means at most one direction in eta space is informed
# by data and the rest is prior, which is why so many panels collapse onto zero.
# plot_eta_sir below is where that is visible rather than inferred: a posterior
# no taller than its prior is a parameter the sample never spoke about. Change
# map_data and this paragraph stops being true, so check it when you do.


# --- helper: the prior covariance for one fit ----------------------------
# Takes a model name and the eta vector that fit returned, and gives back that
# model's Omega cut down to those etas, in that order. Indexing by NAME rather
# than position matters: lotri labels Omega from the model definition
# (ETA_CL, ETA_V1, ...) and fit$eta carries the same labels, so name-indexing
# guarantees each eta is paired with its own variance. Were the two orders ever
# to disagree, positional indexing would silently pair ETA_CL with the variance
# of ETA_V1; name-indexing errors instead.
.omega_of <- function(model_name, eta) {
  # The one way name-indexing can still fail QUIETLY: an unnamed eta makes
  # om[NULL, NULL] a 0x0 matrix, which yields an empty tibble rather than an
  # error, and the patient disappears from the table unnoticed.
  stopifnot(!is.null(names(eta)), all(nzchar(names(eta))), !anyNA(names(eta)))
  om <- cefepime_models[[model_name]]$omega   # full prior covariance matrix
  om[names(eta), names(eta), drop = FALSE]    # rows and cols aligned to eta
}


# --- one row per (patient, model, eta) -----------------------------------
# results2 is a two-level named list: results2[[patient]][[model]] holds the
# posologyr fit. imap() walks a named list handing the function the element AND
# its name -- which is where `id` and `model_name` come from, without having to
# carry them separately -- and list_rbind() stacks the results. (imap_dfr() did
# this in one step but is superseded in purrr >= 1.0.)
eta_tbl <- purrr::imap(results2, function(models_i, id) {
  # compact() drops the NULL entries left behind by models that errored for
  # this patient, so a failed fit simply never reaches the table.
  purrr::imap(purrr::compact(models_i), function(fit, model_name) {
    eta <- fit$eta                      # named vector, one entry per random effect
    om  <- .omega_of(model_name, eta)   # its prior covariance, aligned to it
    tibble(
      ID         = as.integer(id),
      model_name = model_name,
      eta_name   = names(eta),          # ETA_CL, ETA_V1, ...
      eta        = as.numeric(eta),     # the MAP estimate itself
      # diag(om) are the prior VARIANCES; their square roots are the SDs the
      # publication claims for the population.
      omega_sd   = sqrt(diag(om)),
      # Dividing by that SD puts every eta on one common scale, so a value from
      # ETA_CL in one model is comparable with ETA_V1 in another.
      eta_z      = as.numeric(eta) / sqrt(diag(om))
    )
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()


# --- the prior curve, with one vertical line per patient -----------------
# Panel-level constants: the cohort offset the label reports, the axis
# half-width, and the height of the curve at its peak.
eta_panels <- eta_tbl |>
  group_by(model_name, eta_name) |>
  summarise(
    omega_sd = first(omega_sd),
    mean_eta = mean(eta),               # cohort offset, the label's subject
    # Axis half-width. Deliberately not a plain +/-4 SD: this cohort reaches
    # |eta_z| = 7.4 (Seo_2023 ETA_V2), and a fixed axis would crop that patient
    # off the panel, hiding the point precisely because it is extreme. Taking
    # the wider of the two keeps every patient on screen.
    lim      = max(4 * omega_sd, 1.15 * max(abs(eta))),
    # Height of a N(0, sd) density at its mode: the top of the vertical segments.
    peak     = dnorm(0, 0, omega_sd),
    .groups  = "drop"
  )

# The prior curves: for each panel, 400 x values across its own axis and the
# normal density there. map() builds one vector of x per panel, unnest()
# expands those into rows, giving a long frame ggplot can draw directly.
prior_curves <- eta_panels |>
  mutate(x = purrr::map(lim, ~ seq(-.x, .x, length.out = 400))) |>
  tidyr::unnest(x) |>
  mutate(density = dnorm(x, 0, omega_sd))

# The patient lines: every eta, plus its panel's peak height to draw up to, and
# a flag for whether it falls outside the prior's central 95%.
eta_lines <- eta_tbl |>
  left_join(eta_panels |> select(model_name, eta_name, peak),
            by = c("model_name", "eta_name")) |>
  mutate(outside = abs(eta_z) > 1.96)

# The cohort offset, printed in each panel so the number sits next to the shape
# it describes. Used by plot_eta_prior only (plot_eta_sir is per patient, where a
# cohort mean says nothing), kept as its own table for readability.
eta_labs <- eta_panels |>
  mutate(lab = sprintf("mean eta %+.2f", mean_eta))

plot_eta_prior <- ggplot(prior_curves, aes(x, density)) +
  geom_area(fill = "#4575b4", alpha = 0.13) +          # shaded prior
  geom_line(colour = "#4575b4", linewidth = 0.8) +     # its outline
  geom_segment(
    data = eta_lines,
    # y = 0 up to yend = peak, i.e. the FULL height of the curve rather than
    # the density at that eta. Drawing to dnorm(eta) looks neater, but a tail
    # eta has density near zero, so its line would shrink to nothing exactly
    # when it is the one worth seeing.
    aes(x = eta, xend = eta, y = 0, yend = peak, colour = outside),
    linewidth = 0.55, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +  # prior mean
  # geom_label rather than geom_text: the top-left corner is not reliably empty
  # (the curve often peaks there when the cohort is shifted), and a white box
  # keeps the number readable wherever it lands.
  geom_label(
    data = eta_labs, aes(x = -lim, y = peak, label = lab),
    hjust = 0, vjust = 1, size = 3.1, colour = "grey20",
    fill = "white", alpha = 0.75, linewidth = 0,
    lineheight = 0.95, inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#1b7837", `TRUE` = "#b2182b"),
    labels = c(`FALSE` = "within the prior's central 95%", `TRUE` = "outside it"),
    # drop = FALSE keeps both legend keys even when every patient falls on one
    # side, so the panels stay comparable with each other.
    name = NULL, drop = FALSE
  ) +
  # scales = "free" because Omega differs by orders of magnitude between
  # parameters: a shared axis would flatten the tight ones into a spike.
  # Three breaks, not the default six: free y scales in a short panel put the
  # default labels on top of each other. The y axis carries nothing here anyway
  # (every segment is drawn to the same peak), so fewer is strictly better.
  scale_y_continuous(n.breaks = 3) +
  facet_wrap(~ model_name + eta_name, scales = "free") +
  labs(
    title    = "Estimated etas against the distribution each model assumes",
    subtitle = "Curve = the N(0, Omega_jj) declared by the model. Each vertical line is one patient's MAP eta.",
    x = "eta", y = "density"
  ) +
  theme_metrics +
  theme(legend.position = "top")


# --- the same panels, with each patient's POSTERIOR instead of its mode ---
# plot_eta_prior draws one vertical line per patient: the MAP eta, which is the
# MODE of p(eta | data). A mode carries no width, so that figure cannot tell
# "the data put this patient here" apart from "the data said nothing and MAP
# left the patient at the population value". Under a one-sample design that is
# the whole question. This plot replaces each line with the posterior it is the
# mode of, taken from the SIR draws.
#
# READ THE HEIGHT, NOT ONLY THE POSITION
# Prior and posteriors are both TRUE densities on one axis, so within a panel
# their heights are directly comparable:
#   posterior as broad and as low as the prior  the sample carried no
#     information about this parameter. The MAP "individual" estimate is the
#     population value under another name, and the matching line in
#     plot_eta_prior means far less than it looks.
#   posterior clearly taller and narrower       the sample was informative
#     here, and the MAP line is worth reading as an individual estimate.
# Same comparison the shrinkage table used to make, but per patient rather than
# per cohort, and needing no calibration: it is read straight off the picture.
#
# Colour uses the SAME rule as plot_eta_prior -- the MAP eta outside the prior's
# central 95% -- rather than anything computed from the draws, so a patient
# keeps one colour across the two figures and they can be read side by side.
#
# Panels where SIR did not run show the prior curve alone. That is deliberate:
# an absent posterior should look absent rather than like a missing patient.
#
# A DEGENERATE POSTERIOR IS NARROW FOR THE WRONG REASON, so it is drawn dashed.
# If the importance weights collapse onto a handful of proposals, the surviving
# draws are few and close together and the density comes out TALL AND NARROW --
# which this plot's whole reading would score as "the sample was informative".
# It is the opposite: it means SIR could not resolve the posterior. Measured on
# a 12-patient trial sample, the degenerate rows had a median peak ratio of 1.87
# against 1.15 for the healthy ones, and owned the largest ratio in the set, so
# this is the dominant confound here rather than a corner case. n_unique is
# recomputed here rather than taken from sir_post because that block runs after
# this one; .SIR_MIN_UNIQUE is shared so the two cannot drift apart.
.SIR_MIN_UNIQUE <- 50

# ONE FIGURE PER PATIENT, saved to Cefepime/ETA/SIR/<ID>.png, same pattern as
# the Concentration/MAP loop. Overlaying all 39 posteriors in one panel made
# every individual curve unreadable, and a single patient's curve is the unit
# the reading above is about. The prior curves are the same prior_curves for
# every patient, so heights stay comparable ACROSS files as well as within one.
#
# In each panel: the shaded prior, this patient's SIR posterior (dashed when
# its weights degenerated), and a vertical line at this patient's MAP eta. The
# label gives the MAP value and sd_ratio = posterior SD / prior SD -- the same
# quantity sir_post reports -- so the number sits next to the shape.
#
# The cohort label eta_labs is NOT used here: "mean eta" over 39 patients says
# nothing about the one patient in the file.
#
# plot_eta_sir becomes a named list, one ggplot per patient, so any one of them
# can still be looked at interactively: plot_eta_sir[["13"]].
if (nrow(sir_eta_tbl) == 0) {
  message("sir_eta_tbl is empty (sir_run was FALSE): plot_eta_sir not built.")
  plot_eta_sir <- NULL
} else {
  sir_dens <- sir_eta_tbl |>
    group_by(ID, model_name, eta_name) |>
    mutate(weights_ok = dplyr::n_distinct(eta) >= .SIR_MIN_UNIQUE) |>
    ungroup() |>
    mutate(outside = abs(eta_map / omega_sd) > 1.96)

  # One row per panel per patient: the MAP line and the label.
  sir_marks <- sir_dens |>
    group_by(ID, model_name, eta_name) |>
    summarise(
      eta_map    = first(eta_map),
      outside    = first(outside),
      weights_ok = first(weights_ok),
      sd_ratio   = sd(eta) / first(omega_sd),
      .groups    = "drop"
    ) |>
    mutate(lab = sprintf("MAP %+.2f\nSD ratio %.2f%s", eta_map, sd_ratio,
                         if_else(weights_ok, "", "\n(weights degenerated)")))

  out_dir <- "Cefepime/ETA/SIR"
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  plot_eta_sir <- list()

  for (i in sort(unique(sir_dens$ID))) {
    dens_i  <- sir_dens  |> filter(ID == i)
    marks_i <- sir_marks |> filter(ID == i)

    # Only the panels this patient has a posterior for. A model that failed
    # for this patient would otherwise leave a prior-only panel that reads as
    # "no information" when it actually means "no fit".
    prior_i <- prior_curves |>
      semi_join(marks_i, by = c("model_name", "eta_name"))

    p <- ggplot(prior_i, aes(x, density)) +
      geom_area(fill = "#4575b4", alpha = 0.13) +          # shaded prior
      geom_line(colour = "#4575b4", linewidth = 0.8) +     # its outline
      # y defaults to after_stat(density), the same units dnorm() put the prior
      # curve in -- that is what makes the two heights comparable. Do not switch
      # this to ..scaled.. or the informativeness reading silently stops being
      # true.
      geom_density(
        data = dens_i,
        aes(x = eta, colour = outside, linetype = weights_ok),
        linewidth = 0.7, fill = NA, inherit.aes = FALSE
      ) +
      geom_vline(data = marks_i, aes(xintercept = eta_map, colour = outside),
                 linewidth = 0.5, show.legend = FALSE) +       # this patient's MAP
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +  # prior mean
      geom_label(
        data = marks_i, aes(x = -Inf, y = Inf, label = lab),
        hjust = -0.05, vjust = 1.1, size = 3, colour = "grey20",
        fill = "white", alpha = 0.75, linewidth = 0,
        lineheight = 0.95, inherit.aes = FALSE
      ) +
      scale_linetype_manual(
        values = c(`TRUE` = "solid", `FALSE` = "22"),
        labels = c(`TRUE`  = "SIR weights healthy",
                   `FALSE` = "weights degenerated: width is noise"),
        name = NULL, drop = FALSE
      ) +
      scale_colour_manual(
        values = c(`FALSE` = "#1b7837", `TRUE` = "#b2182b"),
        labels = c(`FALSE` = "MAP eta within the prior's central 95%",
                   `TRUE`  = "MAP eta outside it"),
        name = NULL, drop = FALSE
      ) +
      # Three breaks: free y scales in short panels overlap the default labels,
      # and here the height IS the message.
      scale_y_continuous(n.breaks = 3) +
      facet_wrap(~ model_name + eta_name, scales = "free", ncol = 4) +
      labs(
        title    = sprintf("Patient %s: posterior eta against the prior it was drawn from", i),
        subtitle = sprintf(paste(
          "Filled curve = N(0, Omega_jj) declared by the model. Line = this patient's SIR posterior (%d draws); vertical line = MAP eta.",
          "\nA posterior no taller than the prior means the sample (%d observation(s)) carried no information about that parameter."),
          dplyr::n_distinct(dens_i$draw), first(dens_i$n_seen)),
        x = "eta", y = "density"
      ) +
      theme_metrics +
      theme(legend.position = "top", legend.box = "vertical")

    plot_eta_sir[[as.character(i)]] <- p

    # Height scales with the number of panel rows so the panels keep one size
    # whether a patient has 4 fits or 20.
    n_rows <- ceiling(nrow(marks_i) / 4)
    ggsave(file.path(out_dir, paste0(i, ".png")), p,
           width = 14, height = 2.6 * n_rows + 2, dpi = 150, bg = "white")
  }
}


# --- the same reading, as a number ---------------------------------------
# plot_eta_sir asks the reader to compare two curve HEIGHTS by eye. That works
# for a panel or two and stops working at 9 models x 39 patients, so the same
# comparison is computed here.
#
# THE ONE QUANTITY THIS TABLE IS FOR
#   sd_ratio = post_sd / omega_sd
# the width of this patient's posterior as a fraction of the prior it started
# from. 1 means the observation moved nothing: the posterior IS the prior, and
# this patient's MAP eta for this parameter is the population value wearing an
# individual label. 0 would mean the data pinned the parameter exactly. On a
# 1-sample design expect values near 1 for every parameter except the one the
# sample speaks to -- that is the finding, not a failure.
#
# THIS IS NOT ETA-SHRINKAGE, AND MUST NOT BE REPORTED AS IT
# The familiar diagnostic, 1 - sd(eta_MAP over patients)/omega, is a COHORT
# statistic: it compares the spread of the point estimates ACROSS patients with
# the prior SD. sd_ratio here is a WITHIN-patient quantity -- the width of one
# patient's own posterior -- and it needs no cohort to be defined. They answer
# different questions and can disagree: a cohort of patients all genuinely near
# the population mean gives high classical shrinkage with informative
# individual posteriors. Keep the names apart in anything written up.
#
# sd_ratio > 1 is possible and is not a bug. The posterior is wider than the
# prior when the observation conflicts with the covariate-predicted typical
# value hard enough to push probability into both tails, and it also appears as
# pure noise when n_eff below is small. Read the two columns together.
#
# n_unique IS THE CREDIBILITY OF EVERY OTHER COLUMN IN THE ROW
# SIR resamples n_resample draws WITH REPLACEMENT from n_sample proposals, so
# the row count is fixed at n_resample no matter how badly the importance
# weights behaved. When one proposal dominates, the same eta is drawn hundreds
# of times: 1000 rows, a handful of distinct values, and a post_sd that is an
# artefact of which few proposals happened to land near the likelihood rather
# than a description of the posterior. Counting distinct values costs nothing
# and is the only warning available, the weights themselves not being returned.
# Draws are resampled as whole vectors, so this count is per patient-model, not
# per eta -- it repeats across a model's etas by construction.
if (nrow(sir_eta_tbl) == 0) {
  message("sir_eta_tbl is empty (sir_run was FALSE): SIR posterior summaries skipped.")
  sir_post <- NULL
  sir_info <- NULL
} else {
  sir_post <- sir_eta_tbl |>
    group_by(ID, model_name, eta_name) |>
    summarise(
      n_draw   = n(),
      n_unique = dplyr::n_distinct(eta),   # surviving proposals; see above
      n_seen   = first(n_seen),
      omega_sd = first(omega_sd),
      eta_map  = first(eta_map),           # constant within the group
      post_mean = mean(eta),
      post_sd   = sd(eta),
      # Quantiles of the DRAWS, i.e. a credible interval, not +/-1.96 SD: the
      # posterior under one observation is routinely skewed and a symmetric
      # interval would misplace both ends.
      # unname(): quantile() returns a "2.5%"-named scalar, and dplyr keeps that
      # name on the column, which then travels into every join and plot label.
      q025 = unname(quantile(eta, 0.025)),
      q500 = unname(quantile(eta, 0.500)),
      q975 = unname(quantile(eta, 0.975)),
      .groups = "drop"
    ) |>
    mutate(
      sd_ratio  = post_sd / omega_sd,
      # Does the mode sit inside its own posterior's central 95%? Normally yes,
      # and a FALSE is worth looking at: it means the MAP optimiser and the SIR
      # sample disagree about where the mass is, which happens when the
      # posterior is multimodal or when n_unique has collapsed.
      map_in_ci = eta_map >= q025 & eta_map <= q975
    )

  # Per parameter across the cohort. The median is the right centre here: with
  # 39 patients a single patient whose SIR degenerated would drag a mean.
  sir_info <- sir_post |>
    group_by(model_name, eta_name) |>
    summarise(
      n_pat           = dplyr::n_distinct(ID),
      median_sd_ratio = median(sd_ratio),
      iqr_sd_ratio    = IQR(sd_ratio),
      min_n_unique    = min(n_unique),
      n_map_outside   = sum(!map_in_ci),
      .groups = "drop"
    ) |>
    arrange(median_sd_ratio)

  # Printed rather than left in the environment: a parameter whose posterior
  # never narrows is the single thing that invalidates reading its MAP eta as
  # individual, and it should not depend on someone thinking to inspect a table.
  uninformed <- sir_info |> filter(median_sd_ratio > 0.9)
  if (nrow(uninformed) > 0) {
    cat(sprintf(
      "\nSIR: %d parameter(s) keep >90%% of their prior width -- their MAP etas are prior, not data:\n",
      nrow(uninformed)))
    print(uninformed |> select(model_name, eta_name, median_sd_ratio, min_n_unique))
  }

  # Same idea, other direction: rows where the sample itself is too thin to
  # trust whatever the row says.
  # .SIR_MIN_UNIQUE, not a second literal 50: plot_eta_sir dashes exactly these
  # rows, and a threshold defined twice is a threshold that will disagree.
  degenerate <- sir_post |> filter(n_unique < .SIR_MIN_UNIQUE)
  if (nrow(degenerate) > 0)
    warning(sprintf(
      "SIR importance weights degenerated for %d patient-model-eta row(s) (<%d distinct draws of %d). Their post_sd is noise and plot_eta_sir dashes them; raise sir_n_sample.",
      nrow(degenerate), .SIR_MIN_UNIQUE, unique(sir_post$n_draw)[1]), call. = FALSE)
}

# Metrics -----------------------------------------------------------------

m_out_tot <- bind_rows(m_out)

scoring_rows <- bind_rows(m_score) %>%
  select(dataset, model_name, ID, TIME, n_seen, used_in_fit, observed, predicted)

#scoring_rows %>% count(dataset, model_name, name = "n_scored") %>% print(n = Inf)
metrics_scope <- if (any(!scoring_rows$used_in_fit)) "out-of-sample" else "in-sample"
metrics_input <- scoring_rows %>%
  filter(if (metrics_scope == "out-of-sample") !used_in_fit else used_in_fit,
         is.finite(observed), is.finite(predicted), observed > 0)

.ci <- function(m, field, side)
  if (is.null(m[[field]])) NA_real_ else unname(m[[field]][[side]])

metrics_df <- metrics_input %>%
  split(metrics_input$model_name) %>%
  imap_dfr(function(d, nm) {
    m <- calculate_metrics(d, n_bootstrap = 10000, conf_level = 0.95, id_col = "ID")
    tibble(
      Model    = nm,
      n_scored = nrow(d),
      n_pts    = n_distinct(d$ID),
      rBias    = m$rBias, rBias_lo = .ci(m, "rBias_CI", "lower"),
      rBias_hi = .ci(m, "rBias_CI", "upper"),
      rRMSE    = m$rRMSE, rRMSE_lo = .ci(m, "rRMSE_CI", "lower"),
      rRMSE_hi = .ci(m, "rRMSE_CI", "upper"),
      MDPE     = m$MDPE, MDAPE = m$MDAPE,
      F20      = m$F20, F20_lo = .ci(m, "F20_CI", "lower"),
      F20_hi   = .ci(m, "F20_CI", "upper"),
      F30      = m$F30, F30_lo = .ci(m, "F30_CI", "lower"),
      F30_hi   = .ci(m, "F30_CI", "upper")
    )
  }) %>%
  mutate(dataset = paste(unique(scoring_rows$dataset), collapse = ", "),
         scope   = metrics_scope,
         across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(MDAPE)

print(metrics_df, n = Inf)


# MAP plots ---------------------------------------------------------------

# Same three panels as the population section, on the MAP fits. Reuses
# theme_metrics and pass_cols defined there, and reuses model_levels so a model
# occupies the same row in both sets of figures - that is what makes the
# a-priori-versus-MAP comparison readable at a glance. Models absent from one
# section are appended at the end rather than dropped to NA.
map_levels <- if (exists("model_levels")) {
  c(intersect(model_levels, metrics_df$Model),
    setdiff(metrics_df$Model, model_levels))
} else {
  metrics_df %>% arrange(desc(rRMSE)) %>% pull(Model)
}

metrics_plot_df <- metrics_df %>% 
  mutate(
    Model = factor(Model, levels = map_levels),
    bias_ok = abs(rBias) <= 20 & rBias_lo <= 0 & rBias_hi >= 0
  )

# Which design and which scope produced these numbers. Worth carrying on the
# figures themselves: in-sample MAP metrics are optimistically biased because
# the fit already saw the sample it is being scored on, and a figure that does
# not say so is easy to misread later.
map_caption <- sprintf(
  "%s | %d samples from %d patients",
  map_label, nrow(metrics_input), n_distinct(metrics_input$ID)
)

plot_rBias_map <- ggplot(metrics_plot_df, aes(x = Model, y = rBias)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -20, ymax = 20,
           fill = "grey60", alpha = 0.18) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey20") +
  geom_errorbar(aes(ymin = rBias_lo, ymax = rBias_hi, colour = bias_ok),
                width = 0.25, linewidth = 0.7) +
  geom_point(aes(colour = bias_ok), size = 3.2) +
  geom_text(aes(y = rBias_hi, label = sprintf("%+.1f", rBias)),
            hjust = -0.4, size = 4, colour = "grey20") +
  scale_colour_manual(
    values = pass_cols,
    breaks = c("TRUE", "FALSE"),
    labels = c("meets criterion", "fails criterion"),
    name   = NULL,
    drop   = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.20))) +
  coord_flip() +
  labs(
    title    = "MAP bias: relative prediction error",
    subtitle = paste0(
      map_caption,
      "\nPoint = rBias, bar = 95% CI (cluster bootstrap).",
      " Shaded band = +/-20% acceptability;",
      "\na model passes when it lies inside the band AND its CI covers zero."
    ),
    x = NULL,
    y = "rBias (%)"
  ) +
  theme_metrics

plot_rRMSE_map <- ggplot(metrics_plot_df, aes(x = Model, y = rRMSE)) +
  geom_col(fill = "#4575b4", width = 0.7, alpha = 0.8) +
  geom_point() +
  geom_errorbar(aes(ymin = rRMSE_lo, ymax = rRMSE_hi),
                width = 0.25, linewidth = 0.6, colour = "grey20") +
  geom_text(aes(y = rRMSE_hi, label = sprintf("%.1f", rRMSE)),
            hjust = -0.3, size = 4, colour = "grey20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  coord_flip() +
  labs(
    title    = "MAP precision: relative RMSE",
    subtitle = paste0(map_caption,
                      "\nBar = rRMSE, whisker = 95% CI (cluster bootstrap)."),
    x = NULL,
    y = "rRMSE (%)"
  ) +
  theme_metrics

f_df_map <- metrics_plot_df %>% 
  transmute(
    Model,
    F20_est = F20, F20_lo = F20_lo, F20_hi = F20_hi,
    F30_est = F30, F30_lo = F30_lo, F30_hi = F30_hi
  ) %>% 
  pivot_longer(
    cols      = -Model,
    names_to  = c("Band", ".value"),
    names_sep = "_"
  )

plot_F_map <- ggplot(f_df_map, aes(x = Model, y = est, fill = Band)) +
  geom_point(position = position_dodge(width = 0.75))+
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = lo, ymax = hi),
                position = position_dodge(width = 0.75),
                width = 0.25, linewidth = 0.6, colour = "grey20") +
  geom_text(aes(y = hi, label = sprintf("%.0f", est)),
            position = position_dodge(width = 0.75),
            hjust = -0.3, size = 3.5, colour = "grey20") +
  scale_fill_manual(
    values = c(F20 = "#a6bddb", F30 = "#2b8cbe"),
    labels = c(F20 = "within +/-20%", F30 = "within +/-30%"),
    name   = NULL
  ) +
  scale_y_continuous(limits = c(0, 100),
                     expand = expansion(mult = c(0, 0.12))) +
  coord_flip() +
  labs(
    title    = "MAP accuracy: samples predicted within +/-20% and +/-30%",
    subtitle = paste0(map_caption,
                      "\nBar = percentage of samples, whisker = 95% CI (cluster bootstrap). Higher is better."),
    x = NULL,
    y = "Samples (%)"
  ) +
  theme_metrics


# Model averaging ---------------------------------------------------------
# Weights are w_m proportional to exp(-0.5 * OFV_m) (likelihood weights; the
# AIC variant penalises by the number of random effects, the SSE variant uses
# the fit error instead). Two things have to happen and neither did before:
#
#   1. NORMALISE. sum(Cc * w) without dividing by sum(w) is a weighted SUM, not
#      an average, so the result had no meaningful scale.
#   2. SHIFT BEFORE EXPONENTIATING. exp(-0.5 * OFV) with OFV in the hundreds is
#      denormal, and exactly 0 once OFV passes ~1420, so every weight collapsed
#      to 0 and MA_* came out 0. Subtracting the maximum log-weight first is
#      algebraically identical after normalisation (the constant cancels) and
#      keeps the largest weight at exp(0) = 1.
#
# Weights are formed within (ID, time), which holds one row per model, so each
# patient's models are compared against each other and against nothing else.
# Models that failed for a patient are simply absent and the remaining weights
# renormalise over what is there.
ma_weights <- function(log_w) {
  z <- log_w - max(log_w, na.rm = TRUE)
  w <- exp(z)
  w / sum(w, na.rm = TRUE)
}

df_final <- m_out_tot %>%
  group_by(ID, time) %>%
  mutate(
    w_LL   = ma_weights(ll_exp),
    w_AIC  = ma_weights(aic_exp),
    w_SSE  = ma_weights(sse_exp),
    MA_LL  = sum(Cc * w_LL,  na.rm = TRUE),
    MA_AIC = sum(Cc * w_AIC, na.rm = TRUE),
    MA_SSE = sum(Cc * w_SSE, na.rm = TRUE)
  ) %>%
  ungroup()

# m_out_tot holds predicted profiles only; the observations now live in
# scoring_rows, and the plot below draws them from obs_df.
df_final_MA <- df_final |>
  select(ID, time, starts_with("MA_")) |>
  distinct() |>
  pivot_longer(cols = -c(ID, time)) |>
  rename(model_name = name, Cc = value)

df_final <- bind_rows(df_final,df_final_MA )

# group = ID everywhere: df_final now holds every patient, so without it the
# lines connect across patients. Facet on model_name + ID (or filter to one ID)
# if you want one panel per patient rather than all profiles overlaid.

for (i in sort(unique(df_final$ID))) {
ggplot(df_final |> filter(ID == i), aes(x = time, y = Cc, group = ID)) +
  geom_line() +
  geom_point(data = obs_df |> filter(ID == i), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  facet_wrap(~ model_name, scales = "free") 
  
  ggsave(paste0("Cefepime/Concentration/MAP/",i,".png"))
}
