library(tidyverse);library(readxl);library(data.table);library(ggrepel)

source("~/PhD-thesis/Utils/Functions.R")
source("~/PhD-thesis/Utils/med.R")
source("~/PhD-thesis/Ceftazidime/ceftazidime_models.R")
set.seed(456)


# data --------------------------------------------------------------------

total_covariates <- lapply(models_ceftazidime, function(x) x$covariates) %>% unlist() %>% unique()
df <- read_excel("Ceftazidime/data/nantes.xlsx") |> 
  filter_out(is.na(ID)) |> 
  mutate(across(where(is.character), ~na_if(.x, "."))) |> 
  readr::type_convert() |> 
  mutate(
    DV = Y,
    DUR = INFUSION_DURATION,
    DIAL = DIALYSE,
    CVVH = 0,
    TRAUMA = 0, 
    HEMAT = 0, 
    ADM = 3, 
    VM = VENTILATION,
    gender = if_else(SEX == 0, "male", "female", missing = NA_character_),
    CLCR_CG = CG,
    SEX = if_else(SEX == 0, 1, 0), # sex == 0 for male in original dataseta (see Conil model)
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
  filter_out(.tdm_i >= 3)
# Population evaluation ---------------------------------------------------
population_df <- df
pop_n_simul   <- 1000       # 0 = typical profile only; > 0 also simulates
pop_pred_stat <- "typical"  # "typical" | "median" | "mean" -> which becomes Cc
pop_grid_n    <- 400        # cap on plotting-grid points when simulating

results <- list()
m_out <- list()

for (i in sort(unique(population_df$ID))) {
  message("Processing item ", i, " of ", length(unique(population_df$ID)))
  
  df_i <- population_df |> filter(ID == i)
  max_time_i <- max(df_i$TIME, na.rm = TRUE)
  
  by_i <- if (pop_n_simul > 0) max(0.1, max_time_i / pop_grid_n) else 0.1
  et_i <- rxode2::as.et(df_i)
  et_i$add.sampling(seq(0, max_time_i, by = by_i))
  
  obs_i <- df_i |> filter(EVID == 0, !is.na(DV)) |> select(TIME, DV)
  
  id_key <- as.character(i)
  results[[id_key]] <- list()
  
  for (model_name in names(models_ceftazidime)) {
    message("Testing model ", model_name, " for ID ", i)
    
    result_i <- tryCatch({
      model <- models_ceftazidime[[model_name]]
      
      solve_pop <- function(K) {
        s <- posologyr::poso_simu_pop(df_i, model, n_simul = K)
        s$model <- posologyr::poso_replace_et(
          target_model = s$model,
          prior_model  = model,
          event_table  = et_i
        )
        as.data.frame(s$model)
      }
    
      pred_i <- solve_pop(0) |>
        select(time, Cc_typical = Cc, AUC_typical = AUC) |>
        distinct(time, .keep_all = TRUE) |>
        arrange(time)
      
      if (pop_n_simul > 0) {
        sim_raw <- solve_pop(pop_n_simul)
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
    
    results[[id_key]][[model_name]] <- result_i
  }
}

combined_results <- list()
for (model_name in names(models_ceftazidime)) {
  combined_model_data <- do.call(rbind, lapply(results, `[[`, model_name))
  combined_results[[model_name]] <- combined_model_data
}

results_metrics <- lapply(combined_results, function(model_data) {
  model_data <- model_data %>% rename(predicted = Cc, observed = DV)
  model_data <- model_data %>% filter(!is.na(predicted) & !is.na(observed))
  calculate_metrics(model_data, n_bootstrap = 10000, conf_level = 0.95,
                    id_col = "ID")   # cluster bootstrap: resample patients, not samples
})


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
  
  ggsave(paste0("Ceftazidime/Concentration/POP/",i,".png"))
  
}

# MAP estimation ----------------------------------------------------------
map_data <- df_tdm1   ; map_label <- "1st TDM only"

results2 <- list()   # raw posologyr MAP fits, [[ID]][[model]]
m_out   <- list()    # predicted profiles, flat so bind_rows() still works
m_score <- list()    # one row per observation -- this is what metrics use
m_sir   <- list()    # SIR posterior eta DRAWS, one long tibble per (ID, model)

# --- SIR settings --------------------------------------------------------
sir_run        <- TRUE
sir_n_sample   <- 1e4   # etas proposed from the prior
sir_n_resample <- 1e3    # etas kept, resampled by importance weight

obs_df <- map_data |> 
  filter(EVID %in% c(0, 2), !is.na(DV)) |> 
  transmute(ID, TIME, OBSERVED = DV, used_in_fit = EVID == 0) |> 
  distinct()
#sort(unique(map_data$ID))

models_to_use <- names(models_ceftazidime)

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
      model <- models_ceftazidime[[model_name]]
      
      sim_results <- posologyr::poso_estim_map(patient_data_for_map,
                                               model,
                                               return_ofv   = TRUE,
                                               return_model = TRUE,
                                               nocb = FALSE)
      
      profile <- sim_results$model |> 
        as_tibble() |> 
        select(time, Cc, AUC) |> 
        mutate(ID = i)
      
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
          
          
          sir_fit <- posologyr::poso_estim_sir(patient_data_for_map, 
                                               model,
                                               n_sample     = sir_n_sample,
                                               n_resample   = sir_n_resample,
                                               return_model = FALSE,
                                               nocb = FALSE)
          
          eta_names_i <- names(sim_results$eta)
          stopifnot(all(eta_names_i %in% names(sir_fit$eta)))
          
          om_i    <- model$omega[eta_names_i, eta_names_i, drop = FALSE]
          om_sd_i <- setNames(sqrt(diag(om_i)), eta_names_i)
          
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

.omega_of <- function(model_name, eta) {
  stopifnot(!is.null(names(eta)), all(nzchar(names(eta))), !anyNA(names(eta)))
  om <- models_ceftazidime[[model_name]]$omega   # full prior covariance matrix
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
    aes(x = eta, xend = eta, y = 0, yend = peak, colour = outside),
    linewidth = 0.55, alpha = 0.9, inherit.aes = FALSE
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey35") +  
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
.SIR_MIN_UNIQUE <- 50

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
  
  out_dir <- "Ceftazidime/ETA/SIR"
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


# Metrics -----------------------------------------------------------------

m_out_tot <- bind_rows(m_out)

scoring_rows <- bind_rows(m_score) %>%
  select(dataset, model_name, ID, TIME, n_seen, used_in_fit, observed, predicted)

# Model averaging ---------------------------------------------------------
ma_weights <- function(log_w) {
  z <- log_w - max(log_w, na.rm = TRUE)
  w <- exp(z)
  w / sum(w, na.rm = TRUE)
}

# One ofv / SSE per fit, so the weights are per (patient, model) and constant
# over time. Models that failed for a patient are simply absent, and the
# weights renormalise over the ones that are left.
ma_fit_weights <- m_out_tot |>
  distinct(ID, model_name, ll_exp, aic_exp, sse_exp) |>
  group_by(ID) |>
  mutate(
    w_LL  = ma_weights(ll_exp),
    w_AIC = ma_weights(aic_exp),
    w_SSE = ma_weights(sse_exp)
  ) |>
  ungroup() |>
  select(ID, model_name, w_LL, w_AIC, w_SSE)

# Averaged prediction at each scored observation: the weighted mean of the
# single-model predictions there. Dividing by sum(w) renormalises when a model
# has no prediction at that time (observation outside its grid).
ma_scoring_rows <- scoring_rows |>
  filter(is.finite(predicted)) |>
  inner_join(ma_fit_weights, by = c("ID", "model_name")) |>
  group_by(dataset, ID, TIME, n_seen, used_in_fit, observed) |>
  summarise(
    MA_LL  = sum(predicted * w_LL)  / sum(w_LL),
    MA_AIC = sum(predicted * w_AIC) / sum(w_AIC),
    MA_SSE = sum(predicted * w_SSE) / sum(w_SSE),
    .groups = "drop"
  ) |>
  pivot_longer(starts_with("MA_"), names_to = "model_name", values_to = "predicted")

# Scored exactly like any other model from here on.
scoring_rows <- bind_rows(scoring_rows, ma_scoring_rows)

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


# Model averaging profiles ------------------------------------------------
# ma_weights() is defined in the Metrics section above.

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

df_final <- bind_rows(df_final,df_final_MA)

# group = ID everywhere: df_final now holds every patient, so without it the
# lines connect across patients. Facet on model_name + ID (or filter to one ID)
# if you want one panel per patient rather than all profiles overlaid.

for (i in sort(unique(df_final$ID))) {
  ggplot(df_final |> filter(ID == i), aes(x = time, y = Cc, group = ID)) +
    geom_line() +
    geom_point(data = obs_df |> filter(ID == i), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
    theme_minimal() +
    facet_wrap(~ model_name, scales = "free") 
  
  ggsave(paste0("Ceftazidime/Concentration/MAP/",i,".png"))
}
