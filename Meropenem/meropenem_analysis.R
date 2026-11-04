library(tidyverse); library(readxl)
library(data.table); library(ggrepel)
library(dtplyr); library(npde)
source("~/PhD-thesis/Utils/Functions.R")
source("~/PhD-thesis/Meropenem/meropenem_models.R")

total_covariates <- lapply(meropenem_models, function(x) x$covariates) %>% unlist() %>% unique()


df <- read_excel("~/PhD-thesis/Meropenem/data/nantes.xlsx") %>%
  mutate(across(where(is.character), ~na_if(.x, "."))) |> 
  readr::type_convert() |> 
  mutate(
    gender = if_else(SEX == 0, "female", "male", missing = NA_character_),
    WT     = WEIGHT, 
    height = HEIGHT, 
    TBW = WT,
    CREA_micromol_l = SCR,
    creatinine_mg_dl = CREA_micromol_l/88.4,
    ALB = ALB/10,
    DV = Y,
    DUR = INFUSION_DURATION,
    RRT = CRRT
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
    eGFR_CKD_EPI_CYSTATIN = calculate_egfr_cystatin_c_safe(
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender,
      serum_cystatin_c_mg_l  = cystatin_c
    ),
    eGFR_MDRD = calculate_egfr_mdrd_safe(     
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    RFS     = as.integer(eGFR_MDRD >= 120)
  ) %>%
  ungroup() |> 
  group_by(`Center ID`, ID) |> 
  mutate(ID = cur_group_id()) |> 
  ungroup()


# Extend each patient's event table past their last TDM ------------------------
# posologyr builds the solver grid from the events in the data, so a profile can
# stop at the last dose -- which in this dataset is ALWAYS before the last TDM
# (23/23 patients). Appending a placeholder event after each patient's final
# sample guarantees the grid reaches it.
#
# EVID = 1 with AMT = 0 and DUR = 0 is a no-op for the PK: it adds no drug, and
# it stays invisible to everything that counts real doses, because those all
# test AMT > 0 (the OCC index) or DUR > 0 (the reset guard). The last row is
# copied first so time-varying covariates carry forward (LOCF).
.n_before <- nrow(df)

df <- df %>%
  group_by(ID) %>%
  group_modify(function(g, key) {
    last_tdm <- suppressWarnings(max(g$TIME[g$EVID == 0 & !is.na(g$DV)]))
    if (!is.finite(last_tdm)) return(g)          # no TDM -> nothing to extend
    pad <- g %>%
      slice_max(TIME, n = 1, with_ties = FALSE) %>%
      mutate(TIME = last_tdm + 1, EVID = 1, AMT = 0, DUR = 0, DV = NA_real_)
    bind_rows(g, pad) %>% arrange(TIME)
  }) %>%
  ungroup()

# Count the difference, not `EVID == 1 & AMT == 0`: the source data already
# contains 4 such rows (ID 4, ID 5 x2, ID 18), so that signature would over-count.
message("Righe segnaposto aggiunte: ", nrow(df) - .n_before,
        " (una per paziente con almeno un TDM)")
rm(.n_before)


# patient_A <- data.frame(
#   ID  = 1,
#   OCC = 1,
#   AGE    = 65,     
#   WT     = 80,     
#   height = 175,    
#   gender = "male",
#   creatinine_mg_dl = 1.2, 
#   cystatin_c       = 0.8, 
#   ALB              = 3.5, 
#   WBC              = 12.0,
#   RRT    = 0,     
#   CRRT   = 0,     
#   CVVH   = 0,     
#   CVVHDF = 0,     
#   DIA_C  = 0,     
#   DIA_SC = 0,     
#   Qfd    = 0,     
#   Qf_ml_h  = 0,   
#   Sc_siev  = 0.8, 
#   CF       = 1.0, 
#   t_Dial   = 0,   
#   SEPSIS = 1, 
#   SEPTIC = 1, 
#   shock  = 0, 
#   ACLF   = 0, 
#   BURN   = 0, 
#   TBSA   = 0, 
#   ECMO = 0,    
#   LPM  = 1,    # ECMO blood flow, L/min; valid range > 0.733 L/min
#   DIUR = 1200,  # residual diuresis, mL/24h
#   HTA = 1,   # hypertension 
#   VAN = 0,   # concomitant vancomycin binary 
#   COL = 0,   # concomitant colistimethate   
#   IL6_CSF = 100,  # IL-6 in CSF, ng/L (= pg/mL)  
#   MIC  = 2,   # meropenem MIC, mg/L           
#   MRSA = 0,   # empirical anti-MRSA binary    
#   TIME = c(0,    1.5),
#   DV   = c(NA,  40.0),
#   AMT  = c(2000,  NA),
#   EVID = c(1,      0),
#   DUR  = c(0.5,   NA)
# ) %>%
#   rowwise() %>%
#   mutate(
#     TBW = WT,                            
#     IBW = calculate_ibw(
#       height_cm = height,
#       gender    = gender
#     ),
#     LBW = calculate_lbw(                  
#       height_cm = height,
#       weight_kg = WT,
#       gender    = gender
#     ),
#     LBM = LBW,                            
#     ABW = calculate_abw(
#       height_cm = height,
#       weight    = WT,
#       gender    = gender
#     ),
#     BSA = calculate_bsa(
#       weight_kg = WT,
#       height_cm = height
#     ),
#     CREA_micromol_l = convert_creatinine(
#       serum_creatinine_mg_dl = creatinine_mg_dl
#     ),
#     CLCR_CG = calculate_crcl(            
#       age                    = AGE,
#       weight_kg              = WT,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender
#     ),
#     CLCR_TBW    = CLCR_CG,               
#     CLCR        = CLCR_CG,               
#     CLCR_CG_LBW = calculate_crcl(
#       age                    = AGE,
#       weight_kg              = LBW,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender
#     ),
#     CLCR_CG_IBW = calculate_crcl(       
#       age                    = AGE,
#       weight_kg              = IBW,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender
#     ),
#     CLCR_J = calculate_crcl_jelliffe(
#       age              = AGE,
#       weight           = WT,
#       serum_creatinine = creatinine_mg_dl,
#       gender           = gender
#     ),
#     eGFR_CKD_EPI = calculate_egfr(
#       age                    = AGE,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender
#     ),
#     eGFR_CKD_EPI_CYSTATIN = calculate_egfr_cystatin_c(
#       age                    = AGE,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender,
#       serum_cystatin_c_mg_l  = cystatin_c
#     ),
#     eGFR_MDRD = calculate_egfr_mdrd(     #
#       age                    = AGE,
#       serum_creatinine_mg_dl = creatinine_mg_dl,
#       gender                 = gender
#     ),
#     RFS     = as.integer(eGFR_MDRD >= 120),   
#     CL_CVVH = calculate_cl_cvvh(            
#       Qf_ml_h = Qf_ml_h,
#       Sc      = Sc_siev,
#       CF      = CF
#     )
#   ) %>%
#   ungroup()


# How many tdm ------------------------------------------------------------

df_tdm <- df |> 
  group_by(ID) |> 
  summarise(
    tdm = sum(EVID == 0 & !is.na(DV), na.rm = TRUE),
    time = max(TIME[EVID == 0], na.rm = TRUE),
    time2 = mean(TIME[EVID == 0], na.rm = TRUE),
    time3 = median(TIME[EVID == 0], na.rm = TRUE)
  )

# Hide the TDM samples beyond the first `n_keep` of each patient.
#
# EVID 0 -> 2 keeps the row (and therefore the solver's time grid) but removes
# the point from the MAP likelihood. DV is deliberately left intact, so the true
# value is still there to score the prediction against.
# Do NOT mask by setting DV = NA (posologyr errors) or by deleting rows
# (the grid shortens and the hidden time can no longer be read back).
mask_tdm_after <- function(data, n_keep) {
  data |>
    arrange(ID, TIME) |>
    group_by(ID) |>
    mutate(
      .tdm_i = cumsum(EVID == 0 & !is.na(DV)),   # 1,2,3... over TDM rows only
      EVID   = replace(EVID, EVID == 0 & !is.na(DV) & .tdm_i > n_keep, 2)
    ) |>
    ungroup() |>
    select(-.tdm_i)
}

df_tdm1 <- mask_tdm_after(df, 1)   # only the first TDM visible
df_tdm2 <- mask_tdm_after(df, 2)   # first two visible


# Occasion index for the IOV models --------------------------------------------
# Ehmann_2019, Eisert_2021 and Roberts_2009 carry a pi_matrix, so posologyr needs
# OCC. An occasion is a MONITORED dosing interval: Ehmann_2019 (Table 2, note g)
# defines it as "monitored meropenem infusion, i.e. in total six occasions".
# Incrementing at every dose would open 278 occasions here (210 of them without a
# single TDM sample) and apply the IOV variance far more often than estimated.
df <- df %>%
  arrange(ID, TIME) %>%
  group_by(ID) %>%
  mutate(.interval = cumsum(EVID == 1 & !is.na(AMT) & AMT > 0)) %>%
  group_by(ID, .interval) %>%
  mutate(.sampled = any(EVID == 0 & !is.na(DV))) %>%
  group_by(ID) %>%
  mutate(OCC = pmax(cumsum(EVID == 1 & !is.na(AMT) & AMT > 0 & .sampled), 1L)) %>%
  ungroup() %>%
  select(-.interval, -.sampled)

message("OCC: ", sum(df$EVID == 1 & df$AMT > 0, na.rm = TRUE), " dosing intervals -> ",
        sum(tapply(df$OCC, df$ID, function(x) length(unique(x)))), " occasions across ",
        length(unique(df$ID)), " patients")


# a priori evaluation -----------------------------------------------------

doses <- df %>% 
  filter(EVID == 1) |>  
  select(
    ID, TIME, AMT, DUR
  ) 

population_df <- df |>  
  group_by(ID) |>  
  filter(
    TIME == 0
  ) %>% 
  mutate(
    AMT = 0,
    EVID = 0,
    DUR = NA
  )

results <- list()
m_out <- list() 

for (i in sort(unique(population_df$ID))) {
  message("Processing item ", i, " of ", length(unique(population_df$ID)))
  
  patient_data_for_pop <- population_df %>% 
    filter(ID == i) |> 
    mutate(DUR = as.numeric(NA))
  
  df_i <- df |>  filter(ID == i)
  max_time_i <- max(df_i$TIME, na.rm = TRUE)
  et_i <- rxode2::as.et(df_i)
  et_i$add.sampling(seq(0, max_time_i, by = 0.1)) 
  
  results[[as.character(i)]] <- list()
  
  for (model_name in names(meropenem_models)) {
    message("Testing model ", model_name, " for ID ", i)
    
    result_i <- tryCatch({
      model <- meropenem_models[[model_name]]
      
      sim_results <- posologyr::poso_simu_pop(
        patient_data_for_pop, 
        model, 
        n_simul = 0
      )
      
      sim_results$model <- posologyr::poso_replace_et(
        target_model = sim_results$model,
        prior_model  = model,
        event_table  = et_i
      )
      
      sim_results$model |>  
        select(time, Cc, AUC) |> 
        left_join(df_i %>% select(TIME, DV, ID), by = c("time" = "TIME")) |>  
        fill(ID)
      
    }, error = function(e) {
      message("  Model ", model_name, " failed for ID ", i, ": ", conditionMessage(e))
      NULL
    })
    
    results[[i]][[model_name]] <- result_i
  }
}

combined_results <- list()
for (model_name in names(meropenem_models)) {
  combined_model_data <- do.call(rbind, lapply(results, `[[`, model_name))
  combined_results[[model_name]] <- combined_model_data
}

results_metrics <- lapply(combined_results, function(model_data) {
  model_data <- model_data %>% rename(predicted = Cc, observed = DV)
  model_data <- model_data %>% filter(!is.na(predicted) & !is.na(observed))
  calculate_metrics(model_data, n_bootstrap = 1000, conf_level = 0.95,
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
      MDAPE_Upper_CI = metrics[["MDAPE_CI"]][["upper"]]
    )
  })
)

rownames(results_df) <- NULL

results_df <-results_df %>% 
  mutate(
    across(
      rBias:MDAPE_Upper_CI, ~ round(., 2)
    )
  )
  # ) |> 
  # filter_out(Model %in% c("Lee_2021", "Troisi_2024"))

plot_rBias <- ggplot(results_df, aes(x = Model, y = rBias, fill = Model)) +
  geom_point(show.legend = FALSE) +
  geom_errorbar(aes(ymin = rBias_Lower_CI, ymax = rBias_Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 30, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = -30, linetype = "dotted", color = "blue") +
  labs(title = "Relative Bias for Each Model",
       x = "",
       y = "rBias (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0,vjust = 0.5),
        axis.title.y = element_text(size = 18)) 


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
  facet_wrap(~model, ncol = 7)

obs_df <- df %>%
  filter(EVID %in% c(0), !is.na(DV)) %>%
  transmute(ID, TIME, OBSERVED = DV) %>%
  distinct()

plot_cc <- ggplot(combined_results_tot |> filter(ID == 20), aes(x = time, y = Cc, group = ID))+
  geom_line() +
  geom_point(data = obs_df |> filter(ID == 20), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  theme_minimal() +
  facet_wrap(~ model, ncol = 7, scales = "free") 
# MAP results -------------------------------------------------------------

# Swap these two lines to run a different design; everything below follows.
#   map_data <- df        ; map_label <- "all TDM"
   map_data <- df_tdm1   ; map_label <- "1st TDM only"
#   map_data <- df_tdm2   ; map_label <- "first 2 TDM"

results2 <- list()   # raw posologyr MAP fits, [[ID]][[model]]
m_out   <- list()    # predicted profiles, flat so bind_rows() still works
m_score <- list()    # one row per observation -- this is what metrics use

# EVERY true observation, flagged by whether the fit was allowed to see it.
# In a masked dataset the hidden samples keep their DV but carry EVID = 2, so
# posologyr leaves them out of the likelihood while the true value stays here
# for scoring. This flag is what makes a fair comparison possible downstream:
#   used_in_fit = TRUE   the fit saw it -> in-sample, optimistically biased
#   used_in_fit = FALSE  it was hidden  -> out-of-sample, the honest number
obs_df <- map_data %>%
  filter(EVID %in% c(0, 2), !is.na(DV)) %>%
  transmute(ID, TIME, OBSERVED = DV, used_in_fit = EVID == 0) %>%
  distinct()

meropenem_models_to_consider <- setdiff(names(meropenem_models), c("Roberts_2009", "Ehmann_2019","Eisert_2021","Onichimowski_2020"))

for (i in sort(unique(map_data$ID))) {
  patient_data_for_map <- map_data %>%
    filter(ID == i) %>%
    arrange(TIME)

  obs_i <- obs_df %>%
    filter(ID == i) %>%
    select(TIME, OBSERVED, used_in_fit) %>%
    distinct()

  results2[[as.character(i)]] <- list()

  # Iterate over each model in meropenem_models_to_consider
  for (model_name in meropenem_models_to_consider) {
    message("Testing model ", model_name, " for ID ", i)

    result2_i <- tryCatch({
      model <- meropenem_models[[model_name]]

      sim_results <- posologyr::poso_estim_map(patient_data_for_map, model,
                                               return_ofv   = TRUE,
                                               return_model = TRUE)

      # Full predicted profile, for the model-averaging plots.
      profile <- sim_results$model %>%
        as_tibble() %>%
        select(time, Cc, AUC) %>%
        mutate(ID = i)

      # One prediction per observation, by INTERPOLATION onto the solver grid.
      # Do NOT left_join on time: observation times need not land on a grid
      # point (ID 13's samples at 107.25, 179.25, 251.58, 353.92 fall between
      # the 0.1 h steps), and such rows silently became OBSERVED = NA and were
      # dropped from every metric -- losing the whole patient without a warning.
      # Duplicate grid times (dose + sample at the same instant) are collapsed to
      # the first, i.e. the pre-dose value, which is what a trough sample means.
      grid <- profile %>%
        distinct(time, .keep_all = TRUE) %>%
        arrange(time)

      scored_i <- obs_i %>%
        mutate(
          predicted = approx(grid$time, grid$Cc, xout = TIME, rule = 1)$y,
          ID = i, model_name = model_name, dataset = map_label,
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
      sse <- scored_i %>%
        filter(!is.na(predicted), used_in_fit) %>%
        summarise(sse = sum((predicted - observed)^2)) %>%
        pull(sse)

      list(
        fit    = sim_results,
        scored = scored_i,
        tidy   = profile %>%
          mutate(
            ofv        = sim_results$ofv,
            LL         = exp(-0.5 * ofv),
            AIC        = exp(-0.5 * ofv - nrow(model$omega)),
            SSE        = exp(-0.5 * sse),   # squared prediction errors
            # log-scale versions of the three above: exp() underflows to 0 for
            # OFV/SSE above ~1420, which would make every weight NaN
            ll_exp     = -0.5 * ofv,
            aic_exp    = -0.5 * ofv - nrow(model$omega),
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
  }
}

m_out_tot <- bind_rows(m_out)

# Keep each design's output before re-running with the next dataset, e.g.
#   m_out_all  <- m_out_tot   # after map_data <- df
#  m_out_tdm1 <- m_out_tot   # after map_data <- df_tdm1
# then compare with bind_rows(m_out_all, m_out_tdm1, ...) on the `dataset` column.

# Rows to score: exactly one per observation per (model, dataset).
# `used_in_fit == FALSE` are the predictions the fit never saw -- the only ones
# that give an honest out-of-sample metric.
scoring_rows <- bind_rows(m_score) %>%
  select(dataset, model_name, ID, TIME, n_seen, used_in_fit, observed, predicted)

# Sanity check: every model must score the SAME number of observations, or the
# metrics are not comparable across models.
scoring_rows %>% count(dataset, model_name, name = "n_scored") %>% print(n = Inf)

# Metrics ----------------------------------------------------------------------
# calculate_metrics() expects columns literally named `observed` / `predicted`,
# which scoring_rows already provides, plus id_col for the cluster bootstrap
# (resample PATIENTS, not observations -- several samples per patient are
# correlated and a row-level bootstrap gives CIs that are far too narrow).
#
# Which rows to score depends on the design:
#   masked dataset (df_tdm1, df_tdm2) -> out-of-sample rows exist, use those
#   full dataset   (df)               -> nothing was hidden, only in-sample
# The fallback is explicit so an in-sample number can never be mistaken for a
# validation result.
metrics_scope <- if (any(!scoring_rows$used_in_fit)) "out-of-sample" else "in-sample"
message("Metriche calcolate su punti ", metrics_scope,
        " (dataset: ", paste(unique(scoring_rows$dataset), collapse = ", "), ")")
if (metrics_scope == "in-sample")
  warning("Nessun TDM nascosto: queste metriche sono in-sample e ottimistiche, ",
          "non sono una validazione. Usa mask_tdm_after() per una stima onesta.")

metrics_input <- scoring_rows %>%
  filter(if (metrics_scope == "out-of-sample") !used_in_fit else used_in_fit,
         is.finite(observed), is.finite(predicted), observed > 0)

# CI fields are absent when a model has a single scored row (no bootstrap), so
# pull them defensively rather than indexing blind.
.ci <- function(m, field, side)
  if (is.null(m[[field]])) NA_real_ else unname(m[[field]][[side]])

metrics_df <- metrics_input %>%
  split(metrics_input$model_name) %>%
  imap_dfr(function(d, nm) {
    m <- calculate_metrics(d, n_bootstrap = 1000, conf_level = 0.95, id_col = "ID")
    tibble(
      Model    = nm,
      n_scored = nrow(d),
      n_pts    = n_distinct(d$ID),
      rBias    = m$rBias, rBias_lo = .ci(m, "rBias_CI", "lower"),
                          rBias_hi = .ci(m, "rBias_CI", "upper"),
      rRMSE    = m$rRMSE, rRMSE_lo = .ci(m, "rRMSE_CI", "lower"),
                          rRMSE_hi = .ci(m, "rRMSE_CI", "upper"),
      MDPE     = m$MDPE, MDAPE = m$MDAPE, F20 = m$F20, F30 = m$F30
    )
  }) %>%
  mutate(dataset = paste(unique(scoring_rows$dataset), collapse = ", "),
         scope   = metrics_scope,
         across(where(is.numeric), ~ round(.x, 2))) %>%
  arrange(MDAPE)

print(metrics_df, n = Inf)


# Plots --------------------------------------------------------------------
# Same grammar as the population plots above so the two sets can sit side by
# side. The subtitle carries dataset + scope on every panel: an in-sample figure
# must never be mistaken for a validation result once it leaves the script.
map_subtitle <- paste0(unique(metrics_df$dataset), "  |  ", unique(metrics_df$scope),
                       "  |  ", n_distinct(scoring_rows$ID), " patients ")

map_theme <- theme_classic() +
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y   = element_text(size = 18),
        plot.title    = element_text(size = 22, hjust = 0, vjust = 0.5),
        plot.subtitle = element_text(size = 13, colour = "grey30"),
        axis.title.y  = element_text(size = 18),
        legend.text   = element_text(size = 14),
        legend.title  = element_blank())

# rBias: dashed line at 0, dotted guides at +/-30% (as in the population plot)
plot_map_rBias <- metrics_df %>%
  mutate(Model = fct_reorder(Model, rBias)) %>%
  ggplot(aes(x = Model, y = rBias)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = c(-30, 30), linetype = "dotted", colour = "blue") +
  geom_errorbar(aes(ymin = rBias_lo, ymax = rBias_hi), width = 0.2) +
  geom_point(aes(fill = Model), shape = 21, size = 3, show.legend = FALSE) +
  labs(title = "Relative bias (MAP)", subtitle = map_subtitle,
       x = "", y = "rBias (%)") +
  map_theme

# rRMSE: strictly positive, so only the upper guide is meaningful
plot_map_rRMSE <- metrics_df %>%
  mutate(Model = fct_reorder(Model, rRMSE)) %>%
  ggplot(aes(x = Model, y = rRMSE)) +
  geom_hline(yintercept = 30, linetype = "dotted", colour = "blue") +
  geom_errorbar(aes(ymin = rRMSE_lo, ymax = rRMSE_hi), width = 0.2) +
  geom_point(aes(fill = Model), shape = 21, size = 3, show.legend = FALSE) +
  expand_limits(y = 0) +
  labs(title = "Relative RMSE (MAP)", subtitle = map_subtitle,
       x = "", y = "rRMSE (%)") +
  map_theme

# F20 / F30: share of predictions within +/-20% and +/-30% of the observation.
# Ordered by F30 so the best-performing models sit to the right.
plot_map_F <- metrics_df %>%
  select(Model, F20, F30) %>%
  pivot_longer(c(F20, F30), names_to = "metric", values_to = "pct") %>%
  mutate(Model = fct_reorder(Model, pct, .fun = max)) %>%
  ggplot(aes(x = Model, y = pct, fill = metric)) +
  geom_col(position = position_dodge(width = 0.5), width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dotted", colour = "blue") +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, .05))) +
  scale_fill_manual(values = c(F20 = "#2c7fb8", F30 = "#a6bddb")) +
  labs(title = "Predictions within +/-20% and +/-30% (MAP)", subtitle = map_subtitle,
       x = "", y = "% delle predizioni") +
  map_theme + theme(legend.position = "top")

# ggsave("map_metrics.png", combined_map_metrics, width = 14, height = 16, dpi = 300)


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
ggplot(df_final |> filter(ID == 1), aes(x = time, y = Cc, group = ID)) +
  geom_line(colour  = "blue") +
  geom_point(data = obs_df |> filter(ID == 1), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  facet_wrap(~ model_name, ncol = 5, scales = "free") 