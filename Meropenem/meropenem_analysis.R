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
    DUR = INFUSION_DURATION
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




# Event-table guard: resets inside a running infusion --------------------------
# An EVID 3/4 reset clears the compartments, but the infusion's scheduled
# off-event still fires afterwards and subtracts its rate from an already-zero
# rate -- leaving a NEGATIVE infusion running for the rest of the profile.
# ID 4: reset at t=309 sits inside the 120->312 infusion (83.3 mg/h), so from
# t=336 the net rate is 41.7 - 83.3 = -41.7 mg/h and Cc falls to -68 mg/L.
# A reset cannot physically occur mid-infusion, so drop those rows. Resets that
# fall between infusions are legitimate and are kept.
local({
  infusions <- df %>%
    dplyr::filter(EVID == 1, !is.na(DUR), DUR > 0) %>%
    transmute(ID, START = TIME, END = TIME + DUR)

  bad <- df %>%
    mutate(.rid = row_number()) %>%
    filter(EVID %in% c(3, 4)) %>%
    select(.rid, ID, TIME) %>%
    inner_join(infusions, by = "ID", relationship = "many-to-many") %>%
    filter(TIME > START, TIME < END) %>%
    distinct(.rid, ID, TIME)

  if (nrow(bad)) {
    message("Dropping ", nrow(bad), " reset event(s) falling inside a running infusion:")
    message(paste0("  ID ", bad$ID, " at t = ", bad$TIME, collapse = "
"))
    df <<- df %>% mutate(.rid = row_number()) %>%
      filter(!.rid %in% bad$.rid) %>% select(-.rid)
  }
})


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
  
  max_time_i <- max(patient_data_for_pop$TIME, na.rm = TRUE)
  
  df_i <- df |>  filter(ID == i)
  et_i <- rxode2::as.et(df_i)
  et_i$add.sampling(seq(0, max_time_i, by = 0.1)) 
  
  results[[i]] <- list()
  
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

plot_cc <- ggplot(combined_results_tot |> filter(ID == 21), aes(x = time, y = Cc, group = ID))+
  geom_smooth() +
  theme_minimal() +
  facet_wrap(~ model, ncol = 7, scales = "free") 
# MAP results -------------------------------------------------------------

results2 <- list()   # raw posologyr MAP fits, [[ID]][[model]]
m_out   <- list()    # tidy per-(ID, model) frames, flat so bind_rows() still works

obs_df <- df %>%
  filter(!is.na(DV)) %>%
  select(ID, TIME, OBSERVED = DV)

meropenem_models_to_consider <- setdiff(names(meropenem_models), c("Roberts_2009", "Ehmann_2019","Eisert_2021","Onichimowski_2020"))

for (i in sort(unique(df$ID))) {
  patient_data_for_map <- df %>%
    filter(ID == i)

  obs_i <- obs_df %>%
    filter(ID == i) %>%
    select(TIME, OBSERVED) %>%
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

      model_with_obs <- sim_results$model %>%
        select(time, Cc, AUC) %>%
        left_join(obs_i, by = c("time" = "TIME")) %>%
        mutate(ID = i)

      sse <- model_with_obs %>%
        filter(!is.na(OBSERVED)) %>%
        summarise(sse = sum((Cc - OBSERVED)^2)) %>%
        pull(sse)

      list(
        fit  = sim_results,
        tidy = model_with_obs %>%
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
            model_name = model_name
          )
      )

    }, error = function(e) {
      message("  Model ", model_name, " failed for ID ", i, ": ", conditionMessage(e))
      NULL
    })

    results2[[as.character(i)]][[model_name]] <- result2_i$fit
    m_out[[paste(i, model_name, sep = "_")]]  <- result2_i$tidy
  }
}


# .Machine
m_out_tot <- bind_rows(m_out)

# Softmax over the per-model exponents, computed WITHIN each patient: the fits
# are now per-ID, so LL from different patients must not be pooled.
softmax_by_model <- function(model_name, x) {
  key <- !duplicated(model_name)
  z   <- exp(x[key] - max(x[key], na.rm = TRUE))
  w   <- setNames(z / sum(z, na.rm = TRUE), model_name[key])
  unname(w[model_name])
}

df_final <- m_out_tot %>%
  group_by(ID) %>%
  mutate(
    weight_LL  = softmax_by_model(model_name, ll_exp),
    weight_AIC = softmax_by_model(model_name, aic_exp),
    weight_SSE = softmax_by_model(model_name, sse_exp)
  ) %>%
  group_by(ID, time) %>%
  mutate(
    MA_LL = sum(Cc * weight_LL),
    MA_AIC = sum(Cc * weight_AIC),
    MA_SSE = sum(Cc * weight_SSE)
  ) %>% 
  ungroup()  

df_final_MA <- df_final |> 
  select(ID,time,OBSERVED, starts_with("MA_")) |> 
  distinct() |>
  pivot_longer(cols = -c(ID, time, OBSERVED)) |> 
  rename(model_name = name, Cc = value)

df_final <- bind_rows(df_final,df_final_MA )

# group = ID everywhere: df_final now holds every patient, so without it the
# lines connect across patients. Facet on model_name + ID (or filter to one ID)
# if you want one panel per patient rather than all profiles overlaid.
ggplot(df_final |> filter(ID ==1), aes(x = time, y = Cc, group = ID)) +
  geom_line(colour  = "blue") +
  geom_point(data = obs_df |> filter(ID == 1), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  facet_wrap(~ model_name, ncol = 5, scales = "free") 
