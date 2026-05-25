library(tidyverse)
library(data.table)
library(dtplyr)
source("~/PhD-thesis/Utils/Functions.R")
source("~/PhD-thesis/Meropenem/meropenem_models.R")

total_covariates <- lapply(meropenem_models, function(x) x$covariates) %>% unlist() %>% unique()

patient_A <- data.frame(
  ID  = 1,
  OCC = 1,
  AGE    = 65,     
  WT     = 80,     
  height = 175,    
  gender = "male",
  creatinine_mg_dl = 1.2, 
  cystatin_c       = 0.8, 
  ALB              = 3.5, 
  WBC              = 12.0,
  RRT    = 0,     
  CRRT   = 0,     
  CVVH   = 0,     
  CVVHDF = 0,     
  DIA_C  = 0,     
  DIA_SC = 0,     
  Qfd    = 0,     
  Qf_ml_h  = 0,   
  Sc_siev  = 0.8, 
  CF       = 1.0, 
  t_Dial   = 0,   
  SEPSIS = 1, 
  SEPTIC = 1, 
  shock  = 0, 
  ACLF   = 0, 
  BURN   = 0, 
  TBSA   = 0, 
  ECMO = 0,    
  LPM  = 1,    # ECMO blood flow, L/min; valid range > 0.733 L/min
  DIUR = 1200,  # residual diuresis, mL/24h
  HTA = 1,   # hypertension 
  VAN = 0,   # concomitant vancomycin binary 
  COL = 0,   # concomitant colistimethate   
  IL6_CSF = 100,  # IL-6 in CSF, ng/L (= pg/mL)  
  MIC  = 2,   # meropenem MIC, mg/L           
  MRSA = 0,   # empirical anti-MRSA binary    
  TIME = c(0,    1.5),
  DV   = c(NA,  40.0),
  AMT  = c(2000,  NA),
  EVID = c(1,      0),
  DUR  = c(0.5,   NA)
) %>%
  rowwise() %>%
  mutate(
    TBW = WT,                            
    IBW = calculate_ibw(
      height_cm = height,
      gender    = gender
    ),
    LBW = calculate_lbw(                  
      height_cm = height,
      weight_kg = WT,
      gender    = gender
    ),
    LBM = LBW,                            
    ABW = calculate_abw(
      height_cm = height,
      weight    = WT,
      gender    = gender
    ),
    BSA = calculate_bsa(
      weight_kg = WT,
      height_cm = height
    ),
    CREA_micromol_l = convert_creatinine(
      serum_creatinine_mg_dl = creatinine_mg_dl
    ),
    CLCR_CG = calculate_crcl(            
      age                    = AGE,
      weight_kg              = WT,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_TBW    = CLCR_CG,               
    CLCR        = CLCR_CG,               
    CLCR_CG_LBW = calculate_crcl(
      age                    = AGE,
      weight_kg              = LBW,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_CG_IBW = calculate_crcl(       
      age                    = AGE,
      weight_kg              = IBW,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    CLCR_J = calculate_crcl_jelliffe(
      age              = AGE,
      weight           = WT,
      serum_creatinine = creatinine_mg_dl,
      gender           = gender
    ),
    eGFR_CKD_EPI = calculate_egfr(
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    eGFR_CKD_EPI_CYSTATIN = calculate_egfr_cystatin_c(
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender,
      serum_cystatin_c_mg_l  = cystatin_c
    ),
    eGFR_MDRD = calculate_egfr_mdrd(     #
      age                    = AGE,
      serum_creatinine_mg_dl = creatinine_mg_dl,
      gender                 = gender
    ),
    RFS     = as.integer(eGFR_MDRD >= 120),   
    CL_CVVH = calculate_cl_cvvh(            
      Qf_ml_h = Qf_ml_h,
      Sc      = Sc_siev,
      CF      = CF
    )
  ) %>%
  ungroup()

results <- list()
m_out <- list() 

obs_df <- patient_A %>%
  filter(!is.na(DV)) %>%
  select(TIME, OBSERVED = DV)

for (i in names(meropenem_models)) {
  results[[i]] <- posologyr::poso_estim_map(
    dat = patient_A, 
    prior_model  = meropenem_models[[i]],
    return_ofv = TRUE,
    return_model = TRUE
  )
  
  model_with_obs <- results[[i]]$model %>%
    select(time, Cc, AUC) %>%
    left_join(obs_df, by = c("time" = "TIME"))
  
  sse <- model_with_obs %>%
    filter(!is.na(OBSERVED)) %>%
    summarise(sse = sum((Cc - OBSERVED)^2)) %>%
    pull(sse)
  
  m_out[[i]] <- model_with_obs %>% 
    mutate(
      ofv = results[[i]]$ofv,
      LL = exp(-0.5 * ofv),
      AIC = exp(-0.5 * ofv - nrow(meropenem_models[[i]]$omega)),
      SSE = exp(-0.5 * sse), # squared prediction errors
      model_name = i
    )
}

# .Machine
m_out_tot <- bind_rows(m_out)

df_final <- m_out_tot %>%
  mutate(
    weight_LL = LL / sum(unique(LL)),
    weight_AIC = AIC / sum(unique(AIC)),
    weight_SSE = SSE / sum(unique(SSE))
  ) %>%
  group_by(time) %>%
  mutate(
    MA_LL = sum(Cc * weight_LL),
    MA_AIC = sum(Cc * weight_AIC),
    MA_SSE = sum(Cc * weight_SSE)
  ) %>% 
  ungroup()  

ggplot(df_final, aes(x = time, y = Cc)) +
  geom_line(colour  = "blue") +
  geom_line(aes(y=MA_LL),size = 1, colour  = "red") +
  geom_line(aes(y=MA_AIC),size = 1, colour  = "green") +
  geom_line(aes(y=MA_SSE),size = 1, colour  = "orange") +
  geom_point(data = obs_df, aes(x = TIME, y = OBSERVED), size = 3, colour  = "brown", na.rm = TRUE) +
  theme_minimal() +
  facet_wrap(~ model_name, ncol = 5, scales = "free") 
