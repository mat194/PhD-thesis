library(tidyverse)
library(readxl)

# Build a posologyr/rxode2-ready event table from the "dataset pk singolo
# pz (rimossi3" sheet: a loading dose + continuous-infusion piperacillin
# regimen and one steady-state TDM sample per patient.
#
# Column mapping (verified against the sheet's own pre-computed clearance,
# drug_cl_l_hour == dose_pipera*1000/24 / conc_first_tdm for every patient):
#   dose_pipera     -> continuous infusion rate, g/24h (maintenance dose)
#   posology_daily  -> day-1 total dose, g. Exactly dose_pipera*1.125 for
#                      every patient in the cohort -> a fixed loading bolus
#                      of dose_pipera/8 g on top of the maintenance rate.
#                      (`posology` is free text for dose-adjusted patients,
#                      e.g. "18, poi 13.5" -> use posology_daily instead.)
#   conc_first_tdm  -> total plasma concentration (mg/L) at the TDM sample (DV)
#   cl_crea_tdm1    -> Cockcroft-Gault CrCL (mL/min) measured at the TDM, i.e.
#                      contemporary with the sample -> CLCR_CG covariate
#   start_therapy   -> t = 0, loading-dose start
#   data_first_tdm  -> sampling time (day-level resolution only; no time-of-day
#                      is recorded, so TIME is in whole days converted to hours)
#   end_therapy     -> used to size the maintenance-infusion AMT so it covers
#                      the whole treatment course (falls back to the sampling
#                      time when end_therapy is missing)
#
# Dosing structure per patient (mirrors the 0.5h load / 23.5h+ maintenance
# split used elsewhere in this project): a loading bolus of dose_pipera/8 g
# infused over 30 minutes at t=0, immediately followed by the continuous
# maintenance infusion at dose_pipera g/24h for the remainder of the course.
# The 30-minute loading duration is not in the data and is assumed from
# standard protocol -- adjust LOADING_DUR_H below if that's wrong.
LOADING_DUR_H <- 0.5

raw <- read_excel(
  "Piperacillin/Udine_paper/data_udine_pip.xlsx",
  sheet = "dataset pk singolo pz (rimossi3"
) %>%
  janitor::clean_names() %>%
  filter(!is.na(idx)) # drop blank padding rows

stopifnot(!anyDuplicated(raw$id_paz))

pk_data <- raw %>%
  transmute(
    ID     = id_paz,
    idx,
    AGE    = floor(age),
    SEX    = sex_1_m, # 1 = male, 0 = female
    WT     = weight,
    HT     = height,
    CLCR_CG = cl_crea_tdm1,
    start_therapy,
    end_therapy,
    data_first_tdm,
    dose_pipera,
    conc_first_tdm
  ) %>%
  mutate(
    loading_g       = 4,
    loading_rate_h  = loading_g * 1000 / LOADING_DUR_H,
    maint_rate_h    = dose_pipera * 1000 / 24,
    obs_time_h      = as.numeric(difftime(data_first_tdm, start_therapy, units = "hours")),
    dur_h           = as.numeric(difftime(end_therapy, start_therapy, units = "hours")),
    dur_h           = if_else(is.na(dur_h) | dur_h < obs_time_h, obs_time_h, dur_h),
    maint_dur_h     = pmax(dur_h - LOADING_DUR_H, 0)
  )

loading_rows <- pk_data %>%
  transmute(
    ID, idx, AGE, SEX, WT, HT, CLCR_CG,
    TIME = 0,
    EVID = 1,
    OCC  = 1,
    AMT  = loading_g * 1000,
    RATE = loading_rate_h,
    DV   = NA_real_
  )

maint_rows <- pk_data %>%
  transmute(
    ID, idx, AGE, SEX, WT, HT, CLCR_CG,
    TIME = LOADING_DUR_H,
    EVID = 1,
    OCC  = 1,
    AMT  = maint_rate_h * maint_dur_h,
    RATE = maint_rate_h,
    DV   = NA_real_
  )

obs_rows <- pk_data %>%
  transmute(
    ID, idx, AGE, SEX, WT, HT, CLCR_CG,
    TIME = obs_time_h,
    EVID = 0,
    OCC  = 1,
    AMT  = NA_real_,
    RATE = NA_real_,
    DV   = conc_first_tdm
  )

population_df <- bind_rows(loading_rows, maint_rows, obs_rows) %>%
  arrange(ID, TIME, desc(EVID)) %>%
  as.data.frame()


saveRDS(population_df, "Piperacillin/Udine_paper/population_df.rds")
