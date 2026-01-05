library(posologyr)
library(tidyverse)
library(readxl)
library(ggrepel)
library(flextable)
library(gtsummary)

# Funzioni ----------------------------------------------------------------

 calculate_egfr <- function(age, serum_creatinine_mg_dl, gender) {
   if (gender == "male") {
     kappa <- 0.9
     alpha <- -0.302
     min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
     max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
     egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age)
   } else {
     kappa <- 0.7
     alpha <- -0.241
     min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
     max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
     egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age) * 1.012
   }
   return(egfr)
 }
 
df <- read_excel("Vancomycin/Data/data_bloomy.xlsx") %>% filter(!is.na(LD))


calculate_crcl <- function(age, weight_kg, serum_creatinine_mg_dl, gender) {
  CrCl <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
  if (gender == "female") {
    CrCl <- CrCl * 0.85
  }
  return(CrCl)
}


df_patients <- df %>% 
  rowwise() %>%  # Apply the function row-wise
  mutate(
    CLCREAT = ifelse(
      is.na(creatinina) | creatinina == 0, NA,  # Handle missing or zero creatinine values
      calculate_crcl(
        age =age,  # Replace with actual age if available
        weight_kg = peso,
        serum_creatinine_mg_dl = creatinina,
        gender = sex
      )
    )
  ) %>%
  ungroup() %>% 
  group_by(record_id) %>% 
  mutate(
    ID =cur_group_id() , #dplyr way to get the index column working easily: gives a unique numeric identifier for the current group.
  ) %>% 
  ungroup() %>% 
  rename(
    WT = peso,
    DIAL = dialisi
  )


# Statistiche -------------------------------------------------------------

patients <- df_patients %>% 
  filter(!is.na(CLCREAT)) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull()

ci <- df_patients %>% 
  filter(DUR == 24) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull()
 
ci2 <- df_patients %>% 
  filter(DUR == 24) %>% 
  select(record_id) %>% 
  distinct() %>% 
  pull()


# Modelli -----------------------------------------------------------------

models_vancomycin <- list(
  Goti2018 = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
      TVVc  = THETA_Vc*(WT/70)*(0.5^DIAL);
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.5, THETA_Vc=58.4, THETA_Vp=38.4,THETA_Q=6.5),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(0.147,
          0, 0.510,
          0, 0, 0.282,
          0, 0, 0, 0)}),
    sigma = c(additive_a = 3.4, proportional_b = 0.227),
    covariates = c("CLCREAT","WT","DIAL")
  ),
  
  Roberts2011 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/100);
      TVVc  = THETA_Vc*(WT);
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.58, THETA_Vc=1.53),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c(0.375379779^2,
          0, 0.361827976^2)}),
    sigma = c(additive_a = 2.4, proportional_b = 0.227),
    covariates = c("CLCREAT","WT")
  ),
  
  Buelga2015 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl * CLCREAT * 0.06;# CrCL is in L/h not ml/min
      TVVc  = THETA_Vc * WT;
      Cl    = TVCl * exp(ETA_Cl);
      Vc    = TVVc * exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=1.08, THETA_Vc=0.98),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((28.16/ 100)^2 + 1),
           0, log((37.15/ 100)^2 + 1))}),
    sigma = c(additive_a = 3.52, proportional_b = 0),
    covariates = c("CLCREAT","WT")
  ),

  Adane2015 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl * (CLCREAT/125);
      TVVc  = THETA_Vc * WT;
      Cl    = TVCl * exp(ETA_Cl);
      Vc    = TVVc * exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl= 6.54, THETA_Vc= 0.51),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((26.7/ 100)^2 + 1),
           0, log((23.9/ 100)^2 + 1))}),
    sigma = c(additive_a = 0, proportional_b = 0.189),
    covariates = c("CLCREAT","WT")
  ),
  Lopis2006 = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl1 * CLCREAT + THETA_Cl2 * WT ;
      TVVc  = THETA_Vc * WT;
      TVVp  = THETA_Vp * WT;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl1=0.034, THETA_Cl2=0.015, THETA_Vc=0.414, THETA_Vp=1.32,THETA_Q=7.48),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp ~
        c( log((29.2/ 100)^2 + 1),
           0, log((36.4/ 100)^2 + 1),
           0, 0, log((39.8/ 100)^2 + 1))}),
    sigma = c(additive_a = 23.9, proportional_b = 0.185),
    covariates = c("CLCREAT","WT")
  ),
  Okada2018  = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/113)^0.78;
      TVVc  = THETA_Vc*(WT/59.4)^0.70;
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.25, THETA_Vc=39.2, THETA_Vp=56.1,THETA_Q=1.95),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp  ~
        c( log((25.2/ 100)^2 + 1),
           0, log((14.2/ 100)^2 + 1),
           0, 0, log((66.9/ 100)^2 + 1))}),
    sigma = c(additive_a = 0, proportional_b = 0.172),
    covariates = c("CLCREAT","WT")
  ),
    Donadello2014 = list(
      ppk_model   = rxode2::rxode2({
        centr(0) = 0;
        if (CRRT == 1){
          TVCl  = 0.6;
        } else {
          TVCl  = 1;
        }
        TVVc  = THETA_Vc;
        TVVp  = THETA_Vp;
        TVQ   = THETA_Q;
        Cl    = THETA_Cl * TVCl *exp(ETA_Cl);
        Vc    = TVVc*exp(ETA_Vc);
        Vp    = TVVp*exp(ETA_Vp);
        Q     = TVQ;
        ke    = Cl/Vc;
        k12   = Q/Vc;
        k21   = Q/Vp;
        Cc    = centr/Vc;
        d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
        d/dt(periph) = + k12*centr - k21*periph;
        d/dt(AUC)    = Cc;
      }),
      error_model = function(f,sigma){
        g <- sigma[1] + sigma[2]*f
        return(g)
      },
      theta = c(THETA_Cl = 3.7, THETA_Vc=31.8, THETA_Vp=57.1,THETA_Q=3.6),
      omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp  ~
          c( log((16.4/ 100)^2 + 1),
             0, log((47/ 100)^2 + 1),
             0, 0, log((101/ 100)^2 + 1))}),
      sigma = c(additive_a = 0, proportional_b = 0.085),
      covariates = c("CRRT")
    ),
  Medellin2017 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      if (vm == 0){
        TVCl  = (CLCREAT/100)^0.75;
      } else {
        TVCl  = (CLCREAT/100)^0.75 * 0.8;
      }
      TVVc  = THETA_Vc * WT;
      Cl    = THETA_Cl * TVCl *exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=2.86, THETA_Vc=1.03),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((28.4/ 100)^2 + 1),
           0, log((49.1/ 100)^2 + 1))}),
    sigma = c(additive_a = 4.3, proportional_b = 0),
    covariates = c("CLCREAT","WT","vm")
  )
  )
  # Thomson2009 = list(
  #   ppk_model   = rxode2::rxode2({
  #     centr(0)  = 0;
  #     periph(0) = 0;
  #     AUC(0)    = 0;
  #     TVCl  = THETA_Cl * (1 + 0.0154 * (CLCREAT - 66));
  #     TVVc  = THETA_Vc;
  #     TVVp  = THETA_Vp;
  #     TVQ   = THETA_Q;
  #     Cl    = TVCl*exp(ETA_Cl);
  #     Vc    = TVVc*exp(ETA_Vc);
  #     Vp    = TVVp*exp(ETA_Vp);
  #     Q     = TVQ*exp(ETA_Q);
  #     ke    = Cl/Vc;
  #     k12   = Q/Vc;
  #     k21   = Q/Vp;
  #     Cc    = centr/Vc;
  #     d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
  #     d/dt(periph) = + k12*centr - k21*periph;
  #     d/dt(AUC)    = Cc;
  #   }),
  #   error_model = function(f,sigma){
  #     g <- sigma[1] + sigma[2]*f
  #     return(g)
  #   },
  #   theta = c(THETA_Cl=2.99, THETA_Vc=0.732, THETA_Vp=0.675,THETA_Q=2.28),
  #   omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
  #       c( log((27/ 100)^2 + 1),
  #          0, log((15/ 100)^2 + 1),
  #          0, 0, log((130/ 100)^2 + 1),
  #          0, 0, 0, log((49/ 100)^2 + 1))}),
  #   sigma = c(additive_a = 1.6, proportional_b = 0.15),
  #   covariates = c("CLCREAT")
  # )


# df_patients <- data.frame(
#   ID = rep(1:100, each = 4),  # Ten patients, each with 4 time points
#   TIME = rep(c(0.0, 13.0, 24.2, 48.0), times = 100),
#   DV = c(NA, 12, NA, 9.5) + rnorm(40, 0, 1),  # Random variation in DV values
#   AMT = rep(c(2000, 0, 1000, 0), times = 100),
#   DUR = rep(c(2, NA, 2, NA), times = 100),
#   EVID = rep(c(1, 0, 1, 0), times = 100),
#   CLCREAT = round(runif(10, 50, 80)),  # Random creatinine clearance values
#   WT = round(runif(10, 60, 90)),  # Random weights
#   DIAL = sample(c(0, 1), 10, replace = TRUE)  # Random dialysis status
# )


# Population simulation ---------------------------------------------------
doses <- df_patients %>% 
  filter(!is.na(CLCREAT)) %>% 
  filter(EVID == 1) %>% 
  select(
    ID, TIME, AMT, DUR
  ) 

population_df <- df_patients %>% 
  group_by(ID) %>% 
  filter(
    TIME == 0
  ) %>% 
  mutate(
    AMT = 0,
    EVID = 0,
    DUR = NA
  )

results <- list()
max_time <- max(df_patients$TIME, na.rm = TRUE)

for (i in unique(population_df$ID)) {
  # Filter data for the current patient
  patient_data_for_pop <- population_df %>% 
    filter(ID == i) %>% 
    mutate(
      DUR = as.numeric(NA)
    )
  
  patient_data_for_map <- df_patients %>% 
    filter(ID == i)
  
  df_i <- df_patients %>% 
    filter(ID == i)
  
  doses_i <- doses %>% 
    filter(ID == i)
  results[[i]] <- list()
  
  # Iterate over each model in models_vancomycin
  for (model_name in names(models_vancomycin)) {
    model <- models_vancomycin[[model_name]]
    
    sim_results <- posologyr::poso_simu_pop(
      patient_data_for_pop, 
      model, 
      n_simul = 0
    )
    
    sim_results$model$time <- seq(0, max_time, b=0.1)
    
    for (j in 1:nrow(doses_i)) {
      dose <- doses_i[j, ]
      
      sim_results$model$add.dosing(
        dose = dose$AMT, 
        start.time = dose$TIME, 
        nbr.doses = 1, 
        rate = dose$AMT / dose$DUR, 
        dosing.interval = 0, 
        dosing.to = 1L
      )
    }
      
    results[[i]][[model_name]] <- sim_results$model %>% select(time, Cc, AUC) %>%
      left_join(df_i %>% select(TIME, DV, ID), by = c("time" = "TIME")) %>% 
      fill(ID)
    
  }
}

# The function being applied is [[, which is used for extracting elements from a 
# list. Here, model_name is provided as an argument to [[, so it extracts the data 
# corresponding to model_name from each patient's sublist.
# As a result, lapply(results, [[, model_name) returns a list where each element is the data frame for the specified model_name for each patient.
combined_results <- list()
for (model_name in names(models_vancomycin)) {
  combined_model_data <- do.call(rbind, lapply(results, `[[`, model_name))
  combined_results[[model_name]] <- combined_model_data
}

# combined_results <- combined_results %>% 
#   map(~filter(., !is.na(DV)))

# RMSE and r BIAS calculation ---------------------------------------------

#Definizioni prese dal papaer di Wicha sulla vancomicina
calculate_metrics <- function(data, n_bootstrap = 0, conf_level = 0.95) {
  
  N <- nrow(data)
  
  # Function to calculate metrics for a single dataset
  calculate_single <- function(data) {
    # Calculate rBias
    rBias <- (1 / N) * sum((data$predicted - data$observed) / ((data$predicted + data$observed) / 2)) * 100
    
    # Calculate rRMSE
    rRMSE <- sqrt((1 / N) * sum(((data$predicted - data$observed)^2) / (((data$predicted + data$observed) / 2)^2))) * 100
    
    # Calculate rPE for each observation
    data <- data %>%
      mutate(rPE = ((predicted - observed) / (observed / 2)) * 100)
    
    # Calculate MPE (median of rPE)
    MPE <- median(data$rPE, na.rm = TRUE)
    
    # Calculate MAPE (median of absolute rPE)
    MAPE <- median(abs(data$rPE), na.rm = TRUE)
    
    # Return calculated metrics
    c(rBias = rBias, rRMSE = rRMSE, MPE = MPE, MAPE = MAPE)
  }
  
  # Calculate metrics for the original data
  original_metrics <- calculate_single(data)
  
  # Bootstrapping for confidence intervals if n_bootstrap > 0
  if (n_bootstrap > 0 && N > 1) { 
    # Perform bootstrapping
    bootstrap_metrics <- replicate(n_bootstrap, {
      sampled_data <- data[sample(1:N, N, replace = TRUE), ]
      calculate_single(sampled_data)
    }, simplify = TRUE)
    
    # Calculate confidence intervals for each metric
    ci_lower <- apply(bootstrap_metrics, 1, function(x) quantile(x, (1 - conf_level) / 2, na.rm = TRUE))
    ci_upper <- apply(bootstrap_metrics, 1, function(x) quantile(x, 1 - (1 - conf_level) / 2, na.rm = TRUE))
    
    # Combine original metrics with confidence intervals
    results <- list(
      rBias = original_metrics["rBias"],
      rBias_CI = c(lower = ci_lower["rBias"], upper = ci_upper["rBias"]),
      rRMSE = original_metrics["rRMSE"],
      rRMSE_CI = c(lower = ci_lower["rRMSE"], upper = ci_upper["rRMSE"]),
      MPE = original_metrics["MPE"],
      MPE_CI = c(lower = ci_lower["MPE"], upper = ci_upper["MPE"]),
      MAPE = original_metrics["MAPE"],
      MAPE_CI = c(lower = ci_lower["MAPE"], upper = ci_upper["MAPE"])
    )
  } else {
    # Return original metrics without confidence intervals
    results <- list(
      rBias = original_metrics["rBias"],
      rRMSE = original_metrics["rRMSE"],
      MPE = original_metrics["MPE"],
      MAPE = original_metrics["MAPE"]
    )
  }
  
  return(results)
}

# Apply this function to each dataset in combined_results
results_metrics <- lapply(combined_results, function(model_data) {
  model_data <- model_data %>% rename(predicted = Cc, observed = DV)
  model_data <- model_data %>% filter(!is.na(predicted) & !is.na(observed))
  
  # Calculate metrics with bootstrapping (set n_bootstrap > 0 to enable)
  calculate_metrics(model_data, n_bootstrap = 1000, conf_level = 0.95)
})

# Convert results to a data frame for plotting or further analysis

# The performance of the models was considered clinically acceptable if the rBias 
# was between −20% and 20%, with the 95% confidence intervals (CI) including zero.
# Additionally, the precision metric (rRMSE) should be as small as possible.

results_df <- bind_rows(
  lapply(names(results_metrics), function(model_name) {
    metrics <- results_metrics[[model_name]]
    
    data.frame(
      Model = model_name,
      rBias = metrics$rBias,
      rBias_Lower_CI = metrics[["rBias_CI"]][["lower.rBias"]],
      rBias_Upper_CI =  metrics[["rBias_CI"]][["upper.rBias"]],
      rRMSE = metrics$rRMSE,
      rRMSE_Lower_CI = metrics[["rRMSE_CI"]][["lower.rRMSE"]],
      rRMSE_Upper_CI =  metrics[["rRMSE_CI"]][["upper.rRMSE"]],
      MPE = metrics$MPE,
      MPE_Lower_CI =  metrics[["MPE_CI"]][["lower.MPE"]],
      MPE_Upper_CI =  metrics[["MPE_CI"]][["upper.MPE"]],
      MAPE = metrics$MAPE,
      MAPE_Lower_CI =  metrics[["MAPE_CI"]][["lower.MAPE"]],
      MAPE_Upper_CI =  metrics[["MAPE_CI"]][["upper.MAPE"]]
    )
  })
)

rownames(results_df) <- NULL

results_df <-results_df %>% 
  mutate(
    across(
      rBias:MAPE_Upper_CI, ~ round(., 2)
    )
  ) 

# Plots -------------------------------------------------------------------

# Plot rBias for each model
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

df_goti <- combined_results$Goti2018 %>% 
  filter(!is.na(DV)) %>%
  mutate(
    Cc = round(Cc)
  )

ggplot(df_goti, aes(x = DV, y = Cc)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 10) +  # Add labels with ID column
  labs(title = "Model Fit for Goti2018",
       x = "Observed",
       y = "Predicted") +
  coord_cartesian(xlim = c(0, max(df_goti$DV, df_goti$Cc)), 
                  ylim = c(0, max(df_goti$DV, df_goti$Cc))) +
  theme_classic()


# For the best model let's finda patients with worts precision ------------

best_model <- results_df %>% 
  filter(rRMSE == min(rRMSE)) %>% 
  pull(Model)

worst_precision <- combined_results[[best_model]] %>% 
  filter(!is.na(DV)) %>%
  mutate(
    Precision = abs(DV - Cc) / DV
  ) %>% 
  filter(Precision > 0.3) %>% 
  select(ID) %>% 
  distinct() %>% 
  pull()

# MAP results -------------------------------------------------------------

results2 <- list()

df_patients <- df_patients %>% 
 filter(!is.na(CLCREAT))

for (i in unique(population_df$ID)) {
  patient_data_for_map <- df_patients %>%
    filter(ID == i)
  
  results2[[i]] <- list()
  
  # Iterate over each model in models_vancomycin
  for (model_name in names(models_vancomycin)) {
    model <- models_vancomycin[[model_name]]
    
    sim_results <- posologyr::poso_estim_map(patient_data_for_map, model)
    
    results2[[i]][[model_name]] <- sim_results$model %>% select( time, Cc, AUC) %>%
      left_join(patient_data_for_map %>% select(ID,TIME, DV),
                by = c("time" = "TIME")) %>% 
      fill(ID)
    
  }
}

combined_results2 <- list()
for (model_name in names(models_vancomycin)) {
  combined_model_data <- do.call(rbind, lapply(results2, `[[`, model_name))
  combined_results2[[model_name]] <- combined_model_data
}

results_metrics2 <- lapply(combined_results2, function(model_data) {
  model_data <- model_data %>% rename(predicted = Cc, observed = DV)
  model_data <- model_data %>% filter(!is.na(predicted) & !is.na(observed))
  
  # Calculate metrics with bootstrapping (set n_bootstrap > 0 to enable)
  calculate_metrics(model_data, n_bootstrap = 1000, conf_level = 0.95)
})

# Convert results to a data frame for plotting or further analysis

# The performance of the models was considered clinically acceptable if the rBias 
# was between −20% and 20%, with the 95% confidence intervals (CI) including zero.
# Additionally, the precision metric (rRMSE) should be as small as possible.

results_df2 <- bind_rows(
  lapply(names(results_metrics2), function(model_name) {
    metrics <- results_metrics2[[model_name]]
    
    data.frame(
      Model = model_name,
      rBias = metrics$rBias,
      rBias_Lower_CI = metrics[["rBias_CI"]][["lower.rBias"]],
      rBias_Upper_CI =  metrics[["rBias_CI"]][["upper.rBias"]],
      rRMSE = metrics$rRMSE,
      rRMSE_Lower_CI = metrics[["rRMSE_CI"]][["lower.rRMSE"]],
      rRMSE_Upper_CI =  metrics[["rRMSE_CI"]][["upper.rRMSE"]],
      MPE = metrics$MPE,
      MPE_Lower_CI =  metrics[["MPE_CI"]][["lower.MPE"]],
      MPE_Upper_CI =  metrics[["MPE_CI"]][["upper.MPE"]],
      MAPE = metrics$MAPE,
      MAPE_Lower_CI =  metrics[["MAPE_CI"]][["lower.MAPE"]],
      MAPE_Upper_CI =  metrics[["MAPE_CI"]][["upper.MAPE"]]
    )
  })
)

df_goti2 <- combined_results2$Goti2018 %>% 
  mutate(
    Cc = round(Cc)
  )

# Plot rBias for each model
plot_rBias2 <- ggplot(results_df2, aes(x = Model, y = rBias, fill = Model)) +
  geom_point(show.legend = FALSE) +
  geom_errorbar(aes(ymin = rBias_Lower_CI, ymax = rBias_Upper_CI), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 30, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = -30, linetype = "dotted", color = "blue") +
  labs(title = "Relative Bias for Each Model",
       x = "Model",
       y = "rBias (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

df_goti <- combined_results$Goti2018 %>% 
  filter(!is.na(DV)) %>%
  mutate(
    Cc = round(Cc)
  )

ggplot(df_goti, aes(x = DV, y = Cc)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 10) + 
  labs(title = "Model Fit for Goti2018 Model",
       x = "Observed",
       y = "Predicted") +
  coord_cartesian(xlim = c(0, max(df_goti$DV, df_goti$Cc)), 
                  ylim = c(0, max(df_goti$DV, df_goti$Cc))) +
  theme_classic()

# Plot combinato ----------------------------------------------------------

df_combinato <-
  bind_rows(results_df %>% mutate(method = "Population"),
            results_df2 %>% mutate(method = "MAP"))

library(patchwork)

plot_rBias3 <- ggplot(df_combinato, aes(x = Model, y = rBias, fill = method)) +
  geom_point(aes(color = method), position = position_dodge(width = 0.5), width = 1, show.legend = TRUE) +
  geom_errorbar(aes(ymin = rBias_Lower_CI, ymax = rBias_Upper_CI, color = method),
                position = position_dodge(width = 0.5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 30, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = -30, linetype = "dotted", color = "blue") +
  labs(title = "",
       x = "",
       y = "rBias (%)",
       fill = "Method",
       color = "Method") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0,vjust = 0.5),
        axis.title.y = element_text(size = 18),
        legend.position = "bottom") 

ggsave("plotBIAS.png", plot = plot_rBias3, width = 10, height = 6, dpi = 300)

plotMAPE <- ggplot(df_combinato, aes(x = Model, y = rRMSE, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), show.legend = TRUE) +
  geom_errorbar(aes(ymin = rRMSE_Lower_CI, ymax = rRMSE_Upper_CI, color = method),
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "",
       x = "",
       y = "rRMSE (%)",
       fill = "Method",
       color = "Method") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0, vjust = 0.5),
        axis.title.y = element_text(size = 18),
        legend.position = "bottom")

ggsave("plotrRMSE.png", plot = plotMAPE, width = 10, height = 6, dpi = 300)


#COmbined plots
plot_rBias3 <- ggplot(df_combinato, aes(x = Model, y = rBias, color = method)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = rBias_Lower_CI, ymax = rBias_Upper_CI),
                position = position_dodge(width = 0.5), width = 0.4, size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 30, linetype = "dotted", color = "blue") +
  geom_hline(yintercept = -30, linetype = "dotted", color = "blue") +
  labs(x = "", y = "rBias (%)") +  # Removed color label from here
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "none")  # Turn off individual legend

plotMAPE <- ggplot(df_combinato, aes(x = Model, y = rRMSE, fill = method, color = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.8) +
  geom_errorbar(aes(ymin = rRMSE_Lower_CI, ymax = rRMSE_Upper_CI),
                position = position_dodge(width = 0.5), width = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "", y = "rRMSE (%)") +  # Removed fill and color labels from here
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        legend.position = "bottom")  # Turn off individual legend

# Combine plots and explicitly collect legends
combined_plot <- (plot_rBias3 + plotMAPE) +
  plot_layout(guides = "collect") 

ggsave("combined.png", plot = combined_plot, width = 10, height = 6, dpi = 300)

# Plot AUC ----------------------------------------------------------------

dv_valle <- combined_results$Goti2018 %>% 
  filter(!is.na(DV)) %>%
  filter(time == 47.5) %>% 
  select(DV, ID)

df_goti <- combined_results$Goti2018 %>% 
  left_join(df_goti2 %>% select(time, ID, AUC_map = AUC), by = c("time", "ID")) %>% 
  filter(time %% 24 == 0) %>%
  group_by(ID) %>%
  mutate(
    AUC = round(AUC),
    AUC24_pop = round(AUC - lag(AUC, default = 0)),
    AUC24_map = round(AUC_map - lag(AUC_map, default = 0)),
    Cc = round(Cc)
  ) %>% 
  filter(
    time == 48
  ) %>% 
  ungroup() %>% 
  select(- DV) %>% 
  left_join(
    dv_valle, by = "ID"
  ) %>% 
  filter(!is.na(DV))

ggplot(df_goti, aes(x = AUC_map, y = AUC)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 10) + 
  labs(title = "Model Fit for Goti2018 Model",
       x = "Observed",
       y = "Predicted") +
  coord_cartesian(xlim = c(0, max(df_goti$AUC_map, df_goti$AUC)), 
                  ylim = c(0, max(df_goti$AUC_map, df_goti$AUC))) +
  theme_classic()

p1 <- ggplot(df_goti, aes(x = DV, y = Cc)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_text_repel(aes(label = ID), size = 3, max.overlaps = 10) + 
  labs(title = "Model Fit for Goti2018 Model",
       x = "Observed",
       y = "Predicted") +
  coord_cartesian(xlim = c(0, max(df_goti$DV, df_goti$Cc)), 
                  ylim = c(0, max(df_goti$DV, df_goti$Cc))) +
  theme_classic()


# improvment --------------------------------------------------------------


improvment <- df_combinato %>% 
  select(method,Model, rBias, rRMSE) %>% 
  pivot_wider(names_from = method, values_from = c(rBias, rRMSE)) %>%
  mutate(
    # Calculate improvement as difference (MAP - Population) for rBias and rRMSE
    rBias_Improvement = round((rBias_Population - rBias_MAP)/rBias_Population * 100),
    rRMSE_Improvement = round((rRMSE_Population - rRMSE_MAP)/rRMSE_Population * 100)
  )


# AUC plot ----------------------------------------------------------------

combined_auc1 <- bind_rows(
  lapply(names(combined_results), function(model_name) {
    combined_results[[model_name]] %>%
      mutate(model_name = model_name,
             method = "POP") %>% 
      filter(time == 48|time ==24) %>% 
      mutate(
        AUC = AUC- lag(AUC)
      ) %>% 
      filter(time == 48)
  })
)

combined_auc2 <- bind_rows(
  lapply(names(combined_results2), function(model_name) {
    combined_results2[[model_name]] %>%
      mutate(model_name = model_name,
             method = "MAP")  %>% 
      filter(time == 48|time == 24) %>% 
      mutate(
        AUC = AUC- lag(AUC)
      ) %>% 
      filter(time == 48)
  })
)

df_tot <- bind_rows(combined_auc1, combined_auc2) %>% 
  filter(AUC != max(AUC, na.rm = TRUE))

df_tot_mean <- df_tot %>% 
  group_by(model_name, method) %>% 
  summarise(AUC = median(AUC))

auc_plot <- ggplot(df_tot, aes(x = model_name, y = AUC, fill = method)) +
  geom_boxplot() +
  geom_hline(yintercept = 400, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 600, linetype = "dashed", color = "black") +
  labs(
    title = "",
    x = "",
    y = expression("AUC"[0-24])  # Creates AUC with 24-48 as subscript
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(size = 22, hjust = 0, vjust = 0.5),
        axis.title.y = element_text(size = 18),
        legend.position = "none")

ggsave("auc_plot.png", plot = auc_plot, width = 10, height = 6, dpi = 300)