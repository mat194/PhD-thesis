library(posologyr); library(tidyverse); library(readxl)
library(ggrepel)  ; library(flextable); library(gtsummary)
library(janitor)

source("~/PhD-thesis/Vancomycin/models.R")

# Load data ---------------------------------------------------------------

df1 <- read_excel("Vancomycin/Data/data_bloomy.xlsx") %>% filter(!is.na(LD)) |> mutate(centro = "Verona")
df_init <- read_excel("Vancomycin/Data/udine.xlsx") %>% 
  janitor::clean_names() %>% 
  dplyr::filter(!is.na(record_id)) %>% 
  ungroup() %>% 
  transmute(
    record_id, 
    age = floor(age_y),
    sex = tolower(sex_male_female),
    TIME = time, 
    EVID = evid_0_osservazione_1_dose,
    AMT = amt,
    DUR = dur,
    DV = conc,
    altezza = altezza_m,
    peso =  peso_kg,
    vm = 0,
    CRRT = crrt, 
    LD = ld,
    dialisi, 
    creatinina ,
    ii_interdose_interval,
    addl,
    centro = "Udine"
  ) |> 
  bind_rows(df1) |>
  group_by(centro, record_id) |>
  mutate(
    ID = cur_group_id()
  ) |>
  ungroup()

dose_expanded  <- df_init %>% 
  dplyr::filter(EVID == 1) %>% 
  mutate(
    ADDL = replace_na(addl, 0L),
    II   = replace_na(ii_interdose_interval, 0)
  ) %>% 
  #distinct(record_id, TIME, .keep_all = TRUE) %>%
  rowwise() %>%
  mutate(
    TIME_seq = list(seq(
      from = TIME,
      by = as.difftime(II, units = "hours"),
      length.out = ADDL + 1
    ))
  ) %>%
  ungroup() %>%
  select(-TIME) %>%
  unnest(TIME_seq) %>%
  rename(TIME = TIME_seq) %>%
  mutate(
    EVID = 1,
    DV = NA,
    TIME = as.numeric(TIME)
  ) %>% 
  arrange(record_id, TIME) %>% 
  tidyr::fill(creatinina, .direction = "updown")

df <- dplyr::bind_rows(df_init %>% filter(EVID == 0), dose_expanded) %>%
  arrange(ID, TIME) |> 
  filter_out(is.na(sex))

na_sex <- df %>% filter(is.na(sex)) %>% pull(record_id)
na_creat <- df %>% filter(is.na(creatinina)) %>% pull(record_id)

to_fix <- df %>% filter(record_id %in% na_sex)
to_fix2 <- df %>% filter(record_id %in% na_creat)

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


calculate_crcl <- function(age, weight_kg, serum_creatinine_mg_dl, gender) {
  CrCl <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
  if (gender == "female") {
    CrCl <- CrCl * 0.85
  }
  return(CrCl)
}


df_patients <- df %>% 
  rowwise() %>%  
  mutate(
    CLCREAT = if_else(
      is.na(creatinina) | creatinina == 0, NA,
      calculate_crcl(
        age = age,  
        weight_kg = peso,
        serum_creatinine_mg_dl = creatinina,
        gender = sex
      )
    )
  ) %>%
  rename(
    WT = peso,
    DIAL = dialisi
  )


# Statistiche -------------------------------------------------------------

ci <- df_patients %>% 
  filter(DUR == 24) %>% 
  select(record_id) %>% 
  distinct() %>% 
  pull()


# Population simulation ---------------------------------------------------
doses <- df_patients %>% 
  filter(!is.na(CLCREAT)) %>% 
  filter(EVID == 1) %>% 
  select(
    ID, TIME, AMT, DUR
  ) 
# farlo su tutte le osservazioni
population_df <- df_patients %>% 
  group_by(ID) %>% 
  filter(
    TIME == 0
  ) %>% 
  mutate(
    AMT = 0,
    EVID = 0,
    DUR = NA
  ) %>% 
  ungroup() %>% 
  dplyr::filter(!is.na(creatinina),!is.na(WT))

results <- list()
max_time <- max(df_patients$TIME, na.rm = TRUE)

models_to_use <- setdiff(names(models_vancomycin), "Medellin2017")

for (i in unique(population_df$ID)) {
  message("Processing item ", i, " of ", length(unique(population_df$ID)))
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
  for (model_name in models_to_use) {
    message("Testing model ", model_name)
    model <- models_vancomycin[[model_name]]
    
    sim_results <- posologyr::poso_simu_pop(
      patient_data_for_pop, 
      model, 
      n_simul = 0
    )
    
    sim_results$model$time <- seq(0, max_time, b = 0.1)
    
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

# save(results, file = "udine.Rdata")

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

for (i in unique(population_df$ID)) {
  patient_data_for_map <- df_patients %>%
    filter(ID == i)
  
  results2[[i]] <- list()
  
  # Iterate over each model in models_vancomycin
  for (model_name in models_to_use) {
    model <- models_vancomycin[[model_name]]
    
    sim_results <- posologyr::poso_estim_map(patient_data_for_map, model)
    
    results2[[i]][[model_name]] <- sim_results$model %>% select( time, Cc, AUC) %>%
      left_join(patient_data_for_map %>% select(ID,TIME, DV),
                by = c("time" = "TIME")) %>% 
      fill(ID)
    
  }
}

#save(results2, file = "udine_map.Rdata" )

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


#Combined plots
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
  filter(!is.na(DV))


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
  coord_fixed() +
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