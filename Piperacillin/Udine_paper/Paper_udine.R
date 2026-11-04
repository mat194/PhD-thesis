library(tidyverse);library(readxl);library(data.table);library(ggrepel)
library(mediation);library(gtsummary);library(gt)
source("~/PhD-thesis/Utils/Functions.R")
set.seed(123)
# Model ------------------------------------------------------------------
paper_udine_models <- list(
  # https://doi.org/10.1128/aac.01253-25
  # Laporte-Amargos et al. 2026, AAC
  # Population PK of piperacillin in 44 hematological patients with febrile
  # neutropenia (BEATLE study); 221 plasma samples, one-compartment model
  # with first-order elimination described the data best (a 2-compartment
  # model was tested but Vd estimation was biased, likely because most
  # samples were collected during the elimination phase). Time-varying
  # creatinine clearance (Cockcroft-Gault, CLCR_CG) was the only significant
  # covariate on CL, normalized to the population median CrCL (99.3 mL/min):
  #   CL_i (L/h) = 12 * (CLCR_CG_i / 99.3)^0.64
  # BOV on CL was estimated (16.4% CV) but had a negligible additional
  # effect on BSV (median intra-subject CrCL variation between occasions
  # was only 6.6%). No covariate was retained on Vd.
  Laporte_Amargos_2026 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0

      CL <- THETA_CL * (CLCR_CG / 99.3)^THETA_CL_2 * exp(ETA_CL + KAPPA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)

      k10 <- CL / V1
      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 12.0,
      THETA_CL_2 = 0.64,
      THETA_V1   = 29.8
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((36.2 / 100)^2 + 1),
          0, log((57.3 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 5.5, proportional_b = 0.20),
    pi_matrix = lotri::lotri({ # Provided as CV, BOV on CL
      KAPPA_CL ~ log((16.4 / 100)^2 + 1)
    })
  )
)




# Population sim ----------------------------------------------------------
df <- readRDS("~/PhD-thesis/Piperacillin/Udine_paper/population_df.rds")
pop_n_simul   <- 1000     
pop_pred_stat <- "typical"
pop_grid_n    <- 400      
results <- list()
m_out <- list()

for (i in sort(unique(df$ID))) {
  message("Processing item ", i, " of ", length(unique(df$ID)))
  
  df_i <- df |> filter(ID == i)
  max_time_i <- max(df_i$TIME, na.rm = TRUE)
  
  by_i <- if (pop_n_simul > 0) max(0.1, max_time_i / pop_grid_n) else 0.1
  et_i <- rxode2::as.et(df_i)
  et_i$add.sampling(seq(0, max_time_i, by = by_i))
  
  obs_i <- df_i |> filter(EVID == 0, !is.na(DV)) |> select(TIME, DV)
  
  id_key <- as.character(i)
  results[[id_key]] <- list()
  
  for (model_name in names(paper_udine_models)) {
    message("Testing model ", model_name, " for ID ", i)
    
    result_i <- tryCatch({
      model <- paper_udine_models[[model_name]]
      
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
for (model_name in names(paper_udine_models)) {
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

combined_results_tot <- bind_rows(combined_results, .id = "model")

fit_ipred <- lm(Cc ~ DV, data = combined_results_tot |> filter_out(is.na(DV)))
r2_ipred  <- summary(fit_ipred)$r.squared

plot_ipred <- ggplot(combined_results_tot |> filter_out(is.na(DV)), aes(x = DV, y = Cc)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  geom_smooth(method = "lm", se = FALSE, colour = "blue", linewidth = 0.7) +
  annotate("text",
           x = 0.05 * max(combined_results_tot$DV, na.rm = TRUE),
           y = 0.95 * max(combined_results_tot$Cc, na.rm = TRUE),
           label = sprintf("R^2 == %.3f", r2_ipred),
           parse = TRUE, hjust = 0, size = 5) +
  labs(title = "Model Fit", x = "Observed", y = "Predicted") +
  coord_fixed(xlim = c(0, max(combined_results_tot$DV, na.rm = TRUE)),
              ylim = c(0, max(combined_results_tot$Cc, na.rm = TRUE))) +
  theme_classic()

ggsave("Piperacillin/Udine_paper/Plots/model_fit.png",
       plot_ipred, width = 7, height = 6, dpi = 300)

obs_df <- df %>%
  filter(EVID %in% c(0), !is.na(DV)) %>%
  transmute(ID, TIME, OBSERVED = DV) %>%
  distinct()

for (i in sort(unique(combined_results_tot$ID))) {
  ggplot(combined_results_tot |> filter(ID == i), aes(x = time, y = Cc, group = ID))+
    geom_ribbon(aes(ymin = Cc_p05, ymax = Cc_p95), alpha = 0.2, fill = "steelblue", colour = NA) +
    geom_line() +
    geom_point(data = obs_df |> filter(ID == i), aes(x = TIME, y = OBSERVED, group = ID), size = 1, colour  = "brown", na.rm = TRUE) +
    theme_minimal() +
    facet_wrap(~ model, scales = "free")

  ggsave(paste0("Piperacillin/Udine_paper/Plots/Concentration/POP/",i,".png"))

}



# mediation analysis ------------------------------------------------------

db <-  read_excel("Piperacillin/Udine_paper/data_udine_pip.xlsx", 
                  sheet = "dataset clinico -3pz") |> 
  transmute(
    Acute = if_else(Type_of_disease == 1, "Acute", "Chronic"),
    CRCL  = `Cl crea_TDM1`,
    Css   = Conc_first_TDM
  )
mediator_model <- lm(CRCL ~ Acute, data = db)          # mediator ~ treatment
outcome_model  <- lm(Css ~ Acute + CRCL, data = db)     # outcome ~ treatment + mediator

results <- mediate(
  mediator_model,
  outcome_model,
  treat    = "Acute",
  mediator = "CRCL",
  boot = TRUE, sims = 1000
)
summary(results)

# Manuscript table: mediator model + outcome model + mediation effects ------
term_labels <- c(
  "AcuteChronic" = "Disease acuity: Chronic vs Acute",
  "CRCL"         = "Creatinine clearance (CRCL)"
)

mediator_tidy <- broom::tidy(mediator_model, conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  mutate(Section = "Mediator model: CRCL ~ Acute",
         Term = recode(term, !!!term_labels))

outcome_tidy <- broom::tidy(outcome_model, conf.int = TRUE) |>
  filter(term != "(Intercept)") |>
  mutate(Section = "Outcome model: Css ~ Acute + CRCL",
         Term = recode(term, !!!term_labels))

med_summary <- tibble(
  Term      = c("ACME (indirect, via CRCL)", "ADE (direct)", "Total Effect", "Proportion Mediated"),
  estimate  = c(results$d.avg, results$z.avg, results$tau.coef, results$n.avg),
  conf.low  = c(results$d.avg.ci[1], results$z.avg.ci[1], results$tau.ci[1], results$n.avg.ci[1]),
  conf.high = c(results$d.avg.ci[2], results$z.avg.ci[2], results$tau.ci[2], results$n.avg.ci[2]),
  p.value   = c(results$d.avg.p, results$z.avg.p, results$tau.p, results$n.avg.p),
  Section   = "Causal mediation effects (Chronic vs Acute)"
)

manuscript_tbl_data <- bind_rows(
  mediator_tidy |> select(Section, Term, estimate, conf.low, conf.high, p.value),
  outcome_tidy  |> select(Section, Term, estimate, conf.low, conf.high, p.value),
  med_summary   |> select(Section, Term, estimate, conf.low, conf.high, p.value)
)

manuscript_tbl <- manuscript_tbl_data |>
  gt(groupname_col = "Section", rowname_col = "Term") |>
  fmt_number(columns = c(estimate, conf.low, conf.high), decimals = 2) |>
  fmt(columns = p.value, fns = function(x) style_pvalue(x, prepend_p = FALSE)) |>
  cols_merge(columns = c(conf.low, conf.high), pattern = "{1}, {2}") |>
  cols_label(estimate = "Estimate", conf.low = "95% CI", p.value = "p-value") |>
  tab_header(
    title = "Mediation of the effect of disease acuity on piperacillin Css through CRCL",
    subtitle = paste0("Nonparametric bootstrap, ", results$sims, " simulations")
  ) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) |>
  tab_source_note(source_note = "ACME: average causal mediation effect; ADE: average direct effect; Css: steady-state concentration; CRCL: creatinine clearance.")

gtsave(manuscript_tbl, "Piperacillin/Udine_paper/Plots/Mediation/manuscript_tbl.docx")