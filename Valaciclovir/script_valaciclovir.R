library(posologyr);library(tidyverse)

Angelini2026 =  list(
  ppk_model = rxode2::rxode({
    cmt(depot);
    cmt(centr);
    
    depot(0) = 0;
    centr(0) = 0;
    AUC(0)   = 0;
    
    TVka = THETA_ka;
    TVCl = THETA_Cl * (CLCR / 124)^0.44;
    TVV  = THETA_V;
    
    ka = TVka * exp(ETA_ka);
    Cl = TVCl * exp(ETA_Cl);
    V  = TVV  * exp(ETA_V);
    
    d/dt(depot) = - ka * depot;
    d/dt(centr) =   ka * depot - Cl / V * centr;
    
    Cc = centr / V;
    
    d/dt(AUC) = Cc;
  }),
  error_model = function(f, sigma) {
    g <- sigma[1] + sigma[2] * f
    return(g)
  },
  theta = c(
    THETA_ka   = 1.73,
    THETA_Cl   = 52.98,
    THETA_V    = 187.79
  ),
  omega = lotri::lotri({
    ETA_ka + ETA_Cl + ETA_V ~
      c(0.31^2, #monolix output is SD, not variance
        0.00, 0.15^2,
        0.00, 0.00, 0.16^2)
  }),
  covariates = c("CLCR"),
  sigma = c(additive_a = 0.19, proportional_b = 0.26)
)


regimens <- tibble(
  regimen = c("4 g q12h", "2 g q6h"),
  dose    = c(4000, 2000),
  tau     = c(12, 6)
)

# ---- Simulation settings ----------------------------------------------------
set.seed(2026)

n_sim     <- 1000
crcl_grid <- seq(30, 150, by = 10)
sim_days  <- 7                 # long enough to reach SS even at CrCl 20
t_end     <- sim_days * 24
t_start   <- t_end - 24        # last 24 h = steady-state dosing day

auc_targets <- c(30*4)  # mg·h/L, efficacy thresholds


# ---- One regimen x one CrCl -------------------------------------------------
simulate_one <- function(dose, tau, crcl) {
  ev <- rxode2::et(amt = dose, ii = tau, addl = t_end / tau - 1, cmt = "depot") %>%
    rxode2::et(seq(t_start, t_end, by = 0.25)) %>%
    rxode2::et(id = 1:n_sim) %>%
    as.data.frame() %>%
    mutate(CLCR = crcl)
  
  rxode2::rxSolve(Angelini2026$ppk_model,
                  params     = Angelini2026$theta,
                  omega      = Angelini2026$omega,
                  events     = ev,
                  returnType = "data.frame") %>%
    group_by(id) %>%
    summarise(
      AUC24 = AUC[which.max(time)] - AUC[which.min(time)],
      Cmin  = min(Cc),
      Cmax  = max(Cc),
      .groups = "drop"
    )
}

# ---- Run all scenarios ------------------------------------------------------
sim_res <- expand_grid(regimens, CLCR = crcl_grid) %>%
  mutate(res = pmap(list(dose, tau, CLCR), simulate_one)) %>%
  unnest(res)

# ---- PTA --------------------------------------------------------------------
z <- qnorm(0.975)
#Wilson score 95% confidence interval for a proportion

pta <- sim_res %>%
  crossing(target = auc_targets) %>%
  group_by(regimen, CLCR, target) %>%
  summarise(
    n       = n(),
    x       = sum(AUC24 >= target),
    .groups = "drop"
  ) %>%
  mutate(
    PTA    = x / n,
    centre = (PTA + z^2 / (2 * n)) / (1 + z^2 / n),
    half   = z * sqrt(PTA * (1 - PTA) / n + z^2 / (4 * n^2)) / (1 + z^2 / n),
    lower  = centre - half,
    upper  = centre + half
  ) %>%
  select(regimen, CLCR, target, n, PTA, lower, upper)

pta %>%
  mutate(target = paste0("AUC24,ss \u2265 ", target, " mg\u00b7h/L")) %>%
  ggplot(aes(CLCR, PTA, colour = regimen)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = regimen),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point() +
  geom_hline(yintercept = 0.9, linetype = "dashed", colour = "grey50") +
  facet_wrap(~ target) +
  scale_x_continuous(breaks = crcl_grid) +
  scale_y_continuous(labels = scales::percent, limits = c(0.75, 1)) +
  labs(x = "CrCl (mL/min)", y = "Probability of target attainment",
       colour = NULL, linetype = NULL, fill = NULL) +
  theme_bw() +
  theme(legend.position = "bottom") #+
  #facet_wrap(~ regimen)

# ---- Cmin vs AUC24,ss by regimen ---------------------------------------------
fit_lines <- sim_res %>%
  group_by(regimen) %>%
  summarise(fit = list(lm(log(AUC24) ~ log(Cmin))),
            cmin_rng = list(range(Cmin)), .groups = "drop") %>%
  mutate(pred = map2(fit, cmin_rng, ~ {
    nd <- tibble(Cmin = seq(.y[1], .y[2], length.out = 200))
    nd %>% mutate(AUC24 = exp(predict(.x, newdata = nd)))
  })) %>%
  select(regimen, pred) %>%
  unnest(pred)

ggplot(sim_res, aes(Cmin, AUC24)) +
  geom_point(aes(colour = CLCR), alpha = 0.15, size = 0.6) +
  geom_line(data = fit_lines, colour = "black", linewidth = 0.8) +
  geom_hline(yintercept = auc_targets, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ regimen) +
  scale_colour_viridis_c(name = "CrCl (mL/min)") +
  labs(x = "Cmin,ss (mg/L)", y = "AUC24,ss (mg\u00b7h/L)") +
  theme_bw() +
  theme(legend.position = "bottom")

# ---- Cmin equivalent to the AUC target, by regimen ---------------------------
sim_res %>%
  group_by(regimen) %>%
  summarise(fit = list(lm(log(AUC24) ~ log(Cmin))), .groups = "drop") %>%
  mutate(
    intercept = map_dbl(fit, ~ coef(.x)[1]),
    slope     = map_dbl(fit, ~ coef(.x)[2]),
    r2        = map_dbl(fit, ~ summary(.x)$r.squared)
  ) %>%
  crossing(target = auc_targets) %>%
  mutate(Cmin_equiv = exp((log(target) - intercept) / slope)) %>%
  select(regimen, target, Cmin_equiv, r2)

# ---- One log-log fit per regimen x CrCl --------------------------------------
fit_lines <- sim_res %>%
  group_by(regimen, CLCR) %>%
  summarise(fit      = list(lm(log(AUC24) ~ log(Cmin))),
            cmin_rng = list(range(Cmin)),
            .groups  = "drop") %>%
  mutate(pred = map2(fit, cmin_rng, ~ {
    nd <- tibble(Cmin = seq(.y[1], .y[2], length.out = 200))
    nd %>% mutate(AUC24 = exp(predict(.x, newdata = nd)))
  })) %>%
  select(regimen, CLCR, pred) %>%
  unnest(pred)

# ---- Plot ---------------------------------------------------------------------
library(ggrepel)

# Label positions: right end of each line
line_labels <- fit_lines %>%
  group_by(regimen, CLCR) %>%
  slice_max(Cmin, n = 1) %>%
  ungroup()

x_max <- max(fit_lines$Cmin)

ggplot(fit_lines, aes(Cmin, AUC24, colour = CLCR, group = CLCR)) +
  geom_hline(yintercept = auc_targets, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 12, y = auc_targets,
           label = paste0("Target: ", auc_targets, " mg\u00b7h/L"),
           hjust = 0, vjust = -0.5, size = 3.2, colour = "grey30") +
  geom_line(linewidth = 0.9) +
  geom_text_repel(data = line_labels,
                  aes(label = CLCR),
                  hjust = 0, nudge_x = 0.07 * x_max,
                  direction = "y", segment.size = 0.2,
                  segment.colour = "grey60", size = 3,
                  min.segment.length = 0, show.legend = FALSE) +
  facet_wrap(~ regimen) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  scale_colour_viridis_c(end = 0.85, guide = "none") +
  labs(
    title    = "Relationship between trough concentration and daily exposure",
    subtitle = "Log-log regression per CrCl value (labels: CrCl, mL/min)",
    x = expression(C[min*","*ss]~"(mg/L)"),
    y = expression(AUC[24*","*ss]~"(mg"%.%"h/L)")
  ) +
  theme_bw(base_size = 12) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.2, "lines"),
    plot.title       = element_text(face = "bold"),
    plot.subtitle    = element_text(colour = "grey35")
  )


cmin_grid <- c(0.5, 1, 2, 3, 5, 7.5, 10)

cmin_grid <- c(1, 1.5, 2, 3, 4, 5, 6)

nomogram <- sim_res %>%
  group_by(regimen, CLCR) %>%
  summarise(fit = list(lm(log(AUC24) ~ log(Cmin))),
            lo  = quantile(Cmin, 0.01),
            hi  = quantile(Cmin, 0.99),
            .groups = "drop") %>%
  crossing(Cmin = cmin_grid) %>%
  mutate(
    AUC24_pred = map2_dbl(fit, Cmin, ~ exp(predict(.x, newdata = tibble(Cmin = .y)))),
    AUC24_pred = if_else(Cmin < lo | Cmin > hi, NA_real_, round(AUC24_pred))
  ) %>%
  select(regimen, CrCl = CLCR, Cmin, AUC24_pred) %>%
  pivot_wider(names_from = Cmin, values_from = AUC24_pred,
              names_glue = "Cmin {Cmin}")

nomogram %>% filter(regimen == "4 g q12h")
nomogram %>% filter(regimen == "2 g q6h")