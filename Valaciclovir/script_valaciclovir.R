library(posologyr);library(tidyverse)

Angelini2026 =  list(
  ppk_model = rxode2::rxode({
    cmt(depot);
    cmt(centr);
    
    depot(0) = 0;
    centr(0) = 0;
    AUC(0)   = 0;
    
    TVka = THETA_ka;
    TVCl = THETA_Cl * (CLCR / 124)^0.56;
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
    THETA_ka   = 1.94,
    THETA_Cl = 53.17,
    THETA_V    = 185.53
  ),
  omega = lotri::lotri({
    ETA_ka + ETA_Cl + ETA_V ~
      c(0.26,
        0.00, 0.13,
        0.00, 0.03, 0.14)
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
crcl_grid <- seq(20, 120, by = 10)
sim_days  <- 7                 # long enough to reach SS even at CrCl 20
t_end     <- sim_days * 24
t_start   <- t_end - 24        # last 24 h = steady-state dosing day

auc_targets <- c(28*4)  # mg·h/L, efficacy thresholds (placeholders)


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
  ggplot(aes(CLCR, PTA, colour = regimen, linetype = regimen)) +
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

