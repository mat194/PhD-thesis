Dalbavancin <- list(
  Gregoire2025 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      TVVc <- THETA_Vc * (ABW / 66.5) ^ 0.6
      TVCl <- THETA_Cl * (AGFR / 88.5) ^ 0.21
      Vc <- TVVc * exp(ETA_Vc)
      Cl <- TVCl * exp(ETA_Cl)
      Vp <- THETA_Vp
      Q <- THETA_Q
      ke <- Cl / Vc
      k12 <- Q / Vc
      k21 <- Q / Vp
      Cc <- centr / Vc
      Cp <- periph / Vp
      d / dt(centr) <- -ke * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_Cl = 0.04,
      THETA_Vc = 5.62,
      THETA_Vp = 8.04,
      THETA_Q = 0.026
    ),
    covariates = c("ABW", "AGFR"),
    omega = lotri::lotri({
      ETA_Cl + ETA_Vc ~
        c(4.096472,
          0, 3.61458)
    }),
    sigma = c(additive_a = 1.67, proportional_b = 0.22)
  ),

  Cojutti2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      TVCl <- THETA_Cl * (eGFR / 93)^0.0043
      V1   <- THETA_V1 * exp(ETA_V1)
      V2   <- THETA_V2 * exp(ETA_V2)
      Q    <- THETA_Q  * exp(ETA_Q)
      Cl   <- TVCl     * exp(ETA_Cl)
      ke   <- Cl / V1
      k12  <- Q / V1
      k21  <- Q / V2
      Cc   <- centr / V1
      Cp   <- periph / V2
      d / dt(centr) <- -ke * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_Cl = 0.029,
      THETA_V1 = 6.14,
      THETA_V2 = 9.52,
      THETA_Q  = 0.026
    ),
    covariates = c("eGFR"),
    omega = lotri::lotri({ # reported as CV
      ETA_Cl + ETA_V1 + ETA_V2 + ETA_Q ~
        c(log((31.76 / 100) ^ 2 + 1),
          log((16.10 / 100) ^ 2 + 1), 0,
          log((37.19 / 100) ^ 2 + 1), 0, 0,
          log((45.06 / 100) ^ 2 + 1), 0, 0, 0)
    }),
    sigma = c(additive_a = 0, proportional_b = 0.3392)
  ),

  Carrothers2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph1(0) <- 0
      periph2(0) <- 0

      Cl <- THETA_Cl * (ALBUMIN / 3.7) ^ THETA_CL_ALB *
        (CLCR / 100) ^ THETA_CL_CLCR *
        (WT / 85.5) ^ THETA_CL_WT * exp(ETA_Cl)

      V1 <- THETA_V1 * (ALBUMIN / 3.7) ^ THETA_V1_ALB *
        (WT / 85.5) ^ THETA_V1_WT * exp(ETA_V1)

      V2 <- THETA_V2 * (AGE / 47) ^ THETA_V2_AGE *
        (ALBUMIN / 3.7) ^ THETA_V2_ALB *
        (WT / 85.5) ^ THETA_V2_WT * exp(ETA_V2)

      V3 <- THETA_V3 * (ALBUMIN / 3.7) ^ THETA_V3_ALB *
        (WT / 85.5) ^ THETA_V3_WT * exp(ETA_V3)

      Q2 <- THETA_Q2
      Q3 <- THETA_Q3

      k10 <- Cl / V1
      k12 <- Q2 / V1
      k21 <- Q2 / V2
      k13 <- Q3 / V1
      k31 <- Q3 / V3

      Cc <- centr / V1
      d / dt(centr)   <- -k10 * centr - k12 * centr + k21 * periph1 - k13 * centr + k31 * periph2
      d / dt(periph1) <- k12 * centr - k21 * periph1
      d / dt(periph2) <- k13 * centr - k31 * periph2
      d / dt(AUC)     <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_Cl      = 0.0531,
      THETA_V1      = 3.04,
      THETA_V2      = 8.78,
      THETA_V3      = 3.28,
      THETA_Q2      = 0.288,
      THETA_Q3      = 2.11,
      THETA_CL_ALB  = -0.477,
      THETA_CL_CLCR = 0.273,
      THETA_CL_WT   = 0.391,
      THETA_V1_ALB  = -0.340,
      THETA_V1_WT   = 0.683,
      THETA_V2_AGE  = 0.486,
      THETA_V2_ALB  = -0.413,
      THETA_V2_WT   = 0.365,
      THETA_V3_ALB  = -0.551,
      THETA_V3_WT   = 0.518
    ),
    covariates = c("ALBUMIN", "CLCR", "WT", "AGE"),
    omega = lotri::lotri({ # used the reported CV
      ETA_Cl + ETA_V2 + ETA_V1 + ETA_V3 ~
        c(log((22 / 100) ^ 2 + 1),
          0, log((41 / 100) ^ 2 + 1),
          0, 0, log((24 / 100) ^ 2 + 1),
          0, 0, 0, log((74 / 100) ^ 2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.0362)
  ),

  Cojutti2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- THETA_CL * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q  * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d/dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) <-  k12 * centr - k21 * periph
      d/dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g = sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 0.106,
      THETA_V1 = 17.40,
      THETA_V2 = 15.10,
      THETA_Q  = 0.103
    ),
    omega = lotri::lotri({ # reported as CV
      ETA_CL + ETA_V1 + ETA_V2 + ETA_Q ~
        c(
          log((36.21 / 100)^2 + 1),
          0, log((44.27 / 100)^2 + 1),
          0, 0, log((62.34 / 100)^2 + 1),
          0, 0, 0, log((35 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.18)
  ),

  Lodise2026 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)     <- 0
      centr(0)   <- 0
      periph1(0) <- 0
      periph2(0) <- 0

      CL <- THETA_CL * (CLCR / 101.26)^THETA_CL_CRCL * exp(ETA_CL)

      V1 <- THETA_V1 * (WT / 83.8)^THETA_V1_WT * exp(ETA_V1)
      V2 <- THETA_V2 * (WT / 83.8)^THETA_V2_WT * (ALBUMIN / 2.8)^THETA_V2_ALB * exp(ETA_V2)
      V3 <- THETA_V3 * (AGE / 56)^THETA_V3_AGE * (WT / 83.8)^THETA_V3_WT * exp(ETA_V3)

      Q2 <- THETA_Q2
      Q3 <- THETA_Q3

      TV_A    <- THETA_A * (ALBUMIN / 2.8)^THETA_A_ALB
      logit_A <- log(TV_A / (1 - TV_A)) + ETA_A
      A       <- 1 / (1 + exp(-logit_A))

      k10 <- CL / V1
      k12 <- Q2 / V1
      k21 <- Q2 / V2
      k13 <- Q3 / V1
      k31 <- Q3 / V3

      Cc      <- centr / V1          # total concentration
      Cc_free <- A * Cc^THETA_K      # unbound concentration

      d / dt(centr)   <- -k10 * centr - k12 * centr + k21 * periph1 - k13 * centr + k31 * periph2
      d / dt(periph1) <-  k12 * centr - k21 * periph1
      d / dt(periph2) <-  k13 * centr - k31 * periph2
      d / dt(AUC)     <- Cc_free
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL      = 0.0658,
      THETA_CL_CRCL = 0.214,
      THETA_V1      = 5.67,
      THETA_V1_WT   = 0.57,
      THETA_V2      = 8.91,
      THETA_V2_WT   = 0.82,
      THETA_V2_ALB  = -0.806,
      THETA_V3      = 11.1,
      THETA_V3_AGE  = 0.628,
      THETA_V3_WT   = 0.559,
      THETA_Q2      = 0.0259,
      THETA_Q3      = 0.921,
      THETA_A       = 0.00136,
      THETA_A_ALB   = -0.782,
      THETA_K       = 1.32
    ),
    covariates = c("CLCR", "WT", "ALBUMIN", "AGE"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V3 + ETA_V2 + ETA_A ~
        c(
          0.0508,
          0.0264,  0.0386,
          0.0416,  0.0478,  0.0861,
          0,       0,       0,    0.0897,
          0,       0,       0,    0,  0.106
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.13),
    auc_is_free = TRUE
  ),

  Baiardi2025 = list(
    ppk_model = rxode2::rxode2({
      centr(0) <- 0
      periph(0) <- 0
      AUC(0) <- 0

      CL <- THETA_CL * (WT / 70)^0.75 * exp(ETA_CL)
      V1 <- THETA_V1 * (WT / 70) * exp(ETA_V1)
      V2 <- THETA_V2 * (WT / 70) * exp(ETA_V2)
      Q  <- THETA_Q  * (WT / 70)^0.75 * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d/dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) <-  k12 * centr - k21 * periph
      d/dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g = sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 0.053,
      THETA_V1 = 3.04,
      THETA_V2 = 8.78,
      THETA_Q  = 0.288
    ),
    covariates = c("WT"),
    omega = lotri::lotri({ # reported as CV
      ETA_CL + ETA_V1 + ETA_V2 + ETA_Q ~
        c(
          log((22 / 100)^2 + 1),
          0, log((24 / 100)^2 + 1),
          0, 0, log((41 / 100)^2 + 1),
          0, 0, 0, log((30 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.0362)
  ),

  Banavent2025 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      V1 <- THETA_V1 * exp(ETA_V1)
      CL <- THETA_CL * exp(ETA_CL)

      k10 <- CL / V1
      Cc  <- centr / V1

      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 0.036,
      THETA_V1 = 17.9
    ),
    omega = lotri::lotri({ # reported as omega^2
      ETA_CL + ETA_V1 ~
        c(
          0.200,
          0, 0.290
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.120)
  ),

  Pai2026 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      periph(0) <- 0

      # CL modelled on log-transformed BMI/25 and absolute eGFR/75 (CKD-EPI × BSA_Mosteller / 1.73)
      CL <- THETA_CL * (BMI / 25)^THETA_CL_BMI * (AGFR_M / 75)^THETA_CL_AGFR * exp(ETA_CL)
      V1 <- THETA_V1                  # IIV on V1 removed in final model
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q  * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2

      Cc <- centr / V1

      d/dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) <-  k12 * centr - k21 * periph
      d/dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL      = 0.038,
      THETA_CL_BMI  = 0.22,
      THETA_CL_AGFR = 0.25,
      THETA_V1      = 6.59,
      THETA_Q       = 0.037,
      THETA_V2      = 13.83
    ),
    covariates = c("BMI", "AGFR_M"),
    omega = lotri::lotri({ # reported as omega SD; stored here as omega^2
      ETA_CL + ETA_Q + ETA_V2 ~
        c(0.19^2,
          0, 0.61^2,
          0, 0, 0.60^2)
    }),
    sigma = c(additive_a = 0, proportional_b = 0.34)
  ),

  Chiriac2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      V1   <- THETA_V1 * exp(ETA_V1)
      V2   <- THETA_V2 * exp(ETA_V2)
      Q    <- THETA_Q
      Cl   <- THETA_Cl * exp(ETA_Cl)
      ke   <- Cl / V1
      k12  <- Q / V1
      k21  <- Q / V2
      Cc   <- centr / V1
      Cp   <- periph / V2
      d / dt(centr)  <- -ke * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_Cl = 0.050,
      THETA_V1 = 6.5,
      THETA_V2 = 15.4,
      THETA_Q  = 0.476
    ),
    omega = lotri::lotri({ # reported as omega/SD
      ETA_Cl + ETA_V1 + ETA_V2 ~
        c(0.230^2,
          0, 0.260^2,
          0, 0, 0.410^2)
    }),
    sigma = c(additive_a = 0, proportional_b = 0.1)
  ),
  Sayadi2026 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      periph(0) <- 0
      # Covariates centred on population medians: eGFR 75 mL/min/1.73m², WT 79 kg
      # WT exponents on V1 and V2 fixed to 1.0 (standard allometric scaling)
      CL <- THETA_CL * (eGFR / 75)^THETA_CL_eGFR * (WT / 79)^THETA_CL_WT * exp(ETA_CL)
      V1 <- THETA_V1 * (WT / 79) * exp(ETA_V1)
      V2 <- THETA_V2 * (WT / 79) * exp(ETA_V2)
      Q  <- THETA_Q  * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2

      Cc <- centr / V1

      d/dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) <-  k12 * centr - k21 * periph
      d/dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL      = 0.035,   
      THETA_CL_eGFR = 0.22,
      THETA_CL_WT   = 0.61,
      THETA_V1      = 4.74,    
      THETA_Q       = 0.033,   
      THETA_V2      = 8.74     
    ),
    covariates = c("eGFR", "WT"),
    omega = lotri::lotri({ # reported as CV (Table 4)
      ETA_CL + ETA_V1 + ETA_Q + ETA_V2 ~
        c(log((17.5 / 100)^2 + 1),
          0, log((20.02 / 100)^2 + 1),
          0, 0, log((52.39 / 100)^2 + 1),
          0, 0, 0, log((35.22 / 100)^2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.21)
  )
)
