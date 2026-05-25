library(posologyr)
# Models ------------------------------------------------------------------
piperacillin_models <- list(
  # https://pubmed.ncbi.nlm.nih.gov/23908259/
  # Pipe_1
  # THETA_CL_3 * RRT---> semplification from the mean(SD) values of piperacillin (should be CLEC= Quf * Sc)
  Asín_Prieto_2013 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL + (CLCR_CG / 1.95)^THETA_CL_2 + THETA_CLEC * RRT) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2^(WT / 75.3)   
      Q  <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    # rxode2 combined1() --> assumes that the additive and proportional differences are on the SD
    # To be checked the log scale for this model
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    # This would be the combined2() -->  assumes that the additive and proportional differences are combined on a variance
    # error_model <- function(f,sigma){
    #   g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
    #   return(sqrt(g))
    # },
    theta = c(
      THETA_CL    = 2.77,   
      THETA_CL_2  = 4.55,   
      THETA_CLEC  = 0.64,   
      THETA_V1    = 20.86,  
      THETA_V2    = 21.41,  
      THETA_Q     = 20.53   
    ),
    covariates = c("CLCR_CG", "RRT", "WT"),
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((44 / 100)^2 + 1),
          log((80 / 100)^2 + 1), log((75 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.19, proportional_b = 0) 
  ),
  # https://pubmed.ncbi.nlm.nih.gov/31978581/
  # Pipe_3
  # Bue2020
  Bue_2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- THETA_CL * (WT / 70)^0.75 * exp(ETA_CL)
      V1 <- THETA_V1 * (WT / 70)
      V2 <- THETA_V2 * (WT / 70) * exp(ETA_V2)
      Q <- THETA_Q * (WT / 70)^0.75 * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 3.3,
      THETA_V1 = 6.77,
      THETA_V2 = 10.0,
      THETA_Q = 15.4
    ),
    covariates = c("WT"),
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_Q + ETA_V2 ~
        c(
          log((10.9 / 100)^2 + 1),
          0, log((54.8 / 100)^2 + 1),
          0, 0, log((26 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.166)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/33569597/
  # Pipe_4
  Fillâtre_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- THETA_CL * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 5.01,
      THETA_V1 = 17.60,
      THETA_V2 = 10.40,
      THETA_Q = 7.01
    ),
    omega = lotri::lotri({ # Provided as omega/SD
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          0.60^2,
          0, 0.33^2,
          0, 0, 0.84^2
        )
    }),
    sigma = c(additive_a = 5.15, proportional_b = 0.04)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC8694146/
  # Pipe_5
  Hahn_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (1 - THETA_CL_2 * ECMO) + (THETA_CL_3 * (CLCR_CG - 54.7))) * exp(ETA_CL)
      V1 <- THETA_V1 * (1 + THETA_V1_2 * CVVHDF) * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 9.4,
      THETA_CL_2 = 0.092,
      THETA_CL_3 = 0.115,
      THETA_V1 = 6.56,
      THETA_V1_2 = 1.46,
      THETA_V2 = 14.2,
      THETA_Q = 17.2
    ),
    covariates = c("CVVHDF", "ECMO", "CLCR_CG"),
    omega = lotri::lotri({ # Provided as omega squared/SD^2/variance/IIV
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          0.0523,
          0, 0.291,
          0, 0, 0.138
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.313)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC4068553/
  # Pipe_6
  Jeon_2014 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (CLCR_CG / 132) + ((-THETA_CL_1) * DAI)) * exp(ETA_CL)
      V1 <- (THETA_V1 + (THETA_V1_2 * SEPSIS)) * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 16.6,
      THETA_CL_1 = 0.0874,
      THETA_V1 = 25.3,
      THETA_V1_2 = 14.8,
      THETA_V2 = 16.1,
      THETA_Q = 0.636
    ),
    covariates = c("CLCR_CG", "DAI", "SEPSIS"), # DAI = day after burn injury
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_Q ~
        c(
          log((35.4 / 100)^2 + 1),
          0.434, log((42.4 / 100)^2 + 1), # The correlation between the interindividual variability for CL and V
          0, 0, log((90.3 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.359, proportional_b = 0.185)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9047688/
  # Pipe_7
  Kim_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (WT / 70)^0.75 * exp(THETA_CL_2 * (eGFR_CKD_EPI_CYSTATIN - 52.77))) * exp(ETA_CL)

      if (ECMO == 1) {
        V1 <- THETA_V1_E1 * (WT / 70) * exp(ETA_V1_E1)
      } else {
        V1 <- THETA_V1_E0 * (WT / 70) * exp(ETA_V1_E0)
      }

      V2 <- THETA_V2 * (WT / 70)
      Q <- THETA_Q * (WT / 70)^0.75

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 5.05,
      THETA_CL_2 = 0.00932,
      THETA_V1_E1 = 7.38,
      THETA_V1_E0 = 16.5,
      THETA_V2 = 6.27,
      THETA_Q = 6.28
    ),
    covariates = c("eGFR_CKD_EPI_CYSTATIN", "WT", "ECMO"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1_E1 + ETA_V1_E0 ~
        c(
          log((33.7 / 100)^2 + 1),
          0, log((45.9 / 100)^2 + 1),
          0, 0, log((65.3 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.269)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC7318020/
  # Pipe_8
  Klastrup_2020 = list(
    ppk_model = rxode2::rxode({
      AUC(0) <- 0
      centr(0) <- 0

      CL <- THETA_Cl * (THETA_Cl_2 * CLCR_CG) * exp(ETA_Cl)
      V1 <- THETA_V1

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
      THETA_Cl = 2.25,
      THETA_Cl_2 = 0.119, 
      THETA_V1 = 35.8
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({ # Provided as CV
      ETA_Cl ~
        c(log((57.4 / 100)^2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.226)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC4604352/
  # Pipe_9
  Öbrink_Hansen_2015 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL - THETA_CL_2 * (CREA_micromol_l - 170)) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 3.6,
      THETA_CL_2 = 0.011,
      THETA_V1 = 7.3,
      THETA_V2 = 6.58,
      THETA_Q = 3.9
    ),
    covariates = c("CREA_micromol_l"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((71.2 / 100)^2 + 1),
          0, log((57.8 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.147)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/20018492/
  # Pipe_10
  Roberts_2009 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (WT / 70)) * exp(ETA_CL + KAPPA_CL)
      V1 <- THETA_V1 * exp(ETA_V1 + KAPPA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q <- THETA_Q * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc

      alag(centr) <- THETA_ALAG * exp(ETA_ALAG) # time lag from dose infuser to patient
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 17.1,
      THETA_V1   = 7.2,
      THETA_V2   = 17.8,
      THETA_Q    = 52.0,
      THETA_ALAG = 0.07
    ),
    covariates = c("WT"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_V2 + ETA_Q + ETA_ALAG ~
        c(
          log((29.8 / 100)^2 + 1),
          0, log((26.4 / 100)^2 + 1),
          0, 0, log((73.2 / 100)^2 + 1),
          0, 0, 0, log((50.2 / 100)^2 + 1),
          0, 0, 0, 0, log((43.7 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 3.2, proportional_b = 0.253),
    pi_matrix = lotri::lotri({ # Provided as CV
      KAPPA_CL + KAPPA_V1 ~
        c(
          log((46.2 / 100)^2 + 1),
          0, log((24.4 / 100)^2 + 1)
        )
    })
  ),
  # https://pubmed.ncbi.nlm.nih.gov/35625262/
  # Pipe_11
  Selig_05_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (CLCR_CG / 130)^THETA_CL_2) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 17.56,
      THETA_CL_2 = 0.65,
      THETA_V1   = 33.59,
      THETA_V2   = 10.5,
      THETA_Q    = 6.8
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({ # Provided as omega squared/SD^2/variance/IIV
      ETA_CL + ETA_V1 ~
        c(
          0.17,
          0, 0.34
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.240)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/30963365/
  # Pipe_12
  Sukarnjanaset_2019 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL + THETA_CL_2 * (CLCR_J - 55) + THETA_CL_3 * (MAP - 68)) * exp(ETA_CL)
      V1 <- (THETA_V1 + THETA_V1_2 * (ABW - 56)) * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 5.37,
      THETA_CL_2 = 0.06,
      THETA_CL_3 = 0.05,
      THETA_V1   = 9.35,
      THETA_V1_2 = 0.26,
      THETA_V2   = 7.77,
      THETA_Q    = 21.3
    ),
    covariates = c("CLCR_J", "MAP", "ABW"), # Mean Arterial Pressure
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((28.5 / 100)^2 + 1),
          0, log((55.4 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.223)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/25632974/
  # Pipe_13
  # took the mean values---> should take the bootstrapped ones?
  Udy_2015 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (CLCR_CG / 100)) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
      
      #alag(centr) <- THETA_ALAG * exp(ETA_ALAG)
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 16.2,
      THETA_V1   = 15.9,
      THETA_V2   = 21.3,
      THETA_Q    = 42.9#,
      #THETA_ALAG = 0.13
    ),
    covariates = c("CLCR_CG"),
    # omega = lotri::lotri({ # Provided as CV
    #   ETA_CL + ETA_V1 + ETA_V2 + ETA_ALAG ~
    #     c(
    #       log((55.3 / 100)^2 + 1),
    #       0, log((23.5 / 100)^2 + 1),
    #       0, 0, log((69.2/100)^2 + 1),
    #       0, 0, 0, log((0.4/100)^2 +1)
    #     )
    # }),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_V2  ~
        c(
          log((55.3 / 100)^2 + 1),
          0, log((23.5 / 100)^2 + 1),
          0, 0, log((69.2/100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.3, proportional_b = 0.01)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/26869692/
  # Pipe_14
  Ulldemolins_2016 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      periph(0) <- 0
      
      if (MEMB == 1.5) {
        CL <- (THETA_CL * (WT / 80)^THETA_CL_2) * exp(ETA_CL)
      } else {
        CL <- (THETA_CL * (WT / 80)^THETA_CL_2 * (1 - THETA_CL_3)) * exp(ETA_CL)
      }
      
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 6.11,
      THETA_CL_2 = 1.39,
      THETA_CL_3 = 0.49,
      THETA_V1   = 19.4,
      THETA_V2   = 12.9,
      THETA_Q    = 9.5
    ),
    covariates = c("WT","MEMB"), # Influence of the AN69 membrane. filters 0.9 AN69 and 1.5 AN69ST
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((17.54 / 100)^2 + 1),
          0, log((52.2 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 13.3, proportional_b = 0.06)
    ),
    # https://pmc.ncbi.nlm.nih.gov/articles/PMC9249689/
    # Pipe_15
    Wallenburg_2022 = list(
      ppk_model = rxode2::rxode2({
        AUC(0)    <- 0
        centr(0)  <- 0
        periph(0) <- 0
        
        CL_renal     <- THETA_CL_renal * eGFR_CKD_EPI * exp(KAPPA_CL_renal)   
        CL_nonrenal  <- THETA_CL_nonrenal * exp(ETA_CL_nonrenal)  
        CL <- CL_renal + CL_nonrenal
        V1 <- THETA_V1 * exp(ETA_V1)
        V2 <- THETA_V2 * exp(ETA_V2)
        Q <- THETA_Q

        k10 <- CL_nonrenal / V1     
        k15 <- CL_renal    / V1     
        k12 <- Q / V1
        k21 <- Q / V2
        

        Cc <- centr / V1

        d / dt(centr) <- -(k10 + k12 + k15) * centr + k21 * periph
        d / dt(periph) <- k12 * centr - k21 * periph 
        d / dt(AUC) <- Cc
      }),
      error_model = function(f, sigma) {
        g <- sigma[1] + sigma[2] * f
        return(g)
      },
      theta = c(
        THETA_CL_renal    = 0.0613,
        THETA_CL_nonrenal = 4.11,
        THETA_V1          = 29.7,
        THETA_V2          = 4.26,
        THETA_Q           = 3.27
      ),
      covariates = c("eGFR_CKD_EPI"), # They actually used 24-h urine creatine clearance
      omega = lotri::lotri({ # Provided as CV
        ETA_CL_nonrenal + ETA_V1 + ETA_V2 ~
          c(
            log((95.8 / 100)^2 + 1),
            0, log((44 / 100)^2 + 1),
            0, 0, log((336 / 100)^2 + 1)
          )
      }),
      sigma = c(additive_a = 0, proportional_b = 0.28),
      pi_matrix = lotri::lotri({ # Provided as CV
        KAPPA_CL_renal ~
          c(log((31.5 / 100)^2 + 1))
      })
    ),
  # https://pubmed.ncbi.nlm.nih.gov/39628093/
  # Pipe_16
  Reeder_2025 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL_renal <- THETA_CL_renal * (CLCR_CG_LBW/51)^THETA_CL_2 
      CL_nonrenal <- THETA_CL_nonrenal * (LBW/57)^THETA_CL_3
      CL_final  <- ((CL_renal + CL_nonrenal) * (1 - RRT) + (THETA_RRT * RRT)) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)

      k10 <- CL_final / V1

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL_renal    = 1.97,
      THETA_CL_2        = 1,
      THETA_CL_nonrenal = 2.58,
      THETA_CL_3        = 0.75,
      THETA_RRT         = 3.32,
      THETA_V1          = 23.5
    ),
    covariates = c("CLCR_CG_LBW", "LBW", "RRT"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((47 / 100)^2 + 1),
          0, log((33 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.326, proportional_b = 0.39)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC11434833/
  # Pipe_17---> DA SCRIVERE
  Carriére_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL1 <- (THETA_CL1 * exp(beta_CL1_RRT * RRT + beta_CL1_RI * RI)) * exp(ETA_CL1 + KAPPA_CL1)
      CL2 <- THETA_CL2 * exp(ETA_CL2)       # peritoneal elimination
      V1 <- (THETA_V1 * exp(beta_V1_age * AGE + beta_V1_weight * WT)) * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL1 / V1          # central → elimination (renal)
      k12 <- Q  / V1           # central → peritoneal
      k21 <- Q  / V2           # peritoneal → central
      k20 <- CL2 / V2          # peritoneal → elimination (VAC)

      Cc <- centr / V1

      d/dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d/dt(periph) <- k12 * centr - k21 * periph - k20 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
        g <- sigma[1] + sigma[2] * f
        return(g)
      },
    theta = c(
      THETA_CL1 = 38.739,
      THETA_CL2 = 5.567,
      THETA_V1 = 21.23,
      THETA_V2 = 5,
      THETA_Q = 23.833,
      beta_V1_age    = -0.016,   
      beta_V1_weight =  0.017,   
      beta_CL1_RRT   = -2.016,   
      beta_CL1_RI    = -1.243    
    ),
    covariates = c("AGE", "WT",
                   "RRT",        # binary: 1 if CRRT (group B), else 0
                   "RI"),        # binary: 1 if moderate RI (group C), else 0,
    omega = lotri::lotri({       # Provided as omegas
      ETA_CL1 + ETA_CL2 + ETA_V1  ~
        c(
          0.661^2,
          0, 0.123^2,
          0, 0, 0.225^2
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.207),
    pi_matrix = lotri::lotri({ # Provided as omega
      KAPPA_CL1 ~ 0.391^2
    })
  ),
  # https://pubmed.ncbi.nlm.nih.gov/39699210/
  # Pipe_18
  Kumta_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)     <- 0
      centr(0)   <- 0
      periph1(0) <- 0
      periph2(0) <- 0

      CL <- THETA_CL * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
      V3 <- THETA_V3
      Q  <- THETA_Q
      Q3 <- THETA_Q3 * exp(ETA_Q3)
      
      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2
      k13 <- Q3 / V1
      k31 <- Q3 / V3

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph1 - k13 * centr + k31 * periph2
      d / dt(periph1) <- k12 * centr - k21 * periph1
      d / dt(periph2) <- k13 * centr - k31 * periph2
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 12.7,
      THETA_V1 = 13.4,
      THETA_V2 = 4.99,
      THETA_V3 = 0.16,
      THETA_Q  = 7.25,
      THETA_Q3 = 0.00024
    ),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_Q3 ~
        c(
          log((29.5 / 100)^2 + 1),
          0, log((27 / 100)^2 + 1),
          0, 0, log((86.3 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 2.10, proportional_b = 0.08)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/22503878/
  # Pipe_19
  Delattre_2012 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0

      CL <- (THETA_CL * (CLCR_CG/100)^THETA_CL_2) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
      Q <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

      Cc <- centr / V1

      d / dt(centr) <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC) <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 0.97,
      THETA_CL_2 = 0.62,
      THETA_V1   = 24,
      THETA_V2   = 8.11,
      THETA_Q    = 5.48
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1  ~
        c(
          log((51.5 / 100)^2 + 1),
          0, log((45.4 / 100)^2 + 1)
        )
     }),
    sigma = c(additive_a = 1.90, proportional_b = 0.307)
  )
)
