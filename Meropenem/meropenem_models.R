library(posologyr)

meropenem_models <- list(
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9872596/
  # Mero_1
  An_2023 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      if (CRRT == 1) {
        CL <- (THETA_CLNR + THETA_CLCRRT) * exp(ETA_CL)
      } else {
        CL <- (THETA_CLNR + THETA_CLR * (CLCR_TBW / 87)) * exp(ETA_CL)
      }
      V <- THETA_V * (TBW / 85) * exp(ETA_V)
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CLNR   = 1.59,   
      THETA_CLR    = 3.69,   
      THETA_CLCRRT = 2.39,   
      THETA_V      = 35.1    
    ),
    covariates = c("CLCR_TBW", "TBW", "CRRT"),
    omega = lotri::lotri({
      ETA_CL + ETA_V ~
        c(
          log((47.1 / 100)^2 + 1),   
          log((47.1 / 100)^2 + 1),                      
          log((58.9 / 100)^2 + 1)   
        )
    }),
    sigma = c(additive_a = 0.209, proportional_b = 0.316)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9664862/
  # Mero_2
  Boonpeng_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL * ((1 + THETA_CL_CG) * (CLCR_CG - 50))) * exp(ETA_CL)
      V1 <- (THETA_V1 * ((1 + THETA_V1_ALB) * (ALB - 2.5))) * exp(ETA_V1)
      V2 <- THETA_V2 + THETA_V2SHOCK * shock
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
    # This would be the combined2() -->  assumes that the additive and proportional differences are combined on a variance
    # error_model <- function(f,sigma){
    #   g <- sigma[1]^2 + (sigma[2]^2)*(f^2)
    #   return(sqrt(g))
    # },
    theta = c(
      THETA_CL      = 5.16,
      THETA_CL_CG   = 0.016,
      THETA_V1      = 11.1,
      THETA_V1_ALB  = -0.255,
      THETA_V2      = 10.3,
      THETA_V2SHOCK = 11.3,
      THETA_Q       = 9.37
    ),
    covariates = c("CLCR_CG", "ALB", "shock"),
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_V1 ~
        c(
          log((60.5 / 100)^2 + 1),
          0, log((47.7 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.26)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/30304491/
  # Mero_3
  Burger_2018 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      if (CRRT == 1) {
        Sc <- TV_Sc * exp(ETA_Sc)
        TV_CL <- THETA_CL_res + Sc * Qfd
      } else {
        TV_CL <- THETA_CL_no_CRRT * (1 + THETA_CL * (CLCR_CG - 44)/44)
      }
      CL <- TV_CL * exp(ETA_CL)
      V1 <- (THETA_V1 * (TBW/72)^THETA_V1_TBW) * exp(ETA_V1)
      
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
      TV_Sc              = 0.76,
      THETA_CL_res       = 3.2, # meropenem residual clearance in CRRT patients
      THETA_CL           = 0.69,
      THETA_CL_no_CRRT   = 5.9,
      THETA_V1           = 16,
      THETA_V1_TBW       = -0.255
    ),     
    covariates = c("CRRT", "Qfd", "TBW", "CLCR_CG"), # Qfd--> FD flow
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_Sc + ETA_V1 ~
        c(
          log((40 / 100)^2 + 1),
          0, log((26 / 100)^2 + 1),
          0, 0, log((51 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 2.7, proportional_b = 0.041)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/31229669/
  # Mero_4
  Ehmann_2019 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      if (CLCR_CG < 154) {
        TV_CL <- THETA_CL * (1 + THETA_CL_CG1 * (CLCR_CG - 80.8))
      } else {
        TV_CL <- THETA_CL * (1 + THETA_CL_CG2 * (154 - 80.8))
      }
      
      CL <- TV_CL * exp(ETA_CL + KAPPA_CL)
      V1 <- (THETA_v1 * (WT/70)^THETA_WT) * exp(ETA_V1)
      V2 <- (THETA_V2 * (ALB - 2.79)) * exp(ETA_V2)
      Q  <- THETA_Q
      
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
      THETA_CL           = 9.25,
      THETA_CL_CG1       = 0.00977,
      THETA_CL_CG2       = 154,
      THETA_v1           = 7.89,
      THETA_V2           = 16.1,
      THETA_WT           = 0.945,
      THETA_Q            = 28.4
    ),     
    covariates = c( "WT", "CLCR_CG", "ALB"), # Qfd--> FD flow
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((27.1 / 100)^2 + 1),
          0, log((31.5 / 100)^2 + 1),
          0, 0, log((16.9 / 100)^2 + 1)
        )
    }),
    pi_matrix = lotri::lotri({ # Provided as CV
      KAPPA_CL ~ log((12.5 / 100)^2 + 1)
    }),
    sigma = c(additive_a = 2.7, proportional_b = 0.041)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/33515688/
  # Mero_5
  Eisert_2021= list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL * (1 + THETA_CL_CG * (CLCR_CG_IBW - 48.8))) * exp(ETA_CL + KAPPA_CL)
      V1 <- THETA_V1  * exp(ETA_V1 + KAPPA_V1)
      V2 <- THETA_V2  
      Q  <- THETA_Q   
      
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
      THETA_CL           = 7.42,
      THETA_CL_CG        = 0.0154,
      THETA_V1           = 45.8,
      THETA_V2           = 9.18,
      THETA_Q            = 2.21
    ),     
    covariates = c("CLCR_CG_IBW"), 
    omega = lotri::lotri({ # provided as CV
      ETA_CL + ETA_V1  ~
        c(
          log((33.1 / 100)^2 + 1),
          0, log((92.9 / 100)^2 + 1)
        )
    }),
    pi_matrix = lotri::lotri({ # Provided as CV
      KAPPA_CL + KAPPA_V1 ~ 
        c(
          log((29.5 / 100)^2 + 1),
          0, log((13.4 / 100)^2 + 1)
          )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.292)
  ),
  # To be checked
  # https://pubmed.ncbi.nlm.nih.gov/25216821/
  # Mero_6
  # Frippiat_2014 = list(
  #   ppk_model = rxode2::rxode2({
  #     AUC(0)   <- 0
  #     centr(0) <- 0
  #     periph(0) <- 0
  #     ELF(0)   <- 0
  #     
  #     CL  <- THETA_CL * (eGFR_CKD_EPI / 80)^THETA_CL_GFR * exp(ETA_CL)
  #     V1  <- THETA_V1 * (WT  / 70)^THETA_V1_WT  * exp(ETA_V1)
  #     V2  <- THETA_V2                           * exp(ETA_V2) # IIV retained per paper
  #     VE  <- THETA_VE                             
  #     QP  <- THETA_QP                             
  #     QE  <- THETA_QE                             
  #     
  #     k10 <- CL  / V1
  #     k12 <- QP  / V1
  #     k21 <- QP  / V2
  #     k13 <- QE  / V1
  #     k31 <- QE  / VE
  #     
  #     Cc     <- centr / V1   # plasma concentration
  #     Cc_ELF <- ELF   / VE   # ELF concentration
  #     
  #     d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph - k13 * centr + k31 * ELF
  #     d / dt(periph) <-  k12 * centr - k21 * periph
  #     d / dt(ELF)    <-  k13 * centr - k31 * ELF
  #     d / dt(AUC)    <- Cc
  #   }),
  #   error_model = function(f, sigma) {
  #     g <- sigma[1] + sigma[2] * f
  #     return(g)
  #   },
  #   theta = c(
  #     THETA_CL      =  10.2,   
  #     THETA_CL_GFR  =  0.73,  
  #     THETA_V1      =  5.2,   
  #     THETA_V1_WT   =  0.60,  
  #     THETA_V2      =  12.1,   
  #     THETA_VE      =  11.3,   
  #     THETA_QP      =  6.9,   
  #     THETA_QE      =  66.5    
  #   ),
  #   covariates = c("eGFR_CKD_EPI", "WT"), 
  #   omega = lotri::lotri({
  #     ETA_CL + ETA_V1 + ETA_V2 ~
  #       c(
  #         NA_real_,            
  #         0, NA_real_,         
  #         0, 0, NA_real_       
  #       )
  #   }),
  #   sigma = c(additive_a = 0.8, proportional_b = 0.2) 
  # ),
  # https://pubmed.ncbi.nlm.nih.gov/35035223/
  # Mero_7
  Gijsen_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * (CLCR_CG / 111.7)^THETA_CL_CG * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL     = 13.7,   
      THETA_CL_CG  = 0.725, # to be checked
      THETA_V1     = 25.5,   
      THETA_V2     = 12.4,   
      THETA_Q      = 8.1    
    ),
    covariates = c("CLCR_CG"), 
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((57.8 / 100)^2 + 1),
          0.185 ,log((57   / 100)^2 + 1),
          0, 0,log((74.4 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.147, proportional_b = 0.311) 
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC7176801/
  # Mero_8
  Grensemann_2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * exp(ETA_CL)  
      V1 <- THETA_V1                 
      
      if (ACLF == 1) {
        V2 <- THETA_V2_ACLF
      } else {
        V2 <- THETA_V2_NLF
      }
      
      Q <- THETA_Q * exp(ETA_Q)
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL      =  5.06,  
      THETA_V1      =  8.31,  
      THETA_V2_ACLF =  35.5,   
      THETA_V2_NLF  =  18.5,   
      THETA_Q       =  7.23   
    ),
    covariates = c("ACLF"), # acute-on-chronic liver failure
    omega = lotri::lotri({
      ETA_CL + ETA_Q ~
        c(
          log((29.8 / 100)^2 + 1),
          0, log((27.6 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.220) 
  ),
  # https://pubmed.ncbi.nlm.nih.gov/29530848/
  # Mero_9
  Hanberg_2018 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL_frac * CLCR_CG * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)      
      Q  <- THETA_Q  * exp(ETA_Q)
      V2 <- THETA_V2                     
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc      <- centr  / V1                      # plasma concentration
      Cc_SCT  <- (periph / V2) * THETA_fu_tissue  # free concentration in SCT
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL_frac  = 0.0460, 
      THETA_V1       = 8.31,   
      THETA_Q        = 8.52,   
      THETA_V2       = 6.99,   
      THETA_fu_tissue = 0.790  
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_Q ~
        c(
          log((18.9 / 100)^2 + 1),
          0, log((25.1 / 100)^2 + 1),
          0, 0, log((64.6 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.199)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/18307371/
  # Mero_10
  Isla_2008 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      if (SEPSIS == 1) {
        CL <- (THETA_CL + THETA_CLCR_SEP  * CLCR_CG) * exp(ETA_CL)
        V1 <- THETA_V1_SEP  * exp(ETA_V1)
      } else {
        CL <- (THETA_CL + THETA_CLCR_POLY * CLCR_CG) * exp(ETA_CL)
        V1 <- THETA_V1_POLY * exp(ETA_V1)
      }
      
      V2 <- THETA_V2    
      Q  <- THETA_Q     
      
      #Sc <- THETA_Sc * exp(ETA_Sc)
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc   <- centr / V1        
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL        =  6.63,  
      THETA_CLCR_SEP  =  0.064, 
      THETA_CLCR_POLY =  0.72,  
      THETA_V1_SEP    =  15.7,   
      THETA_V1_POLY   =  69.5,   
      THETA_V2        =  19.8,   
      THETA_Q         =  15.0#,   
      #THETA_Sc        =  0.72   
    ),
    covariates = c("CLCR_CG", "SEPSIS"),
    omega = lotri::lotri({
      # ETA_Sc
      ETA_CL + ETA_V1 ~
        c(
          log((38 / 100)^2 + 1),
          0, log((45 / 100)^2 + 1)
        )
          #,
          #0, 0, log((18 / 100)^2 + 1)
    }),
    sigma = c(additive_a = 0.22, proportional_b = 0.21)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/25753628/
  # Mero_11
  Jaruratanasirikul_2015 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL <- (THETA_CL_int + THETA_CL_MDRD * eGFR_MDRD) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      
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
      THETA_CL_int  = 3.01, 
      THETA_CL_MDRD = 0.07, 
      THETA_V1      = 23.7
    ),
    covariates = c("eGFR_MDRD"), # eGFR by MDRD 
    omega = lotri::lotri({
      ETA_CL + ETA_V1 ~
        c(
          log((48 / 100)^2 + 1),
          0, log((35 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0) # not reported
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC9693387/
  # Mero_12
  Kang_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      if (CRRT == 1) {
        CL <- THETA_CL * THETA_CRRT_eff * exp(ETA_CL)
      } else {
        CL <- THETA_CL * exp(ETA_CL)
      }
      
      V1 <- THETA_V1             
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q              
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL       =  3.79,  
      THETA_CRRT_eff =  0.44, 
      THETA_V1       =  2.4,   
      THETA_V2       =  8.56,  
      THETA_Q        = 21.3    
    ),
    covariates = c("CRRT"),
    omega = lotri::lotri({
      ETA_CL + ETA_V2 ~
        c(
          log((47.1 / 100)^2 + 1),
          0, log((44.0 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.473)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC8625191/
  # Mero_13
  Lee_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * (1 + THETA_CL_eGFR * (eGFR_CKD_EPI - 91.57)) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q 
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    # Power and proportional residual error model: SD = sigma[1] × f^sigma[2]
    error_model = function(f, sigma) {
      g <- sigma[1] * f^sigma[2]
      return(g)
    },
    theta = c(
      THETA_CL      =  6.37,    
      THETA_CL_eGFR =  0.00925, 
      THETA_V1      =  9.07,    
      THETA_Q       =  10.7,     
      THETA_V2      =  7.91     
    ),
    covariates = c("eGFR_CKD_EPI"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((31.4 / 100)^2 + 1),
          0, log((43.6 / 100)^2 + 1),
          0, 0, log((36.6 / 100)^2 + 1)
        )
    }),
    sigma = c(proportional_a = 0.246, power_b = 0.865) # understand if this is correct
    ),
  # https://pubmed.ncbi.nlm.nih.gov/34790131/
  # Mero_14
  Ha_Lee_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * (1 + THETA_CL_2 * (CLCR_CG - 49.7)) * exp(ETA_CL)
      V1 <- THETA_V1 * (1 + THETA_V1_2 * (LPM   -  3.7))   * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 7.35,    
      THETA_CL_2 = 0.0104,  
      THETA_V1   = 17.3,    
      THETA_V1_2 = 0.337,   
      THETA_Q    = 14.5,    
      THETA_V2   = 12.8     
    ),
    # LPM: ECMO flow rate (L/min)
    covariates = c("CLCR_CG", "LPM"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((39.6 / 100)^2 + 1),            
          0, log((48.5 / 100)^2 + 1),         
          0, 0, log((38.8 / 100)^2 + 1)       
        )
    }),
    sigma = c(additive_a = 0.370, proportional_b = 0.0473)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/29425283/
  # Mero_15
  # it was built using NONMEM's $PRED subroutine on steady-state concentrations only during continuous infusion. The paper explicitly states that V could not be estimated.
  Minichmayr_2018 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL <- THETA_CL * (1 + THETA_COV * (CLCR_CG - 80)) * exp(ETA_CL)
      V  <- THETA_V   # nominal fixed value; NOT estimated from data
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL  = 7.71,    
      THETA_COV = 0.00924, 
      THETA_V   = 22.0     # NOT estimated; nominal Vd assumed from meropenem literature
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({
      ETA_CL ~ log((25.7 / 100)^2 + 1)
    }),
    sigma = c(additive_a = 0, proportional_b = 0.232)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/32049890/
  # Mero_16
  Yoko_2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * (eGFR_CKD_EPI/26.1)^THETA_CL_2 * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 4.22,    
      THETA_CL_2 = 0.25,  
      THETA_V1   = 14.82,    
      THETA_Q    = 7.84,    
      THETA_V2   = 11.75    
    ),
    covariates = c("eGFR_CKD_EPI"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((26.75 / 100)^2 + 1),            
          0, log((27.33 / 100)^2 + 1),         
          0, 0, log((57.05 / 100)^2 + 1)       
        )
    }),
    sigma = c(additive_a = 1.08, proportional_b = 0.1839)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/34403127/
  # Mero_17
  OJeanson_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL1 <- THETA_CL + THETA_GFR * eGFR_MDRD
      CL2 <- CL1 * (1 - DIA_C) + THETA_DIA_C * DIA_C
      CL3 <- CL2 * (1 - RFS) + THETA_RFS * RFS
      CL4 <- CL3 + THETA_RD * (DIUR / 845 - 1)
      CL5 <- CL4 * (1 - DIA_SC) + THETA_DIA_SC * DIA_SC
      
      CL <- CL5 * exp(ETA_CL)
      V  <- THETA_V * exp(ETA_V)
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL     = 1.36,   
      THETA_GFR    = 0.058,  
      THETA_DIA_C  = 6.38,   
      THETA_RFS    = 13.9,   
      THETA_RD     = 0.60,   
      THETA_DIA_SC = 11.0,   
      THETA_V      = 43.9    
    ),
    # eGFR_MDRD---> 186.3 × (creatinine [micromol/L] × 0.0113)^−1.154 × age^−0.203 × 0.742 (if female) × 1.21 (if African)
    # RFS ---> Renal Function Status; 1 if eGFR_MDRD ≥ 120 mL/min/1.73m², else 0
    # DIA_C ---> 1 if Continuous dialysis
    # DIA_S ---> 1 if Semi-continuous dialysis
    # DIUR ---> Residual diuresis volume in mL/24h
    covariates = c("eGFR_MDRD", "DIA_C", "DIA_SC", "RFS", "DIUR"),
    omega = lotri::lotri({
      ETA_CL + ETA_V ~
        c(
          log((30.8 / 100)^2 + 1),           
          0, log((87.6 / 100)^2 + 1)        
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.321)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/32301057/
  # Mero_18
  Onichimowski_2020 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      V1 <- THETA_V1 * (ALB / 24.6)^THETA_BETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      CL <- THETA_CLCRRT * exp(ETA_CL)
      Q  <- THETA_Q   
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_V1      = 27.9,   
      THETA_BETA_V1 = -2.87,  
      THETA_CLCRRT  = 15.1,   
      THETA_Q       = 21.1,   
      THETA_V2      = 33.7    
    ),
    covariates = c("ALB"),
    omega = lotri::lotri({
      ETA_V1 + ETA_CL + ETA_V2 ~
        c(
          log((53.1 / 100)^2 + 1),      
          0, log((43.7 / 100)^2 + 1),   
          0, 0, log((85.6 / 100)^2 + 1) 
        )
    }),
    sigma = c(additive_a = 0.881, proportional_b = 0.241)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/31335959/
  # Mero_19
  Padulles_Zamora_2019 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL  <- THETA_CL  * exp(ETA_CL)   
      V1  <- THETA_V1  * exp(ETA_V1)   
      CLD <- THETA_CLD                
      V2  <- THETA_V2                 
      
      k10 <- CL  / V1
      k12 <- CLD / V1
      k21 <- CLD / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL  =   7.78,   
      THETA_V1  =   24.9,    
      THETA_CLD =   6.49,   
      THETA_V2  =  283.0     
    ),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 ~
        c(
          0.258,          
          0, 0.191        
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.396)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/19398460/
  # Mero_20
  Roberts_2009 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- THETA_CL * (CLCR_CG / 100) * exp(ETA_CL + KAPPA_CL)
      V1 <- THETA_V1 * exp(ETA_V1 + KAPPA_V1)
      V2 <- THETA_V2 * exp(ETA_V2 + KAPPA_V2)
      Q  <- THETA_Q  * exp(ETA_Q  + KAPPA_Q)
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <-  k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 13.6,   
      THETA_V1 =  7.9,   
      THETA_V2 = 14.8,   
      THETA_Q  = 56.3    
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({
      ETA_CL + ETA_V1 + ETA_Q + ETA_V2 ~
        c(
          log((15.3 / 100)^2 + 1),                    
          0, log((44.7 / 100)^2 + 1),                 
          0, 0, log((22.3 / 100)^2 + 1),              
          0, 0, 0, log(( 8.4 / 100)^2 + 1)            
        )
    }),
    pi_matrix = lotri::lotri({
      KAPPA_CL + KAPPA_V1 + KAPPA_Q + KAPPA_V2 ~
        c(
          log((11.8 / 100)^2 + 1),                    
          0, log((27.7 / 100)^2 + 1),                 
          0, 0, log((11.8 / 100)^2 + 1),              
          0, 0, 0, log((28.8 / 100)^2 + 1)            
        )
    }),
    sigma = c(additive_a = 0.43, proportional_b = 0.152)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/34773921/
  # Mero_21
  Selig_2022 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      # CL_CVVH = Qf * Sc * CF, must be pre-calculated and passed as covariate
      if (CVVH == 0) {
        TV_CL <- THETA_CL * (CLCR_CG / 125)^THETA_CL_CrCl
      } else if (BURN == 1) {
        TV_CL <- THETA_CL_BURN_CVVH * (LBM / 61.7)^0.75 + CL_CVVH
      } else {
        TV_CL <- THETA_CL_NOBURN_CVVH * (LBM / 61.7)^0.75 + CL_CVVH
      }
      
      CL <- TV_CL * exp(ETA_CL)
      V1 <- THETA_V1 * (LBM / 61.7) * exp(ETA_V1)
      V2 <- (THETA_V2 + THETA_Vp_TBSA * TBSA) * exp(ETA_V2)
      Q  <- THETA_Q
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
      
      Cc <- centr / V1
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL             = 15.3,    
      THETA_CL_CrCl        = 0.74,    
      THETA_CL_BURN_CVVH   = 12.85,   
      THETA_CL_NOBURN_CVVH = 6.43,    
      THETA_V1             = 33.21,   
      THETA_V2             = 13.45,   
      THETA_Vp_TBSA        = 0.71,    
      THETA_Q              = 14.1     
    ),
    covariates = c(
      "CLCR_CG",  # Cockcroft-Gault CrCl, mL/min (non-CVVH patients), TBW
      "LBM",  # lean body mass --> Janmahasatian's formula
      "CVVH",     
      "BURN",     
      "TBSA",     
      "CL_CVVH"   # pre-calculated extracorporeal CL = Qf * Sc * CF-->correction factor for pre-filter fluid administration, L/h
    ),
    omega = lotri::lotri({ # provided as omega squared (variance)
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          0.31,
          0, 0.35,
          0, 0, 0.61
        )
    }),
    sigma = c(additive_a = 0.0056, proportional_b = 0.22)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/25636084/
  # Mero_22
  Shekar_2014 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      if (RRT == 1) {
        TV_CL <- THETA_CL_RRT
      } else {
        TV_CL <- THETA_CLCRCL * (CLCR_CG * 0.06)
      }
      
      CL <- TV_CL * exp(ETA_CL)
      Vc <- THETA_Vc * exp(ETA_Vc)
      Vp <- THETA_Vp * exp(ETA_Vp)
      Q  <- THETA_Q
      
      k10 <- CL / Vc
      k12 <- Q  / Vc
      k21 <- Q  / Vp
      
      Cc <- centr / Vc
      
      d / dt(centr)  <- -k10 * centr - k12 * centr + k21 * periph
      d / dt(periph) <- k12 * centr - k21 * periph
      d / dt(AUC)    <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL_RRT  = 5.1,   
      THETA_CLCRCL  = 1.89,  
      THETA_Vc      = 18.7,  
      THETA_Vp      = 13.2,  
      THETA_Q       = 21.0   
    ),
    covariates = c(
      "CLCR_CG",  
      "RRT"       
    ),
    omega = lotri::lotri({ # provided as %CV
      ETA_CL + ETA_Vc + ETA_Vp ~
        c(
          log((51.6 / 100)^2 + 1),
          0, log((45.8 / 100)^2 + 1),
          0, 0, log((28.7 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 2.3, proportional_b = 0.137)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/26124172/
  # Mero_23
  Ulldemolins_2015 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL <- (THETA_CL + THETA_DIUR * (DIUR / 100)) * exp(ETA_CL)
      V  <- THETA_V * (WT / 73)^THETA_WT * exp(ETA_V)
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL   = 3.68,   
      THETA_DIUR = 0.22,   
      THETA_V    = 33.00,  
      THETA_WT   = 2.07    
    ),
    covariates = c(
      "DIUR",  # residual diuresis, mL/24h (24-h urine output)
      "WT"     
    ),
    omega = lotri::lotri({ # provided as %CV
      ETA_CL + ETA_V ~
        c(
          log((37 / 100)^2 + 1),
          0, log((45 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.0002, proportional_b = -0.258)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/33818823/
  # Mero_24
  # t_Dial must be provided as a TIME-VARYING covariate in the posologyr
  # patient data (one row per sampling/covariate-change time point):
  #   - RRT patients:     t_Dial = hours elapsed since dialysis session began
  #   - Non-RRT patients: t_Dial = 0 at all time points
  Westermann_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      # Residual diuresis converted to centilitres (1 cL = 10 mL),
      DIUR_cL <- DIUR / 10
      
      # Full CL: creatinine (Eq.1) x diuresis (Eq.2, linear) x time-on-RRT (Eq.3, sigmoidal)
      TV_CL <- THETA_CL *
        (CREA_micromol_l / 191)^THETA_Crea *
        (1 + (DIUR_cL - 4.30) * THETA_diuresis) *
        (1 + t_Dial / (THETA_D50 + t_Dial))
      
      CL <- TV_CL * exp(ETA_CL)
      V  <- THETA_V
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL       = 3.52,   
      THETA_Crea     = -0.45,  
      THETA_diuresis = 0.02,   
      THETA_D50      = 11.6,   
      THETA_V        = 11.3    
    ),
    covariates = c(
      "CREA_micromol_l",    # serum creatinine, μmol/L (pre-dialysis value)
      "DIUR",    # residual diuresis prior to dialysis, mL/24h
      "t_Dial"   # TIME-VARYING: hours elapsed since current RRT session started;
                 # = 0 for non-RRT patients and at the moment dialysis begins
    ),
    omega = lotri::lotri({ # IIV on CL only (%CV); V not identifiable under CI steady-state
      ETA_CL ~ log((32.6 / 100)^2 + 1)
    }),
    sigma = c(additive_a = 2.47, proportional_b = 0.118)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/39082803/
  # Mero_25
  Wicha_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      csf(0)   <- 0
      
      V  <- THETA_V  * (WT / 70)       * exp(ETA_V)
      CL <- THETA_CL * (WT / 70)^0.75  * exp(ETA_CL)
      KP <- (THETA_KP + (log10(IL6_CSF) - 2.5) * THETA_COVKP_IL6) * exp(ETA_KP)
      QCSF <- THETA_FQCSF * exp(ETA_FQCSF) * 0.57 * 0.12 * 15 * WT^0.74
      VCSF <- 0.15
      
      Cc     <- centr / V
      Cc_csf <- csf   / VCSF
      
      k_in  <- QCSF / V           # plasma → CSF direction 
      k_out <- QCSF / (KP * VCSF) # CSF → plasma direction 
      
      d / dt(centr) <- -(CL / V) * centr - k_in * centr + k_out * csf
      d / dt(csf)   <-              k_in * centr - k_out * csf
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {  # plasma error model 
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_V         = 36.2,   
      THETA_CL        = 11.1,   
      THETA_FQCSF     = 2.49e-5,
      THETA_KP        = 0.0602, 
      THETA_COVKP_IL6 = 0.0154  
    ),
    covariates = c(
      "WT",      
      "IL6_CSF"  # IL-6 concentration in CSF, ng/L (= pg/mL)
    ),
    omega = lotri::lotri({ # provided as %CV
      ETA_CL + ETA_V + ETA_KP + ETA_FQCSF ~
        c(
          log((41.2 / 100)^2 + 1),
          0, log((49.9 / 100)^2 + 1),
          0, 0, log((84.0 / 100)^2 + 1),
          0, 0, 0, log((76.7 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.705, proportional_b = 0.226)
  ),
  # https://pmc.ncbi.nlm.nih.gov/articles/PMC11278387/
  # Mero_26
  Rancic_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      TV_CL <- (THETA_CL * CREA_micromol_l^THETA_CRE * WBC^THETA_WBC +
                  THETA_HTA * HTA +
                  THETA_VAN * VAN +
                  THETA_COL * COL)
      
      CL <- TV_CL * exp(ETA_CL)
      V  <- THETA_V             
      
      k10 <- CL / V
      Cc  <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL  = 5.29,       
      THETA_CRE = 0.000001,   
      THETA_WBC = -0.165,     
      THETA_HTA = 0.000001,   
      THETA_VAN = 0.825,      
      THETA_COL = 1.28,       
      THETA_V   = 2.05        
    ),
    covariates = c(
      "CREA_micromol_l",  # serum creatinine, µmol/L
      "WBC",  # white blood cell count, ×10⁹/L
      "HTA",  # hypertension: binary (1 = yes, 0 = no)
      "VAN",  # concomitant vancomycin: binary (1 = yes, 0 = no)
      "COL"   # concomitant colistimethate: binary (1 = yes, 0 = no)
    ),
    omega = lotri::lotri({ # IIV on CL only; omega² = 0.0215 
      ETA_CL ~ 0.0215
    }),
    sigma = c(additive_a = 0, proportional_b = 0.44) # to check the proportional
  ),
  # https://pubmed.ncbi.nlm.nih.gov/39455501/
  # Mero_27
  Troisi_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0) <- 0
      centr(0) <- 0
      
      # V fixed at 20 L (unfeasible to estimate under CI steady-state)
      CL  <- THETA_CL * exp(THETA_eGFR * eGFR_CKD_EPI) * exp(ETA_CL)
      V   <- 20                     
      
      k10 <- CL / V
      Cc  <- centr / V            
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
      
      # --- PD: C-RP indirect turnover model ---
      # IC50 modified by MIC (mg/L) and DRUG (binary: 1 = empirical anti-MRSA)
      # IC50_i  <- THETA_IC50 * exp(THETA_MIC * MIC + THETA_DRUG * DRUG) * exp(ETA_IC50)
      # kout_i  <- THETA_kout * exp(ETA_kout)
      # CRP0_i  <- THETA_CRP0 * exp(ETA_CRP0)
      # kin_i   <- CRP0_i * kout_i     # baseline production rate (ensures CRP(0) = CRP0)
      
      # Full inhibition of production model:
      # dCRP/dt = kin × (1 - Cc/(Cc + IC50)) - kout × CRP
      # crp(0) <- CRP0_i
      # d / dt(crp) <- kin_i * (1 - Cc / (Cc + IC50_i)) - kout_i * crp
    }),
    error_model = function(f, sigma) {   # meropenem PK error model
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL    = 2.72,    # CL at eGFR = 0 (intercept of exponential), L/h
      THETA_eGFR  = 0.011#,   # exponent of eGFR effect on CL (per mL/min/1.73m²)
      #THETA_CRP0  = 22.58,   # baseline C-RP, mg/dL
      # THETA_IC50  = 0.79,    # IC50 at MIC = 0, no anti-MRSA, mg/L
      # THETA_MIC   = 0.44,    # exponent of MIC on IC50 (per mg/L)
      # THETA_DRUG  = 2.46,    # exponent of empirical anti-MRSA on IC50
      # THETA_kout  = 0.012    # C-RP degradation rate constant, h⁻¹
    ),
    covariates = c(
      "eGFR_CKD_EPI"#,  
      #"MIC",   # pathogen MIC for meropenem, mg/L 
      #"DRUG"   # binary: 1 = empirical anti-MRSA co-therapy, 0 = monotherapy 
    ),
    # omega = lotri::lotri({ # Table 3 reports ω as SD → omega matrix uses ω²
    #   ETA_CL + ETA_CRP0 + ETA_IC50 + ETA_kout ~
    #     c(
    #       0.59^2,                
    #       0, 0.51^2,             
    #       0, 0, 1.59^2,          
    #       0, 0, 0, 1.08^2        
    #     )
    # }),
    omega = lotri::lotri({ # Table 3 reports ω as SD → omega matrix uses ω²
      ETA_CL ~
        c(
          0.59^2
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.37)
  )
)
