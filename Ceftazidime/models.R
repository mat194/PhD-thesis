library(posologyr)

# https://doi.org/10.3390/antibiotics10060612
models_ceftazidime <- list(
  #Ceftaz_1
  #https://pubmed.ncbi.nlm.nih.gov/34063815/
  WerumeusBuning_2021 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      if (CVVH == 1) {
        CL <- THETA_CL_CVVH
      } else {
        CL <- THETA_CL *
          (eGFR / 73)^THETA_CL_eGFR *
          THETA_CL_HEMAT^HEMAT *
          THETA_CL_TRAUMA^TRAUMA *
          exp(ETA_CL)
      }
      V <- THETA_V * exp(ETA_V)

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
      THETA_CL         = 3.42,   
      THETA_CL_CVVH    = 2.9,    
      THETA_CL_eGFR    = 0.772,  
      THETA_CL_HEMAT   = 1.57,   
      THETA_CL_TRAUMA  = 1.99,   
      THETA_V          = 46.8    
    ),
    covariates = c("eGFR", "CVVH", "HEMAT", "TRAUMA"),
    omega = lotri::lotri({ 
      ETA_CL + ETA_V ~     
        c(log((36 / 100)^2 + 1),
          0, log((102.8 / 100)^2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.281)
  ),
  #Ceftaz_2
  # https://doi.org/10.1111/j.1365-2125.2007.02857.x
  Conil_2007 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL_INT + THETA_CL_CLCR * CLCR_CG) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 *
        (1 + THETA_V2_SEX  * SEX) *
        (1 + THETA_V2_VM * VM) *
        (1 + THETA_V2_CLCR * CLCR_CG)
      Q  <- THETA_Q

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

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
      THETA_CL_INT  = 1.08,     
      THETA_CL_CLCR = 0.0536,   
      THETA_V1      = 18.8,     
      THETA_Q       = 6.88,     
      THETA_V2      = 2.69,     
      THETA_V2_SEX  = 1.43,     
      THETA_V2_VM   = 2.44,     
      THETA_V2_CLCR = 0.00414   
    ),
    covariates = c("CLCR_CG", "SEX", "VM"),
    omega = lotri::lotri({ 
      ETA_CL + ETA_V1 ~    
        c(log((16 / 100)^2 + 1),
          0, log((13 / 100)^2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.38)
  ),
  #Ceftaz_3
  #https://pmc.ncbi.nlm.nih.gov/articles/PMC2764213/
  Georges_2009 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0

      CL <- (THETA_CL_INT + THETA_CL_MDRD * MDRD) * exp(ETA_CL)

      if (VM == 1) {
        TVV1 <- THETA_V1_VENT
      } else {
        TVV1 <- THETA_V1_NOVENT
      }
      V1 <- TVV1 * exp(ETA_V1)

      if (ADM == 1) {            # polytrauma
        TVV2 <- THETA_V2_TRAUMA
      } else {
        if (ADM == 2) {          # postoperative
          TVV2 <- THETA_V2_POSTOP
        } else {                 # medical
          TVV2 <- THETA_V2_MED
        }
      }
      V2 <- TVV2 * exp(ETA_V2)

      Q <- THETA_Q * exp(ETA_Q)

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

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
      THETA_CL_INT    =  2.24,   
      THETA_CL_MDRD   =  0.024,  
      THETA_V1_NOVENT = 18.90,   
      THETA_V1_VENT   =  9.02,   
      THETA_Q         = 15.20,   
      THETA_V2_TRAUMA = 57.10,   
      THETA_V2_POSTOP = 25.70,   
      THETA_V2_MED    = 13.60    
    ),
    # ADM: reason for admission, 1 = polytrauma, 2 = postoperative, 3 = medical
    covariates = c("MDRD", "VM", "ADM"),
    omega = lotri::lotri({ # $OMEGA of the final model, already variances
      ETA_CL + ETA_V1 + ETA_Q + ETA_V2 ~
        c(0.09,
          0, 0.12,
          0, 0, 0.50,
          0, 0, 0, 0.11)
    }),
    # $SIGMA of the final model is 0.05, i.e. a proportional RUV of sqrt(0.05)
    sigma = c(additive_a = 0, proportional_b = sqrt(0.05))
  ),
  #Ceftaz_4
  # https://doi.org/10.3390/antibiotics13080756
  Launay_2024 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0

      CL <- THETA_CL * (eGFR / 73.90)^THETA_CL_eGFR * exp(ETA_CL)
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
      THETA_CL      =  4.45,   
      THETA_CL_eGFR =  0.9,    
      THETA_V       = 88.0     
    ),
    covariates = c("eGFR"),
    omega = lotri::lotri({ # Table 2 reports IIV as standard deviations
      ETA_CL + ETA_V ~    
        c(0.46^2,
          0, 0.57^2)
    }),
    sigma = c(additive_a = 0.39, proportional_b = 0)
  ),
  #Ceftaz_5
  # https://doi.org/10.1016/j.clinbiochem.2012.03.030
  Delattre_2012 = list(
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0

      CL <- THETA_CL * (WT / 70)^0.75 * (CLCR_70 / 100)^THETA_CL_RF * exp(ETA_CL)
      V1 <- THETA_V1 * (WT / 70) * exp(ETA_V1)
      V2 <- THETA_V2 * (WT / 70) * exp(ETA_V2)
      Q  <- THETA_Q  * (WT / 70)^0.75   

      k10 <- CL / V1
      k12 <- Q / V1
      k21 <- Q / V2

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
      THETA_CL    =  3.91,  
      THETA_CL_RF =  0.29,  
      THETA_V1    = 20.3,   
      THETA_V2    = 14.9,   
      THETA_Q     = 17.9    
    ),
    covariates = c("WT", "CLCR_70"),
    omega = lotri::lotri({ # Table 1 reports IIV as CV%: 
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(log((58.5 / 100)^2 + 1),
          0, log((63.0 / 100)^2 + 1),
          0, 0, log((42.2 / 100)^2 + 1))
    }),
    sigma = c(additive_a = 0, proportional_b = 0.218)
  )
)
