library(posologyr)

cefepime_models <- list(
  
  #Cefep_1
  #https://pubmed.ncbi.nlm.nih.gov/16943209/
  Roos_2006 = list(
    # 13 ICU patients (11 males), cefepime 2 g q12h as a 3 min iv infusion,
    # sampled on two occasions (day 1 and days 3-6). Serum creatinine
    # < 0.1 mmol/L at inclusion; measured CLCR 2.3-11.7 L/h (median 7.1),
    # TBW 56-128 kg (median 75), age 34-75 years (median 60).
    # Final model: 3-compartment with zero-order input; BSV on CL, V1 and V2
    # (full CL-V1 variance-covariance block, correlation not published ->
    # off-diagonals set to 0 here), BOV on CL and V1. V3 and Q3 were fixed.
    # Covariate model (Equation 4, Table 2): TVCL = THETA_CL * (CLCR/CLCRstd)
    # with CLCRstd = 7.0 L/h = 116.7 mL/min. CLCR is the MEASURED clearance
    # from an 8-h urine collection, not a Cockcroft-Gault estimate.
    ppk_model = rxode2::rxode2({
      AUC(0)     <- 0
      centr(0)   <- 0
      periph1(0) <- 0
      periph2(0) <- 0
      
      CL <- (THETA_CL * (CLCR_urine / 116.7)) * exp(ETA_CL + KAPPA_CL) #reference is 7 L/h--->116.7 ml/min
      V1 <- THETA_V1 * exp(ETA_V1 + KAPPA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      V3 <- THETA_V3
      Q  <- THETA_Q
      Q3 <- THETA_Q3
      
      k10 <- CL / V1
      k12 <- Q  / V1
      k21 <- Q  / V2
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
      THETA_CL = 6.51,   
      THETA_V1 = 5.74,   
      THETA_V2 = 9.57,   
      THETA_V3 = 7.3,    
      THETA_Q  = 33.8,   
      THETA_Q3 = 3.5     
    ),
    covariates = c("CLCR_urine"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((38 / 100)^2 + 1),                         
          0, log((23 / 100)^2 + 1),
          0, 0, log((34 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0.452, proportional_b = 0.08),
    pi_matrix = lotri::lotri({ # Provided as CV
      KAPPA_CL + KAPPA_V1 ~
        c(
          log((15 / 100)^2 + 1),                      
          0, log((30 / 100)^2 + 1)
        )
    })
  ),
  #Cefep_2
  # https://pubmed.ncbi.nlm.nih.gov/27246236/
  Jonckheere_2016 = list(
    # 20 critically ill patients (16 ICU, 4 coronary care), 208 plasma and 51
    # urine samples over a median follow-up of 3 days. Cefepime given per
    # standard of care (Sanford guide, CKD-EPI based) as a 30 min iv infusion,
    # 2 g every 8 h with normal CLCR. Age median 75 years, TBW 76 kg,
    # Cockcroft-Gault 49 mL/min (31-123), 4/20 patients on intermittent
    # haemodialysis (IHD).
    # Final model (Table 4): 2-compartment with zero-order input and parallel
    # first-order elimination split into renal, non-renal and dialysis
    # clearance. In IHD patients renal clearance was assumed to be absent:
    # off dialysis only CL_other applies, on dialysis CL_other + CL_dialysis.
    # Covariate (Equation 2): CL_renal = THETA_CLR * (1 + THETA_CLR_CG *
    # (CLCR - 3.40)), with Cockcroft-Gault expressed in L/h -- the mL/min
    # covariate used elsewhere in this project is converted with * 0.06.
    # The relation is linear and would turn negative below ~11 mL/min,
    # so it is floored at 0.
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CLCR_Lh  <- CLCR_CG * 0.06 # converted to L/h inside the model, it gets ml/min
      TV_CLR   <- THETA_CLR * (1 + THETA_CLR_CG * (CLCR_Lh - 3.40))
      CL_renal <- max(TV_CLR, 0) * exp(ETA_CLR)
      
      if (DIAL == 1) {
        CL <- THETA_CL_other + DIAL * THETA_CL_dial
      } else {
        CL <- CL_renal + THETA_CL_other
      }
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2
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
      THETA_CLR      = 2.88,   
      THETA_CLR_CG   = 0.368,  
      THETA_CL_other = 0.87,   
      THETA_CL_dial  = 5.74,   
      THETA_V1       = 18.3,   
      THETA_V2       = 11.1,   
      THETA_Q        = 6.63    
    ),
    covariates = c("CLCR_CG", "DIAL"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CLR + ETA_V1 ~
        c(
          log((49.5 / 100)^2 + 1),
          0, log((39.7 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.332)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/37882514/
  #Cefep_3
  Barreto_2023 = list(
    # 120 critically ill adults (Mayo Clinic, 2018-2022), 379 serum samples from
    # a richly sampled (n = 20) and an opportunistically sampled (n = 100)
    # cohort. Cefepime 0.5-2 g every 8-24 h as a 30 min infusion; no extended
    # infusions. Age 61 +/- 17 years, weight 87 +/- 26 kg, 82% sepsis.
    # Patients with AKI stage II or higher, on kidney replacement therapy or on
    # ECMO were excluded, so this model does not cover RRT/ECMO.
    # Final model (Table 3): 2-compartment with first-order elimination.
    # Q and V2 were fixed (low observed interindividual variation).
    # Covariates, both time-dependent: V1 = THETA_V1 * (WT / 70) and
    # CL = THETA_CL * (eGFR / 120), i.e. the eGFR expressed in L/h divided by
    # the 7.2 L/h (= 120 mL/min) reference of the paper. Allometric scaling
    # with a fixed 0.75 exponent was tested but did not improve the fit.
    # eGFR is the 2012 CKD-EPI creatinine-cystatin C equation, de-indexed to
    # mL/min (NOT mL/min/1.73 m2) -- it had the lowest BICc and -2LL of the
    # 11 equations compared (Table 2).
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL * (eGFR_CKD_EPI / 120)) * exp(ETA_CL) # should be eGFR_CKD_EPI_CYSTATIN
      V1 <- (THETA_V1 * (WT / 70)) * exp(ETA_V1)
      V2 <- THETA_V2
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
      THETA_CL = 7.84,    # CL1 (L/h) at eGFR = 120 mL/min
      THETA_V1 = 15.56,   # central volume (L) at 70 kg
      THETA_V2 = 10.63,   # peripheral volume (L), fixed
      THETA_Q  = 7.09     # inter-compartmental clearance (L/h), fixed
    ),
    covariates = c("eGFR_CKD_EPI", "WT"),
    omega = lotri::lotri({ # Provided as SD of the random effects (Monolix omega)
      ETA_CL + ETA_V1 ~
        c(
          0.24^2,
          0, 0.41^2
        )
    }),
    sigma = c(additive_a = 0.97, proportional_b = 0.27)
  ),
  # https://pubmed.ncbi.nlm.nih.gov/36450287/
  #Cefep_4
  Seo_2023 = list(
    # 21 critically ill adults (12 men, 9 women) with hospital-acquired or
    # ventilator-associated pneumonia, prospectively enrolled at a Korean
    # tertiary hospital; 144 plasma samples, rich sampling over the first dose
    # period. Cefepime 2 g as a 30 min infusion every 8, 12 or 24 h.
    # Age 67 years (62-76), weight 59.7 +/- 10.5 kg, albumin 2.87 g/dL,
    # Cockcroft-Gault CLCR 77.2 mL/min (45.3-125), APACHE II 17.
    # Final model (Table 2): 2-compartment, IIV on CL, V1, Q and V2.
    # Covariate: CL = THETA_CL * (CLCR_CG / 77.21)^THETA_CL_CG, with
    # Cockcroft-Gault CLCR in mL/min. Including it dropped the IIV on CL from
    # 59.0% to 33.7%. Note the small cohort: weight was NOT retained as a
    # covariate, so the volumes are not weight-scaled despite the low mean
    # body weight (~60 kg) of this population.
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL * (CLCR_CG / 77.21)^THETA_CL_CG) * exp(ETA_CL)
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
      Q  <- THETA_Q  * exp(ETA_Q)
      
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
      THETA_CL    = 6.60,   
      THETA_CL_CG = 0.656,  
      THETA_V1    = 13.3,   
      THETA_V2    = 13.0,   
      THETA_Q     = 16.5    
    ),
    covariates = c("CLCR_CG"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_Q + ETA_V2 ~
        c(
          log((33.7 / 100)^2 + 1),
          0, log((34.1 / 100)^2 + 1),
          0, 0, log((50.8 / 100)^2 + 1),
          0, 0, 0, log((40.1 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.0762)
  ),
  #Cefep_5
  #https://pmc.ncbi.nlm.nih.gov/articles/PMC10353438/
  Bilal_2023 = list(
    # 14 medical-surgical ICU patients (Munich), 344 plasma samples over four
    # consecutive days, median dose 2 g by 30 min infusion q8h to q24h.
    # 4/14 patients on CVVHD, 13/14 septic. Age 62 years (22-94), weight 70 kg
    # (45-120), eCLCR (Cockcroft-Gault) 60.8 mL/min (16.7-217).
    # Final model (Table 3): 2-compartment with linear elimination and a
    # separate clearance for CVVHD patients (CL_D) and non-dialysis patients
    # (CL_ND). IIV on both clearances, V1 and V2; IOV on CL_ND.
    # Covariate (Equation 1, power model): CL_ND = THETA_CL_ND *
    # (eCLCR / median)^0.896, with eCLCR by Cockcroft-Gault in mL/min.
    # NOTE: the normalising constant is not stated explicitly. Equation 1 says
    # "covariate/median" and Table 2 gives a median eCLCR of 60.8 mL/min, while
    # the Discussion quotes the typical CL_ND next to the cohort MEAN eCLCR of
    # 87.1 mL/min. 60.8 is the right one: reproducing the Monte Carlo of
    # Table 4 (90% PTA for 100% fT>2xMIC, MIC 2 mg/L, 20% protein binding,
    # q8h) gives 149/391/783/1385/2264 mg per dose at eCLCR 30/60/90/120/150
    # mL/min, i.e. within 2-7% of the published daily doses of
    # 440/1100/2300/4000/6500 mg divided over three administrations -- and it
    # matches the 1,333 mg q8h quoted in the Discussion for eCLCR 120 mL/min.
    # Normalising to 87.1 would underpredict those doses roughly twofold.
    # Measured creatinine clearance, weight, age and sex did not improve the
    # model and are not part of it.
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      # CVVHD = 1 for patients on continuous venovenous haemodialysis
      if (CVVHD == 1) {
        CL <- THETA_CL_D * exp(ETA_CL_D)
      } else {
        CL <- (THETA_CL_ND * (CLCR_CG / 60.8)^THETA_CL_CG) * exp(ETA_CL_ND + KAPPA_CL_ND)
      }
      V1 <- THETA_V1 * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
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
      THETA_CL_D  = 2.74,    
      THETA_CL_ND = 2.94,    
      THETA_CL_CG = 0.896,   
      THETA_V1    = 8.88,    
      THETA_V2    = 18.7,    
      THETA_Q     = 27.4    
    ),
    covariates = c("CLCR_CG", "CVVHD"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL_D + ETA_CL_ND + ETA_V1 + ETA_V2 ~
        c(
          log((22.7 / 100)^2 + 1),
          0, log((31.3 / 100)^2 + 1),
          0, 0, log((47.3 / 100)^2 + 1),
          0, 0, 0, log((19.3 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.278),
    pi_matrix = lotri::lotri({
      KAPPA_CL_ND ~ log((22.4 / 100)^2 + 1)
    })
  ),
  #Cefep_6
  # https://pubmed.ncbi.nlm.nih.gov/26208469/
  Carlier_2015 = list(
    # 13 septic shock patients on CRRT (CVVHDF or CVVH, AN69 polyacrylonitrile
    # filter), 62 samples pooled from two earlier PK studies. Cefepime 2 g q8h
    # or q12h as a 30 min infusion. Age 59 years (43-70), weight 70 kg (65-75),
    # SOFA 10, blood flow 150 mL/min, ultrafiltration flow rate (UFR)
    # 1750 mL/h (1500-2000), CRRT intensity 36 mL/kg/h.
    # Final model (Table 2): ONE-compartment with zero-order input and combined
    # additive-proportional RUV, BSV on CL and Vd.
    # Covariates (Equations 1 and 2): UFR, normalised to the population median
    # of 1750 mL/h, on BOTH clearance and volume:
    #   TVCL = THETA_CL * (UFR / 1750)   and   TVV = THETA_V * (UFR / 1750).
    # Weight was NOT retained (no correlation with Vd, R2 = 0.0068) and
    # between-occasion variability did not improve the fit.
    # CAUTION: because CL and V scale with UFR identically, k10 is independent
    # of UFR -- changing UFR rescales the whole concentration profile without
    # changing the half-life (~6.3 h). The model is only defined for patients
    # on CRRT; at UFR = 0 both CL and V collapse to zero.
    # Table 2 values are used here; the running text and abstract quote the
    # typical values as 4.4 L/h and 40.9 L with a 30% CV RUV, but the
    # bootstrap CI of the RUV in Table 2 (10.9-27.3%) excludes 30%.
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      # UFR = ultrafiltration flow rate during CRRT, in mL/h
      CL <- (THETA_CL * (UFR / 1750)) * exp(ETA_CL)
      V  <- (THETA_V  * (UFR / 1750)) * exp(ETA_V)
      
      k10 <- CL / V
      
      Cc <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL = 4.5,    
      THETA_V  = 40.8    
    ),
    covariates = c("UFR"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V ~
        c(
          log((37.7 / 100)^2 + 1),
          0, log((21.2 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 3.3, proportional_b = 0.204)
  ),
  
  #Cefep_7
  # https://pubmed.ncbi.nlm.nih.gov/36648390/
  Barreto_2023b = list(
    # Pharmacotherapy 2023;43:1112-1120 -- companion analysis to Barreto_2023
    # (Cefep_3, AAC). 100 critically ill adults at Mayo Clinic (Nov 2018 - Aug
    # 2020), 286 opportunistically collected residual specimens. Cefepime
    # 0.5-2 g every 8-24 h as a 30 min infusion (89% got 2 g, 70% q12h).
    # Age 62 +/- 18 years, weight 86 +/- 28 kg, eGFR_SCr-CysC 86 +/- 40 mL/min,
    # SOFA 7 +/- 4, 82% sepsis. AKI stage II+, kidney replacement therapy and
    # ECMO were exclusions, so this model does not cover RRT/ECMO.
    # CAUTION: same group, same centre, same IRB (#18-004992) and an
    # overlapping cohort as Barreto_2023 -- these patients are the
    # opportunistically sampled arm of that model. The two are NOT independent
    # models of independent data; do not treat them as two votes in a
    # model-averaging or comparative PTA exercise.
    # Final model (supplementary Table S5): ONE-compartment with first-order
    # elimination, combined residual error, BSV on CL and V.
    # Covariates, both log-transformed (Table S4 selected eGFR_SCr-CysC in
    # mL/min over the /1.73 m2 form, BICc 1824 vs 1835):
    #   CL = CLpop * (WT/ref)^0.75 * (eGFR/ref)^0.81   (0.75 fixed, 0.81 estimated)
    #   V  = Vpop  * (WT/ref)^1                        (1 fixed)
    # NOTE: neither the paper nor the supplement states the values the
    # covariates are normalised to. The cohort means (86 kg, 86 mL/min) are
    # used here, on three grounds: (i) reproducing the paper's own %fT>MIC
    # results with the actual dosing mix gives a median %fT>8 of 91% (IQR
    # 63-100) and %fT>32 of 24% (11-44) against the published 85% (66-100) and
    # 21% (9-34), whereas normalising eGFR to 120 mL/min pushes %fT>8 to 99%;
    # (ii) at the cohort mean patient this model then gives CL 5.9 L/h and
    # V 30.1 L against 5.6 L/h and 29.7 L from the AAC model of the same
    # patients, while a 70 kg reference would inflate V to 37 L; (iii) Monolix
    # centres log-transformed covariates on the covariate's central value.
    # If the true references turn out to be different, only the two constants
    # below need changing.
    ppk_model = rxode2::rxode2({
      AUC(0)   <- 0
      centr(0) <- 0
      
      CL <- (THETA_CL * (WT / 86)^THETA_CL_WT * (eGFR_CKD_EPI / 86)^THETA_CL_eGFR) * exp(ETA_CL) # should be eGFR_CKD_EPI_CYSTATIN
      V  <- (THETA_V  * (WT / 86)^THETA_V_WT) * exp(ETA_V)
      
      k10 <- CL / V
      
      Cc <- centr / V
      
      d / dt(centr) <- -k10 * centr
      d / dt(AUC)   <- Cc
    }),
    error_model = function(f, sigma) {
      g <- sigma[1] + sigma[2] * f
      return(g)
    },
    theta = c(
      THETA_CL      = 5.9,    
      THETA_CL_WT   = 0.75,   
      THETA_CL_eGFR = 0.81,   
      THETA_V       = 30.1,   
      THETA_V_WT    = 1       
    ),
    covariates = c("WT", "eGFR_CKD_EPI"),
    omega = lotri::lotri({ # Provided as SD of the random effects (Monolix omega)
      ETA_CL + ETA_V ~
        c(
          0.21^2,
          0, 0.21^2
        )
    }),
    sigma = c(additive_a = 1.47, proportional_b = 0.24)
  ),
  
  #Cefep_8
  #https://pubmed.ncbi.nlm.nih.gov/37071586/
  An_2023 = list(
    # 130 critically ill patients from two ICUs (University of Iowa and
    # Vanderbilt), prospective and opportunistic sampling, 278 plasma samples,
    # Dec 2017 - Dec 2019. Cefepime 1 g (249 doses) or 2 g (518 doses) as
    # intermittent infusions of 0.083-3 h (most 0.5 h). Age 65 years [52-72],
    # TBW 90 kg [70-110], LBW 58 kg [48-67], CLCR_CG by lean body weight
    # 54 mL/min [32-86], 35% AKI, 35% sepsis, 26% septic shock, SOFA 6 [4-9].
    # 12/130 patients on CRRT (10 CVVHDF, 2 CVVH).
    # Final model (Table 2): 2-compartment with zero-order input and
    # first-order elimination, clearance split into renal, non-renal and CRRT
    # components (Equation 2):
    #   CL = THETA_CLR * (CLCR_CG_LBW/54) * (1 - CRRT) + THETA_CLNR
    #        + THETA_CLCRRT * CRRT
    # The model assumes non-renal clearance is the same on and off CRRT, and
    # that CRRT patients have negligible residual renal function -- so the
    # CLCR term is switched off entirely when CRRT = 1. This reproduces the
    # published totals: 2.53 L/h off CRRT at the median CLCR, 2.17 L/h on CRRT.
    # V1 carries total body weight normalised to the median 90 kg with the
    # exponent fixed at 1 (Equation 1); the 0.75 weight exponent mentioned in
    # the methods applies to clearance, where weight was not retained.
    # NOTE: the IIV on V2 is essentially unidentifiable (RSE 285.6%, shrinkage
    # 87.8%) -- it is kept here as published, but consider dropping ETA_V2 if
    # this model is used for MAP estimation.
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      # CRRT = 1 for patients on continuous renal replacement therapy
      # CLCR_CG_LBW = Cockcroft-Gault creatinine clearance using LEAN body weight (mL/min)
      CL <- (THETA_CLR * (CLCR_CG_LBW / 54) * (1 - CRRT) +
             THETA_CLNR + THETA_CLCRRT * CRRT) * exp(ETA_CL)
      V1 <- (THETA_V1 * (WT / 90)) * exp(ETA_V1)
      V2 <- THETA_V2 * exp(ETA_V2)
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
      THETA_CLR    = 2.00,   
      THETA_CLNR   = 0.526,  
      THETA_CLCRRT = 1.64,   
      THETA_V1     = 13.4,   
      THETA_V2     = 7.52,   
      THETA_Q      = 12.0    
    ),
    covariates = c("CLCR_CG_LBW", "WT", "CRRT"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((29.9 / 100)^2 + 1),
          0, log((65.4 / 100)^2 + 1),
          0, 0, log((17.3 / 100)^2 + 1)   # RSE 285.6%, shrinkage 87.8%
        )
    }),
    sigma = c(additive_a = 6.44, proportional_b = 0.224)
  ),
  
  #Cefep_9
  # https://pubmed.ncbi.nlm.nih.gov/22503878/
  Delattre_2012 = list(
    # 88 critically ill patients with severe sepsis or septic shock in four
    # Belgian ICUs, treated with an amikacin loading dose (25 mg/kg) combined
    # with one of four beta-lactams; the cefepime subset is n = 19 (Table 2).
    # Cefepime 2 g every 8 h as a 30 min infusion; FIRST-DOSE pharmacokinetics
    # (samples before and 1, 1.5, 4.5, 6 or 8 and 24 h after start of the
    # infusion), so this model describes the acute phase, not steady state.
    # Final model 1 (Table 2): 2-compartment with first-order elimination.
    # All disposition parameters are allometrically scaled to a 70 kg adult
    # (Equation 5), with the exponent fixed at 1 for volumes and 0.75 for
    # clearances. Renal function enters as a power function (Equations 6-8):
    #   CL = THETA_CL * (WT/70)^0.75 * (CLCR_70 / 100)^THETA_CL_RF
    # CLCR_70 is Cockcroft-Gault computed for a NOMINAL 70 kg subject:
    #   CLCR_70 = (140 - AGE) * 70 / (SCr * 72), x 0.85 if female,
    #   with SCr in mg/dL. The standard weight is used deliberately so that
    #   body weight is not counted twice (it is already in the allometric
    #   term); do NOT substitute an ordinary weight-based Cockcroft-Gault.
    # The reference CLCR is 100 mL/min per 70 kg.
    # No IIV could be estimated on Q ("-" in Table 2), and the residual error
    # for cefepime is proportional only.
    # NOTE: the paper reports two further models (2 and 3) which add amikacin
    # PK parameters or measured amikacin concentrations as covariates on CL
    # and V1. They fit better but require amikacin TDM data, so only model 1
    # is implemented here.
    ppk_model = rxode2::rxode2({
      AUC(0)    <- 0
      centr(0)  <- 0
      periph(0) <- 0
      
      CL <- (THETA_CL * (WT / 70)^0.75 * (CLCR_70 / 100)^THETA_CL_RF) * exp(ETA_CL)
      V1 <- (THETA_V1 * (WT / 70)) * exp(ETA_V1)
      V2 <- (THETA_V2 * (WT / 70)) * exp(ETA_V2)
      Q  <- THETA_Q * (WT / 70)^0.75
      
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
      THETA_CL    = 6.77,   # CL (L/h) for 70 kg at CLCR = 100 mL/min
      THETA_CL_RF = 0.41,   # theta_CL-RF, exponent of the renal function effect
      THETA_V1    = 17.8,   # central volume (L) for 70 kg
      THETA_V2    = 12.2,   # peripheral volume (L) for 70 kg
      THETA_Q     = 11.4    # inter-compartmental clearance (L/h) for 70 kg
    ),
    covariates = c("CLCR_70", "WT"),
    omega = lotri::lotri({ # Provided as CV
      ETA_CL + ETA_V1 + ETA_V2 ~
        c(
          log((47.1 / 100)^2 + 1),
          0, log((39.2 / 100)^2 + 1),
          0, 0, log((47.5 / 100)^2 + 1)
        )
    }),
    sigma = c(additive_a = 0, proportional_b = 0.228)
  )
)
