Fludarabine <- list(
  Langenhorst_CPK2019 = list(
  ppk_model   = rxode2::rxode({
    
    TVCl  = (THETA_Cl+eGFR*(0.78))*(BW/70)^0.75; 
    # THETA_Cl is "Cl70kg-non-renal" from the article
    # 0.78 is the value of Slopepop from the article
    TVV1  = THETA_V1*(BW/70);
    TVV2  = THETA_V2*(BW/70);
    TVV3  = THETA_V3*(BW/70);
    TVQ2  = THETA_Q2*(BW/70)^0.75;
    TVQ3  = THETA_Q3*(BW/70)^0.75;
    
    Cl    = TVCl*exp(ETA_ClQ2Q3 + KAPPA_ClQ2Q3);
    V1    = TVV1*exp(ETA_V1V2V3 + KAPPA_V1V2V3);
    V2    = TVV2*exp(ETA_V1V2V3 + KAPPA_V1V2V3);
    V3    = TVV3*exp(ETA_V1V2V3 + KAPPA_V1V2V3);
    Q2    = TVQ2*exp(ETA_ClQ2Q3 + KAPPA_ClQ2Q3);
    Q3    = TVQ3*exp(ETA_ClQ2Q3 + KAPPA_ClQ2Q3);
    
    ke    = Cl/V1;
    k12   = Q2/V1;
    k21   = Q2/V2;
    k13   = Q3/V1;
    k31   = Q3/V3;
    Cc    = centr/V1;
    
    d/dt(centr)   =  - ke*centr - k12*centr + k21*periph2 - k13*centr + k31*periph3;
    d/dt(periph2) =             + k12*centr - k21*periph2;
    d/dt(periph3) =                                       + k13*centr - k31*periph3;
    d/dt(AUC)     =  Cc;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=3.2,THETA_V1=39,THETA_V2=20,THETA_V3=50,
            THETA_Q2=8.6,THETA_Q3=3.8),
  omega = lotri::lotri({ETA_ClQ2Q3 + ETA_V1V2V3 ~
      c(0.05154826,
        0.00     ,  0.2073393)}),
  pi_matrix = lotri::lotri({KAPPA_ClQ2Q3 + KAPPA_V1V2V3 ~
      c(0.0142973,
        0.00     ,  0.09175843)}),
  covariates  = c("eGFR","BW"),
  sigma       = c(additive_a = 0.0, proportional_b = 0.063))
)
