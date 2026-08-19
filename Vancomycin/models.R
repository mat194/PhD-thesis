models_vancomycin <- list(
  Goti2018 = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/120)^0.8*(0.7^DIAL);
      TVVc  = THETA_Vc*(WT/70)*(0.5^DIAL);
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.5, THETA_Vc=58.4, THETA_Vp=38.4,THETA_Q=6.5),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
        c(0.147,
          0, 0.510,
          0, 0, 0.282,
          0, 0, 0, 0)}),
    sigma = c(additive_a = 3.4, proportional_b = 0.227),
    covariates = c("CLCREAT","WT","DIAL")
  ),
  
  Roberts2011 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/100);
      TVVc  = THETA_Vc*(WT);
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.58, THETA_Vc=1.53),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c(0.375379779^2,
          0, 0.361827976^2)}),
    sigma = c(additive_a = 2.4, proportional_b = 0.227),
    covariates = c("CLCREAT","WT")
  ),
  
  Buelga2015 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl * CLCREAT * 0.06;# CrCL is in L/h not ml/min
      TVVc  = THETA_Vc * WT;
      Cl    = TVCl * exp(ETA_Cl);
      Vc    = TVVc * exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=1.08, THETA_Vc=0.98),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((28.16/ 100)^2 + 1),
           0, log((37.15/ 100)^2 + 1))}),
    sigma = c(additive_a = 3.52, proportional_b = 0),
    covariates = c("CLCREAT","WT")
  ),
  
  Adane2015 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      TVCl  = THETA_Cl * (CLCREAT/125);
      TVVc  = THETA_Vc * WT;
      Cl    = TVCl * exp(ETA_Cl);
      Vc    = TVVc * exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl= 6.54, THETA_Vc= 0.51),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((26.7/ 100)^2 + 1),
           0, log((23.9/ 100)^2 + 1))}),
    sigma = c(additive_a = 0, proportional_b = 0.189),
    covariates = c("CLCREAT","WT")
  ),
  Lopis2006 = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl1 * CLCREAT + THETA_Cl2 * WT ;
      TVVc  = THETA_Vc * WT;
      TVVp  = THETA_Vp * WT;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl1=0.034, THETA_Cl2=0.015, THETA_Vc=0.414, THETA_Vp=1.32,THETA_Q=7.48),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp ~
        c( log((29.2/ 100)^2 + 1),
           0, log((36.4/ 100)^2 + 1),
           0, 0, log((39.8/ 100)^2 + 1))}),
    sigma = c(additive_a = 23.9, proportional_b = 0.185),
    covariates = c("CLCREAT","WT")
  ),
  Okada2018  = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      TVCl  = THETA_Cl*(CLCREAT/113)^0.78;
      TVVc  = THETA_Vc*(WT/59.4)^0.70;
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = TVCl*exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=4.25, THETA_Vc=39.2, THETA_Vp=56.1,THETA_Q=1.95),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp  ~
        c( log((25.2/ 100)^2 + 1),
           0, log((14.2/ 100)^2 + 1),
           0, 0, log((66.9/ 100)^2 + 1))}),
    sigma = c(additive_a = 0, proportional_b = 0.172),
    covariates = c("CLCREAT","WT")
  ),
  Donadello2014 = list(
    ppk_model   = rxode2::rxode2({
      centr(0) = 0;
      if (CRRT == 1){
        TVCl  = 0.6;
      } else {
        TVCl  = 1;
      }
      TVVc  = THETA_Vc;
      TVVp  = THETA_Vp;
      TVQ   = THETA_Q;
      Cl    = THETA_Cl * TVCl *exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      Vp    = TVVp*exp(ETA_Vp);
      Q     = TVQ;
      ke    = Cl/Vc;
      k12   = Q/Vc;
      k21   = Q/Vp;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
      d/dt(periph) = + k12*centr - k21*periph;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl = 3.7, THETA_Vc=31.8, THETA_Vp=57.1,THETA_Q=3.6),
    omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp  ~
        c( log((16.4/ 100)^2 + 1),
           0, log((47/ 100)^2 + 1),
           0, 0, log((101/ 100)^2 + 1))}),
    sigma = c(additive_a = 0, proportional_b = 0.085),
    covariates = c("CRRT")
  ),
  Medellin2017 = list(
    ppk_model   = rxode2::rxode({
      centr(0) = 0;
      if (vm == 0){
        TVCl  = (CLCREAT/100)^0.75;
      } else {
        TVCl  = (CLCREAT/100)^0.75 * 0.8;
      }
      TVVc  = THETA_Vc * WT;
      Cl    = THETA_Cl * TVCl *exp(ETA_Cl);
      Vc    = TVVc*exp(ETA_Vc);
      ke    = Cl/Vc;
      Cc    = centr/Vc;
      d/dt(centr)  = - ke*centr;
      d/dt(AUC)    = Cc;
    }),
    error_model = function(f,sigma){
      g <- sigma[1] + sigma[2]*f
      return(g)
    },
    theta = c(THETA_Cl=2.86, THETA_Vc=1.03),
    omega = lotri::lotri({ETA_Cl + ETA_Vc ~
        c( log((28.4/ 100)^2 + 1),
           0, log((49.1/ 100)^2 + 1))}),
    sigma = c(additive_a = 4.3, proportional_b = 0),
    covariates = c("CLCREAT","WT","vm")
  )
)
# Thomson2009 = list(
#   ppk_model   = rxode2::rxode2({
#     centr(0)  = 0;
#     periph(0) = 0;
#     AUC(0)    = 0;
#     TVCl  = THETA_Cl * (1 + 0.0154 * (CLCREAT - 66));
#     TVVc  = THETA_Vc;
#     TVVp  = THETA_Vp;
#     TVQ   = THETA_Q;
#     Cl    = TVCl*exp(ETA_Cl);
#     Vc    = TVVc*exp(ETA_Vc);
#     Vp    = TVVp*exp(ETA_Vp);
#     Q     = TVQ*exp(ETA_Q);
#     ke    = Cl/Vc;
#     k12   = Q/Vc;
#     k21   = Q/Vp;
#     Cc    = centr/Vc;
#     d/dt(centr)  = - ke*centr - k12*centr + k21*periph;
#     d/dt(periph) = + k12*centr - k21*periph;
#     d/dt(AUC)    = Cc;
#   }),
#   error_model = function(f,sigma){
#     g <- sigma[1] + sigma[2]*f
#     return(g)
#   },
#   theta = c(THETA_Cl=2.99, THETA_Vc=0.732, THETA_Vp=0.675,THETA_Q=2.28),
#   omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Q ~
#       c( log((27/ 100)^2 + 1),
#          0, log((15/ 100)^2 + 1),
#          0, 0, log((130/ 100)^2 + 1),
#          0, 0, 0, log((49/ 100)^2 + 1))}),
#   sigma = c(additive_a = 1.6, proportional_b = 0.15),
#   covariates = c("CLCREAT")
# )
