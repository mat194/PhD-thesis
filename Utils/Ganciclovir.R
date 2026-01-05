# models
mod_ganciclovir_Caldes_AAC2009 <- function() {
  ini({
    THETA_cl  <- 7.49
    THETA_v1  <- 31.90
    THETA_cld <- 10.20
    THETA_v2  <- 32.0
    THETA_ka  <- 0.895
    THETA_baf <- 0.825
    ETA_cl ~ 0.107
    ETA_v1 ~ 0.227
    ETA_ka ~ 0.464
    ETA_baf ~ 0.049
    add.sd <-0.465
    prop.sd <- 0.143
  })
  model({
    TVcl  = THETA_cl*(ClCr/57);
    TVv1  = THETA_v1;
    TVcld = THETA_cld;
    TVv2  = THETA_v2;
    TVka  = THETA_ka;
    TVbaf = THETA_baf;
    
    cl  = TVcl*exp(ETA_cl);
    v1  = TVv1*exp(ETA_v1);
    cld = TVcld;
    v2  = TVv2;
    ka  = TVka*exp(ETA_ka);
    baf = TVbaf*exp(ETA_baf);
    
    k10 = cl/v1;
    k12 = cld / v1;
    k21 = cld / v2;
    Cc = centr/v1;
    
    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k10*centr - k12*centr + k21*periph;
    d/dt(periph) =                         k12*centr - k21*periph;
    d/dt(AUC)    = Cc;
    
    f(depot)=baf*0.72;
    alag(depot)=0.382;
    
    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}

mod_ganciclovir_Lalagkas_CPK2023 <- list(
  ppk_model   = rxode2::rxode({
    
    TVCl  = THETA_Cl*(CKDepi/55)^0.817*(WT/70)^0.75;
    TVVc  = THETA_Vc*(WT/70);
    TVVp  = THETA_Vp*(WT/70);
    TVCld = THETA_Cld*(WT/70)^0.75;
    TVKa  = THETA_Ka;
    TVbaF = THETA_baF;
    
    Cl    = TVCl*exp(ETA_Cl + KaPPA_Cl);
    Vc    = TVVc*exp(ETA_Vc);
    Vp    = TVVp*exp(ETA_Vp);
    Cld   = TVCld;
    Ka    = TVKa *exp(ETA_Ka);
    baF   = TVbaF *exp(ETA_baF);
    
    ke    = Cl/Vc;
    k12   = Cld/Vc;
    k21   = Cld/Vp;
    Cc    = centr/Vc;
    
    d/dt(depot)  = -Ka*depot
    d/dt(centr)  =  Ka*depot - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =                      + k12*centr - k21*periph;
    d/dt(AUC)    =  Cc;
    
    f(depot)=baF;
    alag(depot)=0.331;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=6.93,THETA_Vc=43.1,THETA_Vp=219,THETA_Cld=9.23,THETA_Ka=0.766,THETA_baF=0.699),
  omega = lotri::lotri({ETA_Cl + ETA_Vc + ETA_Vp + ETA_Ka + ETA_baF ~
      c(0.0893,
        0.00     ,  0.130,
        0.00     ,  0.00      ,  1.07,
        0.00     ,  0.00      ,  0.00     ,  0.209,
        0.00     ,  0.00      ,  0.00     ,  0.00    ,  0.0275)}),
  pi_matrix = lotri::lotri({KaPPA_Cl ~
      c(0.025)}),
  covariates  = c("CKDepi","WT"),
  sigma       = c(additive_a = 0.237, proportional_b = 0.282)) #proportional_b = 0.282

mod_ganciclovir_Padulles_TDM2014 <- function() {
  ini({
    THETA_cl  <- 7.49
    THETA_v1  <- 31.90
    THETA_cld <- 10.20
    THETA_v2  <- 32.0
    THETA_ka  <- 0.895
    THETA_baf <- 0.825
    ETA_cl ~ 0.107
    ETA_v1 ~ 0.227
    ETA_ka ~ 0.464
    ETA_baf ~ 1.54
    add.sd <- 0.465
    prop.sd <- 0.143
  })
  model({
    TVcl  = THETA_cl*(ClCr/57);
    TVv1  = THETA_v1;
    TVcld = THETA_cld;
    TVv2  = THETA_v2;
    TVka  = THETA_ka;
    TVbaf = THETA_baf;
    
    cl  = TVcl*exp(ETA_cl);
    v1  = TVv1*exp(ETA_v1);
    cld = TVcld;
    v2  = TVv2;
    ka  = TVka*exp(ETA_ka);
    baf = TVbaf*exp(ETA_baf);
    
    k10 = cl/v1;
    k12 = cld / v1;
    k21 = cld / v2;
    Cc = centr/v1;
    
    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k10*centr - k12*centr + k21*periph;
    d/dt(periph) =                         k12*centr - k21*periph;
    d/dt(AUC)    = Cc;
    
    f(depot)=baf;
    alag(depot)=0.382;
    
    Cc ~ add(add.sd) + prop(prop.sd) + combined1()
  })
}

mod_valganciclovir_Chen_JCP2020 <- function() {
  ini({
    THETA_cl <- 7.09
    THETA_v2 <- 10.8
    THETA_q  <- 3.96
    THETA_v3 <- 174
    THETA_ka <- 0.23
    ETA_cl ~ 0.0713751
    ETA_v2 ~ 1.20624
    ETA_q  ~ 0.3351578
    ETA_v3 ~ 0.7630929
    prop.sd <- 0.429
  })
  model({
    TVcl = THETA_cl*(1+ClCr/68.3*1.08);
    TVv2 = THETA_v2;
    TVq  = THETA_q;
    TVv3 = THETA_v3;
    TVka = THETA_ka;
    
    cl = TVcl*exp(ETA_cl);
    v2 = TVv2*exp(ETA_v2);
    q  = TVq*exp(ETA_q);
    v3 = TVv3*exp(ETA_v3);
    ka = TVka;
    
    k20 = cl/v2;
    k23 = q / v2;
    k32 = q / v3;
    Cc  = centr/v2;
    
    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k20*centr - k23*centr + k32*periph;
    d/dt(periph) =                         k23*centr - k32*periph;
    d/dt(AUC)    = Cc;
    
    alag(depot)=0.93;
    
    Cc ~ prop(prop.sd)
  })
}

mod_valganciclovir_Vezina_BCJP2014 <- function() {
  ini({
    THETA_cl  <- 14.5
    THETA_v2  <- 87.5
    THETA_q   <- 4.80
    THETA_v3  <- 42.6
    THETA_ka  <- 3.0
    ETA_cl ~ 0.1064
    ETA_ka ~ 0.0001
    prop.sd <- 0.327
  })
  model({
    TVcl = THETA_cl*((ClCr/60)*(70/WT))^0.492*(WT/70)^0.75;
    TVv2 = THETA_v2*(70/WT);
    TVq  = THETA_q*(70/WT)^0.75;
    TVv3 = THETA_v3*(70/WT);
    TVka = THETA_ka;
    
    cl = TVcl*exp(ETA_cl);
    v2 = TVv2;
    q  = TVq;
    v3 = TVv3;
    ka = TVka*exp(ETA_ka);
    
    k20 = cl/v2;
    k23 = q/v2;
    k32 = q/v3;
    Cc  = centr/v2;
    
    d/dt(depot)  = -ka*depot
    d/dt(centr)  =  ka*depot - k20*centr - k23*centr + k32*periph;
    d/dt(periph) =                         k23*centr - k32*periph;
    d/dt(AUC)    = Cc;
    
    alag(depot)=0.5;
    
    Cc ~ prop(prop.sd)
  })
}

mod_ganciclovir_Perrottet_TDM2009 <- list(
  ppk_model   = rxode2::rxode({
    
    TVCl  = THETA_Cl*MDRD*(1.21)^FEMALE; # != THETA_Cl according to GRAFT TYPE
    # kidney=1.68
    # lung/liver=1.17
    # heart=0.86
    TVV1  = THETA_V1*(WT/70)*(0.78)^FEMALE;
    TVV2  = THETA_V2;
    TVQ   = THETA_Q;
    TVKa  = THETA_Ka;
    
    Cl    = TVCl*exp(ETA_Cl + KaPPA_Cl);
    V1    = TVV1*exp(ETA_V1);
    V2    = TVV2;
    Q     = TVQ;
    Ka    = TVKa;
    
    ke    = Cl/V1;
    k12   = Q/V1;
    k21   = Q/V2;
    Cc    = centr/V1;
    
    d/dt(depot)  = -Ka*depot
    d/dt(centr)  =  Ka*depot - ke*centr - k12*centr + k21*periph;
    d/dt(periph) =                      + k12*centr - k21*periph;
    d/dt(AUC)    =  Cc;
    
    f(depot)=0.6;
  }),
  error_model = function(f,sigma){
    g <- sigma[1] + sigma[2]*f
    return(g)
  },
  theta = c(THETA_Cl=1.68,THETA_V1=24,THETA_V2=22,THETA_Q=4.1,THETA_Ka=0.56),
  omega = lotri::lotri({ETA_Cl + ETA_V1 ~
      c(0.06541314,
        0.00     ,  0.03922071)}),
  pi_matrix = lotri::lotri({KaPPA_Cl ~
      c(0.04315527)}),
  covariates  = c("MDRD","WT","FEMALE"),
  sigma       = c(additive_a = 0.0, proportional_b = 0.21))

# list of models
model_list <- list(mod_ganciclovir_Lalagkas_CPK2023,
                   mod_valganciclovir_Vezina_BCJP2014,
                   mod_valganciclovir_Chen_JCP2020,
                   mod_ganciclovir_Padulles_TDM2014,
                   mod_ganciclovir_Caldes_AAC2009,
                   mod_ganciclovir_Perrottet_TDM2009)

model_list_names <- list("Lalagkas2023",
                         "Vezina2014",
                         "Chen2020",
                         "Padulles2014",
                         "Caldes2009",
                         "Perrottet2009")

m_out <- vector("list", length(model_list))

# libraries
library(posologyr)
library(readr)
library(data.table)

# import patient
patient <- read_csv("~/Documents/recherche/outils R/applications posologyr/ganciclovir/bouet-leboeuf.csv")

# map estimation for each model
for (i in 1:length(model_list)){
  m_i <- poso_estim_map(patient,model_list[[i]],return_model = TRUE,return_ofv = TRUE)
  
  m_out[[i]] <- cbind(data.table::data.table(m_i$model)[,c("time","Cc","AUC")],
                      ofv = m_i$ofv,
                      LL = exp(-0.5 * m_i$ofv),
                      model_name = model_list_names[[i]],
                      model_number = i)
}

m_dt_out <- data.table::data.table(do.call(rbind.data.frame,m_out))

# calculate weights
m_dt_out[,weight := LL/sum(unique(LL))]

# model averaging
m_dt_out[,MoA_IPRED:=sum(Cc*weight),by=time]

m_dt_out[,MoA_AUC:=sum(AUC*weight),by=time]
# plots
indiv_obs           <- patient[,c("DV","TIME")]
names(indiv_obs)    <- c("Cc","time")

ggplot2::ggplot(m_dt_out,ggplot2::aes(x=time,y=Cc)) +
  ggplot2::geom_line(ggplot2::aes(color=model_name)) +
  ggplot2::ylab("Plasma concentration") +
  ggplot2::geom_point(data=indiv_obs, size= 3, na.rm=TRUE) +
  ggplot2::xlim(0,24) +
  ggplot2::theme_bw() +
  ggplot2::facet_grid(model_name ~ .)

ggplot2::ggplot(m_dt_out,ggplot2::aes(x=time,y=Cc)) +
  ggplot2::geom_line(ggplot2::aes(group=model_name,color=model_name)) +
  ggplot2::geom_line(ggplot2::aes(y=MoA_IPRED),size=1) +
  ggplot2::ylab("Plasma concentration") +
  ggplot2::geom_point(data=indiv_obs, size= 3, na.rm=TRUE) +
  ggplot2::xlim(0,24) +
  ggplot2::theme_bw()

ggplot2::ggplot(m_dt_out[time==0],ggplot2::aes(x=time,y=weight,fill=model_name)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::ylim(0,1.01) +
  ggplot2::xlab("") +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank())

# all AUC
m_AUC24 <- as.data.table(list(model_name=unique(m_dt_out$model_name),
                              AUC_24=m_dt_out[round(time,1) == 24,AUC,
                                              by=model_name]$AUC - m_dt_out[round(time,1) == 0,AUC,
                                                                            by=model_name]$AUC))
# MoA_AUC24 / min AUC24 / max AUC24
unique(m_dt_out[round(time,1) == 24,MoA_AUC] - m_dt_out[round(time,1) == 0,MoA_AUC])
m_AUC24[AUC_24 == min(AUC_24),]
m_AUC24[AUC_24 == max(AUC_24),]
unique(m_dt_out[round(time,1) == 48,MoA_AUC] - m_dt_out[round(time,1) == 24,MoA_AUC])
# full table
m_AUC24