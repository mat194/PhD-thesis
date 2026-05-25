# Conversion serum_creatinine_mg_dl to serum_micromol_l
convert_creatinine <- function(serum_creatinine_mg_dl) {
  serum_micromol_l <- serum_creatinine_mg_dl * 88.4
  return(serum_micromol_l)
}

# Calculate BSA (Du Bois formula)
calculate_bsa <- function(weight_kg, height_cm) {
  bsa <- 0.007184 * (weight_kg^0.425) * (height_cm^0.725)
  return(bsa)
}
# Calculate ABW
calculate_abw <- function(height_cm, weight, gender) {
  height_in <- (height_cm / 100) * 39.37
  if (gender == "male") {
    IBW <- 50 + 2.3 * (height_in - 60)
  } else {
    IBW <- 45.5 + 2.3 * (height_in - 60)
  }
  ABW <- IBW + 0.4 * (weight - IBW)
  return(ABW)
}

# Calculate LBW, James Formula
calculate_lbw <- function(height_cm, weight_kg, gender){
  if (gender =="male"){
    lbw <- 1.1 * weight_kg - 128 * (weight_kg / height_cm)^2
  } else {
    lbw <- 1.07 * weight_kg - 148 * (weight_kg / height_cm)^2
  }
  return(lbw)
}

# Ideal Body Weight Devine formula
calculate_ibw <- function(height_cm, gender) {
  height_in <- (height_cm / 100) * 39.37
  if (gender == "male") {
    50   + 2.3 * (height_in - 60)
  } else {
    45.5 + 2.3 * (height_in - 60)
  }
}

# Calculate CRCL_CG
calculate_crcl <- function(age, weight_kg, serum_creatinine_mg_dl, gender) {
  CrCl <- (140 - age) * weight_kg / (72 * serum_creatinine_mg_dl)
  if (gender == "female") {
    CrCl <- CrCl * 0.85
  }
  return(CrCl)
}
# Calculate eGFR -2021
calculate_egfr <- function(age, serum_creatinine_mg_dl, gender) {
  if (gender == "male") {
    kappa <- 0.9
    alpha <- -0.302
    min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
    max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
    egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age)
  } else {
    kappa <- 0.7
    alpha <- -0.241
    min_sc_k <- min(serum_creatinine_mg_dl / kappa, 1)
    max_sc_k <- max(serum_creatinine_mg_dl / kappa, 1)
    egfr <- 142 * (min_sc_k^alpha) * (max_sc_k^-1.200) * (0.9938^age) * 1.012
  }
  return(egfr)
}

# eGFR – MDRD, race
calculate_egfr_mdrd <- function(age, serum_creatinine_mg_dl,
                                gender, race = "non-african") {
  egfr <- 186.3 * serum_creatinine_mg_dl^(-1.154) * age^(-0.203)
  if (gender == "female")  egfr <- egfr * 0.742
  if (race   == "african") egfr <- egfr * 1.21
  egfr
}


# Calculate CRCL_J (according to Jeliffe equation)
calculate_crcl_jelliffe <- function(age, weight, serum_creatinine, gender) {
  if (gender == "male") {
    crcl <- (98 - 0.8 * (age - 20)) / serum_creatinine * (0.9 / (weight / 70))
  } else {
    crcl <- (98 - 0.8 * (age - 20)) / serum_creatinine * (0.85 * 0.9 / (weight / 70))
  }
  return(crcl)
}

# Calculate CKD_EPI_CYSTATIN_C 2021
calculate_egfr_cystatin_c <- function(age, serum_creatinine_mg_dl, gender, serum_cystatin_c_mg_l) {
  if (gender == "male") {
    if (serum_creatinine_mg_dl <= 0.9 & serum_cystatin_c_mg_l <= 0.8) {
      A <- 0.9
      B <- -0.144
      C <- 0.8
      D <- -0.323
    } else if (serum_creatinine_mg_dl <= 0.9 & serum_cystatin_c_mg_l > 0.8) {
      A <- 0.9
      B <- -0.144
      C <- 0.8
      D <- -0.778
    } else if (serum_creatinine_mg_dl > 0.9 & serum_cystatin_c_mg_l <= 0.8) {
      A <- 0.9
      B <- -0.544
      C <- 0.8
      D <- -0.323
    } else {
      A <- 0.9
      B <- -0.544
      C <- 0.8
      D <- -0.778
    }
    egfr <- 135 * (serum_creatinine_mg_dl / A)^B * (serum_cystatin_c_mg_l / C)^D * 0.9961^age 
  } else {
    if (serum_creatinine_mg_dl <= 0.7 & serum_cystatin_c_mg_l <= 0.8) {
      A <- 0.7
      B <- -0.219
      C <- 0.8
      D <- -0323
    } else if (serum_creatinine_mg_dl <= 0.7 & serum_cystatin_c_mg_l > 0.8) {
      A <- 0.7
      B <- -0.219
      C <- 0.8
      D <- -0.778
    } else if (serum_creatinine_mg_dl > 0.7 & serum_cystatin_c_mg_l <= 0.8) {
      A <- 0.7
      B <- -0.544
      C <- 0.8
      D <- -0.323
    } else {
      A <- 0.7
      B <- -0.544
      C <- 0.8
      D <- -0.778
    }
    egfr <- 135 * (serum_creatinine_mg_dl / A)^B * (serum_cystatin_c_mg_l / C)^D * 0.9961^age * 0.963
  }
  return(egfr)
}

# Sc is a patient-level sieving coefficient input
calculate_cl_cvvh <- function(Qf_ml_h, Sc = 0.8, CF = 1) {
  (Qf_ml_h / 1000) * Sc * CF
}