library(purrr)

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


# Metrics
# All four metrics share the symmetric denominator (predicted + observed)/2, so
# rBias/rRMSE/MDPE/MDAPE are on one common scale.
#   rBias  = mean relative prediction error (%)
#   rRMSE  = root mean square relative prediction error (%)
#   MDPE   = median relative prediction error (%)      - accuracy/bias
#   MDAPE  = median absolute relative prediction error (%) - precision
# F20/F30 are the exception: they use the conventional |(pred-obs)/obs| so they
# stay comparable with published values (see the note in calculate_single).
# Set id_col to a subject identifier to bootstrap by patient (cluster bootstrap)
# instead of by observation. Required whenever a subject contributes more than
# one sample, otherwise the confidence intervals come out too narrow.
calculate_metrics <- function(data, n_bootstrap = 0, conf_level = 0.95,
                              id_col = NULL) {

  keep <- is.finite(data$predicted) & is.finite(data$observed) &
          (data$predicted + data$observed) != 0
  if (any(!keep)) {
    warning(sum(!keep), " row(s) dropped (missing/non-finite values, ",
            "or predicted + observed == 0).")
  }
  data <- data[keep, , drop = FALSE]
  N <- nrow(data)
  if (N == 0) stop("No usable rows in `data`.")

  # Metrics for a single dataset
  calculate_single <- function(d) {
    rPE <- (d$predicted - d$observed) / ((d$predicted + d$observed) / 2)

    # F20 / F30: share of predictions within +/-20% and +/-30% of the observation.
    # NOTE these two deliberately use the CONVENTIONAL denominator (the observed
    # value), not the symmetric one used above. F20/F30 are a named convention
    # whose whole value is comparability with published TDM papers, and every one
    # of those defines them as |(pred - obs)/obs|. Scoring them symmetrically
    # would shift the cut-off (a +/-20% classic band equals +/-18.2% symmetric)
    # and make the numbers silently non-comparable.
    cPE <- abs((d$predicted - d$observed) / d$observed)

    c(rBias = mean(rPE)         * 100,
      rRMSE = sqrt(mean(rPE^2)) * 100,
      MDPE  = median(rPE)       * 100,
      MDAPE = median(abs(rPE))  * 100,
      F20   = mean(cPE <= 0.20) * 100,
      F30   = mean(cPE <= 0.30) * 100)
  }

  original_metrics <- calculate_single(data)

  # Bootstrapping for confidence intervals if n_bootstrap > 0
  if (n_bootstrap > 0 && N > 1) {

    if (is.null(id_col)) {
      # Observation-level bootstrap: assumes independent rows.
      resample <- function() data[sample.int(N, N, replace = TRUE), , drop = FALSE]
    } else {
      # Cluster bootstrap: resample subjects, keeping each subject's rows together.
      if (!id_col %in% names(data)) {
        stop("`id_col` not found in `data`: ", id_col)
      }
      idx_by_id <- split(seq_len(N), data[[id_col]])
      n_id      <- length(idx_by_id)
      resample  <- function() {
        data[unlist(idx_by_id[sample.int(n_id, n_id, replace = TRUE)],
                    use.names = FALSE), , drop = FALSE]
      }
    }

    bootstrap_metrics <- replicate(n_bootstrap, calculate_single(resample()))

    # Percentile confidence intervals: 2 x 4 matrix, columns named by metric
    probs <- c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2)
    ci    <- apply(bootstrap_metrics, 1, quantile, probs = probs, na.rm = TRUE)

    results <- list()
    for (m in names(original_metrics)) {
      results[[m]]                <- unname(original_metrics[[m]])
      results[[paste0(m, "_CI")]] <- c(lower = unname(ci[1, m]),
                                       upper = unname(ci[2, m]))
    }
  } else {
    results <- as.list(original_metrics)
  }

  return(results)
}


convert_creatinine_safe <- possibly(convert_creatinine, otherwise = NA_real_)
calculate_bsa_safe <- possibly(calculate_bsa, otherwise = NA_real_)
calculate_ibw_safe <- possibly(calculate_ibw, otherwise = NA_real_)
calculate_lbw_safe <- possibly(calculate_lbw, otherwise = NA_real_)
calculate_abw_safe <- possibly(calculate_abw, otherwise = NA_real_)
calculate_crcl_safe <- possibly(calculate_crcl, otherwise = NA_real_)
calculate_crcl_jelliffe_safe <- possibly(calculate_crcl_jelliffe, otherwise = NA_real_)
calculate_egfr_safe <- possibly(calculate_egfr, otherwise = NA_real_)
calculate_egfr_mdrd_safe <- possibly(calculate_egfr_mdrd, otherwise = NA_real_)
calculate_egfr_cystatin_c_safe <- possibly(calculate_egfr_cystatin_c, otherwise = NA_real_)