# Multivariate Exact Discrepancy (MED)
# =============================================================================
# Baklouti S, Comets E, Gandia P, Concordet D. Multivariate Exact Discrepancy:
# A New Tool for PK/PD Model Evaluation. Clin Pharmacokinet 2023;62:1599-1609.
# doi:10.1007/s40262-023-01296-6
#
# WHAT IT COMPUTES
# For one patient the model does not predict a concentration, it predicts a
# JOINT distribution over the whole vector of that patient's concentrations,
# (Y_i1, ..., Y_in). MED asks one question of that joint distribution:
#
#     what fraction of the profiles the model can produce are LESS typical
#     than the profile we actually observed?
#
# "Less typical" is measured by the joint density f: a simulated vector sitting
# where f is low is an unusual thing for the model to produce. Writing
# lambda = f(y_obs), the MED for that patient is
#
#     MED_i = P( f(Y) <= lambda )  =  the share of simulated vectors landing in
#                                     a region of density no higher than the
#                                     observed vector's.
#
# MED_i near 1  -> the observed profile sits in the heart of the model.
# MED_i near 0  -> the model can produce this profile only very rarely.
#
# WHY THIS IS WORTH HAVING ON TOP OF VPC / NPDE / rBias
# The theoretical result behind it (paper Sect. 2.2) is that for ANY random
# vector X with density f, the variable F(lambda(X)) is exactly Uniform(0,1) --
# whatever the dimension of X, and with NO Gaussian assumption anywhere. So:
#
#   * ACROSS PATIENTS the MEDs must be U(0,1) if the model is right. One
#     Kolmogorov-Smirnov test of N numbers against U(0,1) is then a global
#     goodness-of-fit test for that model on this cohort.
#   * WITHIN A PATIENT the dependence between that patient's own samples is
#     used rather than assumed away. This is the gap the paper opens against
#     NPDE: NPDE decorrelates the prediction discrepancies and then calls them
#     independent, but decorrelation implies independence only for Gaussian
#     vectors, and a nonlinear mixed-effects model does not produce a Gaussian
#     concentration vector. MED never takes that step.
#   * ONE NUMBER PER PATIENT, on a probability scale, so a patient can be
#     flagged as an outlier as a whole KINETIC PROFILE rather than one
#     suspicious point at a time.
#
# HOW IT RELATES TO THE ETA / MAHALANOBIS SECTION OF THE ANALYSIS SCRIPTS
# The d2 = eta' Omega^-1 eta machinery asks the same shape of question -- how
# improbable is this patient under this model -- but in PARAMETER space and
# conditionally on MAP estimates. That makes it shrinkage-contaminated, which
# is why that section needs an explicit calibration block. MED is the
# OBSERVATION-space analogue and contains no estimation step at all: the etas
# are integrated out by simulation, so there is nothing to shrink and nothing
# to calibrate. The two are complementary rather than redundant -- d2 says
# WHICH PARAMETER is odd, MED says HOW ODD THE PATIENT IS, honestly.
#
# TWO ESTIMATORS ARE PROVIDED
#   method = "kde"         faithful to the paper. Estimates f with a kernel
#                          density estimate over the K simulated vectors
#                          (package ks), then ranks. No shape assumption.
#   method = "mahalanobis" cheap elliptical approximation. Replaces the density
#                          ranking with a Mahalanobis-distance ranking computed
#                          from the SAME simulated sample. Equal to MED exactly
#                          whenever the density contours are ellipsoids, which
#                          log-concentrations are usually close to. Needs no
#                          extra package and is far faster.
#   Both are valid uniform-under-the-null statistics -- any scalar ordering of
#   profiles gives a uniform p-value under the null -- they differ only in how
#   sharply they order. Screen with "mahalanobis", report with "kde".
#
# =============================================================================


# --- one patient's predictive sample -----------------------------------------
# Returns a K x n_i matrix: one row per simulated kinetic profile, one column
# per observation time, on the concentration scale, WITH residual error.
#
# THREE THINGS poso_simu_pop DOES NOT DO, AND THIS FUNCTION HAS TO
#   1. It returns f (the model prediction), not Y. MED needs the full
#      predictive distribution, so the residual error is added here, following
#      posologyr's convention Y = f + g(f, sigma) * eps with eps ~ N(0, 1) and
#      g supplied by the model as error_model().
#   2. It fixes every KAPPA at zero -- read its source: kappa_mat is a matrix
#      of 0s. For a model carrying an IOV prior the simulated spread is
#      therefore TOO NARROW, which pushes MED toward 0 and would make the model
#      look misspecified when it is not. Such models are refused rather than
#      silently mis-scored; allow_iov = TRUE overrides and accepts a result
#      that is conservative (biased against the model).
#   3. It reads covariates from the FIRST row of the patient's data only, so a
#      time-varying covariate is held at baseline for the whole profile. That
#      is a property of posologyr, not of MED, and it is warned about.
#
# nonneg: WHAT TO DO ABOUT SIMULATED CONCENTRATIONS BELOW ZERO, AND WHY IT IS
# NOT A DETAIL. posologyr writes every residual error additively on f, so even
# a purely proportional model, Y = f * (1 + CV*eps), produces Y <= 0 whenever
# eps < -1/CV. At CV = 33% that is P(eps < -3.01) = 0.13% of all simulated
# points. Left alone, every profile containing one gets floored at the log
# guard, lands in the extreme tail of the density, and is therefore counted as
# "less typical" than essentially any observed profile. MED then has a HARD
# FLOOR of roughly n * 0.13%, which for a 3-sample patient is 0.004 -- three
# times LL_alpha at N = 39 (0.0013). The outlier test would be structurally
# incapable of firing, and silently so.
#   "truncate" (default) draws the residual error from a normal truncated at
#   Y > 0 instead, by inverse-CDF, which costs nothing and removes the
#   artefact. It is a mild change to the published error model and should be
#   reported as such; it is also the assumption the published model implies,
#   since a negative concentration is not a thing the drug can do.
#   "none" reproduces posologyr's untruncated behaviour, for comparison.
med_simulate_patient <- function(dat_i, prior_model, K = 1e5,
                                 allow_iov = FALSE,
                                 nonneg = c("truncate", "none")) {

  nonneg <- match.arg(nonneg)

  if (!is.null(prior_model$pi_matrix) && !allow_iov) {
    stop("Model carries an IOV prior (pi_matrix) and poso_simu_pop fixes all ",
         "KAPPAs at 0, so the simulated predictive distribution would be too ",
         "narrow and MED biased downward. Set allow_iov = TRUE to accept a ",
         "conservative result, or simulate the KAPPAs yourself.")
  }

  obs_rows <- dat_i[dat_i$EVID %in% c(0, 2) & !is.na(dat_i$DV), , drop = FALSE]
  obs_rows <- obs_rows[order(obs_rows$TIME), , drop = FALSE]
  if (nrow(obs_rows) < 1L) stop("No usable observation rows for this patient.")

  covs <- prior_model$covariates
  if (length(covs)) {
    varying <- vapply(covs, function(cc) {
      v <- dat_i[[cc]]
      !is.null(v) && length(unique(v[!is.na(v)])) > 1L
    }, logical(1))
    if (any(varying)) {
      warning("Time-varying covariate(s) held at baseline by poso_simu_pop: ",
              paste(covs[varying], collapse = ", "))
    }
  }

  # poso_simu_pop needs the observations VISIBLE as sampling points, so any row
  # masked to EVID = 2 upstream is put back to EVID = 0 here. This leaks
  # nothing: the simulation is drawn from the population prior and never looks
  # at DV. Masking exists to keep MAP fits honest; MED has no fit to keep
  # honest, it is a prior-predictive check by construction.
  sim_dat <- dat_i
  sim_dat$EVID[sim_dat$EVID == 2 & !is.na(sim_dat$DV)] <- 0

  s  <- posologyr::poso_simu_pop(sim_dat, prior_model, n_simul = K)
  md <- as.data.frame(s$model)

  md <- md[md$time %in% obs_rows$TIME, c("sim.id", "time", "Cc"), drop = FALSE]
  md <- md[order(md$sim.id, md$time), , drop = FALSE]

  times <- sort(unique(md$time))
  if (length(times) != nrow(obs_rows)) {
    stop("Solver returned ", length(times), " of the ", nrow(obs_rows),
         " observation times. Duplicated TIMEs in the patient's data are the ",
         "usual cause.")
  }

  f <- matrix(md$Cc, nrow = K, ncol = length(times), byrow = TRUE)

  g <- prior_model$error_model(f, prior_model$sigma)
  if (nonneg == "truncate") {
    # Y = f + g*eps > 0  <=>  eps > -f/g. Inverse-CDF draw from the standard
    # normal truncated to that half-line: one runif and one qnorm, vectorised
    # over the whole K x n matrix, so no rejection loop.
    lo  <- stats::pnorm(-f / g)
    eps <- stats::qnorm(lo + stats::runif(length(f)) * (1 - lo))
  } else {
    eps <- matrix(stats::rnorm(length(f)), nrow = nrow(f))
  }
  Y <- f + g * eps

  attr(Y, "times") <- times
  attr(Y, "y_obs") <- obs_rows$DV[match(times, obs_rows$TIME)]
  Y
}


# --- MED from a simulated sample ---------------------------------------------
# Ysim   K x n matrix of simulated observations (rows = profiles)
# y_obs  length-n observed vector
#
# ON THE SCALE (not cosmetic)
# MED is NOT invariant to a nonlinear change of variable: a monotone transform
# reshapes the density and so reorders which profiles count as "less typical".
# The paper works in log-concentration -- its Figs. 3 and 4 are on log axes --
# and so does the default here, because concentration vectors are strongly
# right-skewed on the natural scale and one global bandwidth matrix cannot
# follow that. Whatever is chosen must be held FIXED across every model being
# compared, otherwise the comparison is not between models.
#
# floor  guards the log. Additive residual error can push a simulated
#        concentration to <= 0; those are censored at `floor` rather than
#        dropped, because dropping them would quietly remove exactly the
#        low-density profiles that set the tail MED cares about. The censored
#        fraction is returned so it can be checked against the assay LLOQ.
med_from_sample <- function(Ysim, y_obs = attr(Ysim, "y_obs"),
                            method = c("kde", "mahalanobis"),
                            scale  = c("log", "identity"),
                            floor  = 0.01,
                            bgridsize = 101L) {

  method <- match.arg(method)
  scale  <- match.arg(scale)
  stopifnot(is.matrix(Ysim), length(y_obs) == ncol(Ysim))

  frac_censored <- mean(Ysim < floor)
  if (scale == "log") {
    Ysim  <- log(pmax(Ysim,  floor))
    y_obs <- log(pmax(y_obs, floor))
  }

  y_obs <- matrix(as.numeric(y_obs), nrow = 1L)
  d <- ncol(Ysim)

  if (method == "mahalanobis") {
    mu <- colMeans(Ysim)
    S  <- stats::cov(Ysim)
    Si <- tryCatch(chol2inv(chol(S)), error = function(e)
      stop("Simulated covariance is not positive-definite: the sample is ",
           "degenerate (K too small, or a time point with no variability)."))
    # Large distance = low density, so the inequality flips relative to the
    # density version. Under elliptical contours the two rankings coincide.
    d2_sim <- stats::mahalanobis(Ysim,  mu, Si, inverted = TRUE)
    d2_obs <- stats::mahalanobis(y_obs, mu, Si, inverted = TRUE)
    med    <- mean(d2_sim >= as.numeric(d2_obs))
    stat   <- as.numeric(d2_obs)

  } else {
    if (!requireNamespace("ks", quietly = TRUE)) {
      stop('method = "kde" needs the ks package: install.packages("ks") or ',
           'renv::install("ks"). Use method = "mahalanobis" meanwhile.')
    }
    if (d > 4L) {
      stop("Binned kernel density estimation is limited to 4 dimensions and ",
           "this patient has ", d, " observations. Use method = ",
           '"mahalanobis".')
    }
    # Normal-scale bandwidth rather than the plug-in Hpi: with K around 1e5,
    # Hpi spends minutes searching for an optimal smoothing parameter -- the
    # practical bottleneck the paper itself flags in its discussion -- while in
    # log space the target is close enough to Gaussian for the normal-scale
    # rule. binned = TRUE is what keeps evaluating f at all K simulated points
    # from being an O(K^2) operation: bin onto a grid, convolve by FFT,
    # interpolate back.
    if (d == 1L) {
      x  <- as.numeric(Ysim)
      fh <- ks::kde(x = x, h = ks::hns(x), binned = TRUE, bgridsize = bgridsize)
      dens_sim <- stats::predict(fh, x = x)
      dens_obs <- stats::predict(fh, x = as.numeric(y_obs))
    } else {
      fh <- ks::kde(x = Ysim, H = ks::Hns(Ysim), binned = TRUE,
                    bgridsize = rep(bgridsize, d))
      dens_sim <- stats::predict(fh, x = Ysim)
      dens_obs <- stats::predict(fh, x = y_obs)
    }
    med  <- mean(dens_sim <= as.numeric(dens_obs))
    stat <- as.numeric(dens_obs)
  }

  list(med = med, statistic = stat, K = nrow(Ysim), n_obs = d,
       method = method, scale = scale, frac_censored = frac_censored)
}


# --- whole cohort, one model -------------------------------------------------
med_cohort <- function(data, prior_model, K = 1e5,
                       method = c("kde", "mahalanobis"),
                       scale  = c("log", "identity"),
                       allow_iov = FALSE, floor = 0.01,
                       nonneg = c("truncate", "none"), verbose = TRUE) {

  method <- match.arg(method)
  scale  <- match.arg(scale)
  nonneg <- match.arg(nonneg)
  ids    <- sort(unique(data$ID))

  out <- lapply(ids, function(i) {
    if (verbose) message("  ID ", i)
    tryCatch({
      Y <- med_simulate_patient(data[data$ID == i, , drop = FALSE],
                                prior_model, K = K, allow_iov = allow_iov,
                                nonneg = nonneg)
      r <- med_from_sample(Y, method = method, scale = scale, floor = floor)
      data.frame(ID = i, n_obs = r$n_obs, med = r$med,
                 frac_censored = r$frac_censored, failed = FALSE)
    }, error = function(e) {
      if (verbose) message("    failed: ", conditionMessage(e))
      data.frame(ID = i, n_obs = NA_integer_, med = NA_real_,
                 frac_censored = NA_real_, failed = TRUE)
    })
  })

  res <- do.call(rbind, out)
  attr(res, "K") <- K
  attr(res, "method") <- method
  attr(res, "scale")  <- scale
  res
}


# --- the two decision rules --------------------------------------------------
# GLOBAL: are the MEDs U(0,1)? A significant KS test invalidates the model on
# this cohort (paper Sect. 2.3). The direction carries the diagnosis: MEDs
# bunched near 0 say the observed profiles are more extreme than the model can
# routinely produce (too narrow, or the wrong structure); bunched near 1 says
# the model is over-dispersed relative to the data.
med_uniformity <- function(med, ...) {
  med <- med[is.finite(med)]
  kt  <- suppressWarnings(stats::ks.test(med, "punif", 0, 1, ...))
  list(n = length(med), D = unname(kt$statistic), p.value = kt$p.value,
       mean = mean(med), test = kt)
}

# PER PATIENT: the threshold is NOT alpha. With N patients the smallest of N
# uniforms is small by construction, so the paper controls the probability that
# the MINIMUM falls below the bound (Eqs. 5a/5b):
#     LL_alpha = 1 - (1 - alpha)^(1/N)
# For N = 39 and alpha = 0.05 that is 0.0013 -- a patient needs a MED below
# 0.13% to be called an outlier. That number also sets the resolution K has to
# deliver: K must be large enough that about K * LL_alpha simulated profiles
# can sit below the observed one at all, so K >= 1e5 here, not 1e3.
med_outlier_threshold <- function(N, alpha = 0.05) 1 - (1 - alpha)^(1 / N)

med_flag_outliers <- function(med_tbl, alpha = 0.05) {
  N  <- sum(!is.na(med_tbl$med))
  LL <- med_outlier_threshold(N, alpha)

  # TWO WAYS THE OUTLIER TEST CAN BE UNABLE TO FIRE, both silent. Checked here
  # because "0 outliers" otherwise reads as a clean bill of health for a model
  # that was never actually given the chance to fail.
  K <- attr(med_tbl, "K")
  if (!is.null(K) && 1 / K > LL) {
    warning(sprintf(
      "K = %s gives a MED resolution of %.2g, coarser than LL_alpha = %.2g: ",
      format(K, scientific = FALSE), 1 / K, LL),
      "no patient can be flagged. Raise K to at least ", ceiling(10 / LL), ".")
  }
  fc   <- max(med_tbl$frac_censored, na.rm = TRUE)
  nmax <- suppressWarnings(max(med_tbl$n_obs, na.rm = TRUE))
  if (is.finite(fc) && is.finite(nmax) && fc * nmax > LL) {
    n_heavy <- sum(med_tbl$frac_censored > 0.01, na.rm = TRUE)
    warning(sprintf(
      "Censoring at the log floor reaches %.2g of simulated values for the ",
      fc), sprintf("worst of %d patient(s) above 1%%. ", n_heavy),
      "With nonneg = \"truncate\" this is not a simulation artefact: it means ",
      "the model predicts sub-floor concentrations where drug was measured. ",
      "Those patients are still flagged correctly, but their MED is a LOWER ",
      "BOUND on how atypical they are, because profiles collapsed onto the ",
      "floor share a coordinate and inflate the simulated spread. Check ",
      "frac_censored per patient before quoting individual MED values.")
  }

  med_tbl$LL_alpha <- LL
  med_tbl$outlier  <- !is.na(med_tbl$med) & med_tbl$med <= LL
  attr(med_tbl, "alpha") <- alpha
  med_tbl
}
