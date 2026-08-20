# ============================================================================
# messieR_core.R
#
# The reconciled canonical core for the thesis and the messieR package.
#
# This file resolves the synergy problem: M2 defined
#   sim_var1_realistic(vars, N, delta, psi, residual_pool, sigma_eps, ...)
# and the earlier bootstrap core defined
#   sim_var1(T_total, phi, Sigma, innov_fn, ...)
# These were two contracts for one job. Below they are folded into a single
# sim_var1() whose error source is chosen by an explicit precedence chain:
#
#   residual_pool  (empirical resampling, iid or block)   -- M2's path
#     else innov_fn  (caller-supplied generator)          -- bootstrap core path
#       else sigma   (parametric MVN)                      -- M2's fallback
#         else error
#
# Everything downstream -- the thesis scripts and the messieR package alike --
# calls these functions and no others. M2 and M3 are retained in archive/ for
# provenance but are superseded by this file.
#
# Register: this is package code, so it is written defensively (stopifnot
# validation, tryCatch around fragile fits, roxygen on the exported API). The
# thesis ANALYSIS SCRIPTS that source this file are written in the leaner
# procedural style; the two registers are deliberate.
#
# Dependencies: MASS (mvrnorm). sn is needed only by callers that build a
# skew-normal innov_fn; it is not required to load this file.
#
# Author: Karthik Tumuluru
# ============================================================================

# Note for packaging: when this becomes messieR/, split each titled section
# into its own R/ file (sim_var1.R, fit_var1.R, boot.R, missingness.R, mcse.R)
# and move the `library(MASS)` call to Imports in DESCRIPTION with @importFrom.
if (requireNamespace("MASS", quietly = TRUE)) {
  # mvrnorm pulled in lazily by sim_var1 when a parametric source is used.
}


# ============================================================================
# Section 1. Stationarity
# ============================================================================

#' Check covariance-stationarity of a VAR(1) transition matrix
#'
#' @param phi A square numeric transition matrix.
#' @param tol Numerical tolerance on the spectral radius.
#' @return Logical; TRUE if all eigenvalues lie strictly inside the unit disc.
#' @export
is_stationary <- function(phi, tol = 1e-8) {
  stopifnot(is.matrix(phi), nrow(phi) == ncol(phi))
  max(Mod(eigen(phi, only.values = TRUE)$values)) < 1 - tol
}


# ============================================================================
# Section 2. Canonical VAR(1) simulator (the reconciliation)
# ============================================================================

#' Simulate a stationary VAR(1) series with a choice of error source
#'
#' One simulator, three error sources, chosen by precedence:
#' \enumerate{
#'   \item \code{residual_pool} non-NULL: innovations are resampled from the
#'     pool (rows are k-dimensional residual vectors, so the contemporaneous
#'     correlation structure is preserved). \code{resample = "iid"} draws rows
#'     with replacement; \code{resample = "block"} draws contiguous blocks.
#'   \item else \code{innov_fn} non-NULL: innovations are \code{innov_fn(n)},
#'     which must return an n-by-k matrix.
#'   \item else \code{sigma} non-NULL: innovations are MVN(0, sigma).
#'   \item else: error.
#' }
#'
#' @param n_obs Number of observations to return (after burn-in).
#' @param phi The k-by-k transition matrix (must be stationary).
#' @param delta Length-k intercept vector; defaults to zeros.
#' @param sigma k-by-k innovation covariance for the parametric source.
#' @param residual_pool An m-by-k matrix of residual vectors to resample.
#' @param innov_fn A function(n) returning an n-by-k innovation matrix.
#' @param resample Resampling scheme for the pool source: "iid" or "block".
#' @param block_size Block length for "block" resampling; defaults to the
#'   Jentsch-Lunsford rule ceiling(5.03 * (n_obs + burn_in)^(1/4)).
#' @param burn_in Number of initial observations discarded.
#' @return An n_obs-by-k numeric matrix.
#' @export
sim_var1 <- function(n_obs, phi,
                     delta = NULL, sigma = NULL,
                     residual_pool = NULL, innov_fn = NULL,
                     resample = c("iid", "block"), block_size = NULL,
                     burn_in = 1000L) {
  resample <- match.arg(resample)
  stopifnot(is.matrix(phi), nrow(phi) == ncol(phi), is_stationary(phi))
  k <- nrow(phi)
  if (is.null(delta)) delta <- rep(0, k)
  stopifnot(length(delta) == k)

  N <- n_obs + burn_in

  # --- error source precedence: pool > innov_fn > sigma -------------------
  if (!is.null(residual_pool)) {
    pool <- as.matrix(residual_pool)
    stopifnot(ncol(pool) == k)
    E <- resample_pool(pool, N, method = resample, block_size = block_size)
  } else if (!is.null(innov_fn)) {
    E <- innov_fn(N)
    stopifnot(is.matrix(E), nrow(E) == N, ncol(E) == k)
  } else if (!is.null(sigma)) {
    stopifnot(is.matrix(sigma), nrow(sigma) == k)
    E <- MASS::mvrnorm(N, mu = rep(0, k), Sigma = sigma)
  } else {
    stop("sim_var1: supply one of residual_pool, innov_fn, or sigma.")
  }

  # --- recursive construction with burn-in --------------------------------
  Y <- matrix(0, N, k)
  for (t in 2:N) {
    Y[t, ] <- delta + phi %*% Y[t - 1L, ] + E[t, ]
  }
  Y[(burn_in + 1L):N, , drop = FALSE]
}


#' Resample residual rows, iid or in fixed blocks
#'
#' Helper used by both the simulator's pool source and the bootstraps.
#' Recycles rows as a last resort so the return always has T_target rows.
#'
#' @param pool An m-by-k residual matrix.
#' @param T_target Number of rows to return.
#' @param method "iid" or "block".
#' @param block_size Block length for "block"; defaults to the J-L rule.
#' @return A T_target-by-k matrix.
#' @export
resample_pool <- function(pool, T_target,
                          method = c("iid", "block"), block_size = NULL) {
  method <- match.arg(method)
  m <- nrow(pool)
  if (method == "iid") {
    return(pool[sample.int(m, T_target, replace = TRUE), , drop = FALSE])
  }
  if (is.null(block_size)) block_size <- default_block_length(T_target)
  if (block_size > m) {
    warning("resample_pool: block_size > pool size; using iid.")
    return(pool[sample.int(m, T_target, replace = TRUE), , drop = FALSE])
  }
  n_blocks <- ceiling(T_target / block_size)
  starts <- sample.int(m - block_size + 1L, n_blocks, replace = TRUE)
  out <- do.call(rbind, lapply(starts, function(s)
    pool[s:(s + block_size - 1L), , drop = FALSE]))
  if (nrow(out) < T_target)
    out <- out[rep(seq_len(nrow(out)), length.out = T_target), , drop = FALSE]
  out[1:T_target, , drop = FALSE]
}


#' Jentsch-Lunsford (2019) block length: ceiling(5.03 * T^(1/4))
#' @param T_use Series length.
#' @export
default_block_length <- function(T_use) {
  max(1L, as.integer(ceiling(5.03 * T_use^(1 / 4))))
}

#' Politis-White (2004)-style stationary-bootstrap probability: 1 / ceil(2.5 T^.25)
#' @param T_use Series length.
#' @export
default_stationary_p <- function(T_use) {
  1 / max(2L, as.integer(ceiling(2.5 * T_use^(1 / 4))))
}


# ============================================================================
# Section 3. OLS estimation (transition-level complete cases)
# ============================================================================

#' Fit a bivariate (or k-variate) VAR(1) by equation-by-equation OLS
#'
#' Uses only transitions whose two endpoints are both observed and
#' adjacent in time: the series is lagged first and incomplete pairs are
#' dropped second, so incomplete series need no separate code path. Returns NULL on samples
#' too small or singular to fit, which the bootstraps treat as a failed draw.
#'
#' @param Y A T-by-k series, possibly with NA rows.
#' @param min_obs Minimum complete transitions required to attempt a fit.
#' @return A list (phi_hat, intercept, Sigma_hat, residuals, se, pvals) or NULL.
#' @export
fit_var1 <- function(Y, min_obs = 10L, full = TRUE) {
  Y <- as.matrix(Y)
  T_use <- nrow(Y)
  if (T_use < min_obs) return(NULL)
  k <- ncol(Y)

  # GAP-AWARE PAIRING. Build the transition pairs from the ORIGINAL series and
  # keep only those where both endpoints are observed AND adjacent in time.
  # Collapsing with complete.cases() first would splice non-adjacent rows into
  # false one-step transitions and attenuate phi badly (see Section 3.5).
  Y_t   <- Y[2:T_use, , drop = FALSE]
  Y_tm1 <- Y[1:(T_use - 1L), , drop = FALSE]
  ok <- stats::complete.cases(Y_t) & stats::complete.cases(Y_tm1)
  n_ok <- sum(ok)
  if (n_ok < min_obs) return(NULL)
  Y_t   <- Y_t[ok, , drop = FALSE]
  Y_tm1 <- Y_tm1[ok, , drop = FALSE]

  # Both equations regress on the same design [1, Y_tm1], so ONE cross-product
  # inverse serves both. Replaces two lm() + two summary.lm() calls; this is
  # the change that took the production grid from ~4 days to ~5 hours.
  X <- cbind(1, Y_tm1)
  XtX_inv <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(XtX_inv)) return(NULL)

  B_hat  <- XtX_inv %*% crossprod(X, Y_t)      # (k+1) x k, column j = eqn j
  resid  <- Y_t - X %*% B_hat
  df_res <- n_ok - (k + 1L)
  if (df_res < 1L) return(NULL)

  sigma2  <- colSums(resid^2) / df_res
  se_full <- sqrt(outer(diag(XtX_inv), sigma2))
  phi_hat <- t(B_hat[-1, , drop = FALSE])      # row = equation i, col = pred j
  se      <- t(se_full[-1, , drop = FALSE])

  # bootstrap refits need only these two; skipping the rest is a large saving
  if (!full) return(list(phi_hat = phi_hat, se = se))

  tstat <- phi_hat / se
  list(phi_hat   = phi_hat,
       intercept = B_hat[1, ],
       Sigma_hat = crossprod(resid) / df_res,
       residuals = resid,
       se        = se,
       pvals     = 2 * stats::pt(-abs(tstat), df = df_res))
}


# ============================================================================
# Section 4. Monte Carlo standard errors (Siepe et al. 2024)
# ============================================================================

#' @rdname mcse
#' @param x Bernoulli-coded vector (rejections, coverage indicators).
#' @param n_sim Replicate count.
#' @export
mcse_proportion <- function(x, n_sim = sum(!is.na(x))) {
  p <- mean(x, na.rm = TRUE); sqrt(p * (1 - p) / n_sim)
}
#' @rdname mcse
#' @param estimates Numeric vector of per-replicate estimates.
#' @export
mcse_mean <- function(estimates, n_sim = sum(!is.na(estimates)))
  stats::sd(estimates, na.rm = TRUE) / sqrt(n_sim)
#' @rdname mcse
#' @param true_value Scalar true value (unused beyond signature symmetry).
#' @export
mcse_bias <- function(estimates, true_value = NULL,
                      n_sim = sum(!is.na(estimates)))
  stats::sd(estimates, na.rm = TRUE) / sqrt(n_sim)
#' Monte Carlo standard errors for simulation summaries
#' @name mcse
#' @export
mcse_mse <- function(estimates, true_value, n_sim = sum(!is.na(estimates)))
  stats::sd((estimates - true_value)^2, na.rm = TRUE) / sqrt(n_sim)


# ============================================================================
# Section 5. Missingness injectors (harmonised from M2)
# ============================================================================

#' Inject MCAR missingness into one column
#' @param Y A T-by-k series.
#' @param rate Marginal missingness rate.
#' @param col Column to thin (default 2).
#' @return Y with selected entries of `col` set to NA.
#' @export
inject_mcar <- function(Y, rate = 0.30, col = 2L) {
  Y <- as.matrix(Y)
  Y[stats::runif(nrow(Y)) < rate, col] <- NA
  Y
}

#' Inject MAR missingness driven by another column's prior value
#'
#' P(Y[t, col] missing | Y[t-1, predictor]) = logit^{-1}(alpha + coef * z),
#' with alpha solved per call to hit the target marginal rate.
#'
#' @param Y A T-by-k series.
#' @param rate Target marginal missingness rate.
#' @param coef Logistic coefficient on the standardised predictor.
#' @param col Column to thin (default 2).
#' @param predictor Column whose lag drives missingness (default 1).
#' @return Y with selected entries of `col` set to NA.
#' @export
inject_mar <- function(Y, rate = 0.30, coef = 0.5, col = 2L, predictor = 1L) {
  Y <- as.matrix(Y); T_use <- nrow(Y)
  lag_pred <- c(0, Y[1:(T_use - 1L), predictor])
  z <- as.numeric(scale(lag_pred))
  z[is.na(z)] <- 0
  f <- function(a) mean(stats::plogis(a + coef * z)) - rate
  alpha <- tryCatch(stats::uniroot(f, c(-10, 10))$root,
                    error = function(e) stats::qlogis(rate))
  drop <- stats::runif(T_use) < stats::plogis(alpha + coef * z)
  Y[drop, col] <- NA
  Y
}


# ============================================================================
# Section 6. Bootstrap inference
# ============================================================================
# Two primary procedures (model-based residual; stationary) and a moving-block
# sensitivity variant. All three share the summary machinery in boot_summary().

#' Stationary-bootstrap resample of a series (Politis & Romano 1994)
#' @param Y A T-by-k series.
#' @param p Per-step block-end probability (geometric block lengths).
#' @return A T-by-k resampled series with circular wrap.
#' @export
stationary_resample <- function(Y, p = NULL) {
  Y <- as.matrix(Y); T_use <- nrow(Y)
  if (is.null(p)) p <- default_stationary_p(T_use)
  # Vectorised: position 1 starts a block, each later position starts a new one
  # with probability p, and within a block the index walks forward circularly
  # from a random start. Same stationary bootstrap, drawn without a loop.
  is_new       <- c(TRUE, stats::runif(T_use - 1L) < p)
  block_id     <- cumsum(is_new)
  starts_all   <- sample.int(T_use, T_use, replace = TRUE)
  start_of_blk <- starts_all[is_new][block_id]
  first_pos    <- which(is_new)[block_id]
  offset       <- seq_len(T_use) - first_pos
  idx          <- ((start_of_blk - 1L + offset) %% T_use) + 1L
  Y[idx, , drop = FALSE]
}

# internal: turn a B-by-4 bootstrap estimate matrix into the summary list.
boot_summary <- function(boot_phi, boot_t, phi, se_orig, method, alpha,
                         block_param) {
  pv <- as.vector(t(phi)); sv <- as.vector(t(se_orig))
  ci_pct <- apply(boot_phi, 2, stats::quantile,
                  probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  t_lo <- apply(boot_t, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  t_hi <- apply(boot_t, 2, stats::quantile, probs = alpha / 2,     na.rm = TRUE)
  pbc  <- 2 * pv - colMeans(boot_phi, na.rm = TRUE)
  ppv  <- apply(boot_phi, 2, function(v) {
    v <- v[!is.na(v)]; if (!length(v)) return(NA_real_)
    2 * min(mean(v <= 0), mean(v >= 0), 0.5)
  })
  m <- function(v) matrix(v, nrow = 2, byrow = TRUE)
  list(method = method, phi_hat = phi, phi_bc = m(pbc),
       ci_percentile = list(lo = m(ci_pct[1, ]), hi = m(ci_pct[2, ])),
       ci_student_t  = list(lo = m(pv - t_lo * sv), hi = m(pv - t_hi * sv)),
       pvals_boot = m(ppv), n_boot_failed = sum(is.na(boot_phi[, 1])),
       block_param = block_param)
}

# internal: run B refits given a resampler that returns a bootstrap series.
boot_engine <- function(Y, phi, se_orig, B, alpha, make_series, method,
                        block_param) {
  pv <- as.vector(t(phi))
  bp <- matrix(NA_real_, B, 4L); bt <- matrix(NA_real_, B, 4L)
  for (b in seq_len(B)) {
    Ys <- make_series()
    fs <- tryCatch(fit_var1(Ys, full = FALSE), error = function(e) NULL)
    if (!is.null(fs)) {
      vs <- as.vector(t(fs$phi_hat)); ss <- as.vector(t(fs$se))
      bp[b, ] <- vs; bt[b, ] <- (vs - pv) / ss
    }
  }
  boot_summary(bp, bt, phi, se_orig, method, alpha, block_param)
}

#' Model-based i.i.d. residual bootstrap for VAR(1) inference
#' @param Y A complete T-by-k series.
#' @param B Bootstrap replicates.
#' @param alpha 1 - confidence level.
#' @return A bootstrap summary list.
#' @export
boot_residual <- function(Y, B = 999L, alpha = 0.05) {
  fit <- fit_var1(Y); if (is.null(fit)) return(NULL)
  phi <- fit$phi_hat; delta <- fit$intercept
  resid <- fit$residuals
  Ymat <- as.matrix(Y); T_use <- nrow(Ymat)
  na_mask <- is.na(Ymat)
  # Start from the first FULLY OBSERVED row so the recursion never ingests an
  # NA. Under missingness, row 1 has a ~30% chance of being incomplete, which
  # propagated NA through the whole bootstrap series and silently drove
  # bootstrap power down in exactly the missingness cells. With complete data
  # this is row 1, so the procedure is unchanged there.
  start <- Ymat[which(stats::complete.cases(Ymat))[1L], ]
  make <- function() {
    e <- resample_pool(resid, T_use - 1L, method = "iid")
    Ys <- matrix(0, T_use, ncol(phi)); Ys[1L, ] <- start
    for (t in 2:T_use) Ys[t, ] <- delta + phi %*% Ys[t - 1L, ] + e[t - 1L, ]
    # Reimpose the observed missingness pattern, so every bootstrap refit faces
    # the same gap structure the data did.
    Ys[na_mask] <- NA
    Ys
  }
  boot_engine(Y, phi, fit$se, B, alpha, make, "residual_iid", NA)
}

#' Stationary bootstrap for VAR(1) inference (Politis & Romano 1994)
#' @inheritParams boot_residual
#' @param p Block-end probability; defaults to the Politis-White rule.
#' @export
boot_stationary <- function(Y, B = 999L, alpha = 0.05, p = NULL) {
  fit <- fit_var1(Y); if (is.null(fit)) return(NULL)
  if (is.null(p)) p <- default_stationary_p(nrow(as.matrix(Y)))
  make <- function() stationary_resample(Y, p = p)
  boot_engine(Y, fit$phi_hat, fit$se, B, alpha, make, "stationary", p)
}

#' Moving-block bootstrap for VAR(1) inference (sensitivity variant)
#' @inheritParams boot_residual
#' @param block_size Fixed block length; defaults to the J-L rule.
#' @export
boot_moving_block <- function(Y, B = 999L, alpha = 0.05, block_size = NULL) {
  fit <- fit_var1(Y); if (is.null(fit)) return(NULL)
  Ym <- as.matrix(Y); T_use <- nrow(Ym)
  if (is.null(block_size)) block_size <- default_block_length(T_use)
  make <- function() resample_pool(Ym, T_use, method = "block",
                                   block_size = block_size)
  boot_engine(Y, fit$phi_hat, fit$se, B, alpha, make, "moving_block",
              block_size)
}


# ============================================================================
# Section 7. Unified per-replicate inference (tidy output for the grid driver)
# ============================================================================

#' Run Wald + both bootstraps on one series; return one tidy row per coef
#'
#' @param Y A T-by-k observed series. Rows containing NA are permitted;
#'   they are resampled like any other row, and the downstream refit's
#'   transition pairing handles the resulting gaps.
#' @param B Bootstrap replicates.
#' @param alpha 1 - confidence level.
#' @param true_phi The k-by-k true matrix, for coverage/bias columns.
#' @param include_moving_block Also run the moving-block sensitivity variant.
#' @return A data.frame: one row per (procedure x coefficient).
#' @export
boot_var1_all <- function(Y, B = 999L, alpha = 0.05, true_phi = NULL,
                          include_moving_block = FALSE) {
  fit <- fit_var1(Y); if (is.null(fit)) return(NULL)
  zc <- stats::qnorm(1 - alpha / 2)

  procs <- list(
    wald = list(est = fit$phi_hat, bc = matrix(NA_real_, 2, 2), se = fit$se,
                lo = fit$phi_hat - zc * fit$se, hi = fit$phi_hat + zc * fit$se,
                p = fit$pvals))
  br <- boot_residual(Y, B, alpha)
  procs[["residual_iid"]] <- list(est = br$phi_hat, bc = br$phi_bc,
    se = matrix(NA_real_, 2, 2), lo = br$ci_student_t$lo,
    hi = br$ci_student_t$hi, p = br$pvals_boot)
  bs <- boot_stationary(Y, B, alpha)
  procs[["stationary"]] <- list(est = bs$phi_hat, bc = bs$phi_bc,
    se = matrix(NA_real_, 2, 2), lo = bs$ci_student_t$lo,
    hi = bs$ci_student_t$hi, p = bs$pvals_boot)
  if (include_moving_block) {
    bm <- boot_moving_block(Y, B, alpha)
    procs[["moving_block"]] <- list(est = bm$phi_hat, bc = bm$phi_bc,
      se = matrix(NA_real_, 2, 2), lo = bm$ci_student_t$lo,
      hi = bm$ci_student_t$hi, p = bm$pvals_boot)
  }

  cn <- c("phi11", "phi12", "phi21", "phi22")
  rows <- list()
  for (pn in names(procs)) {
    P <- procs[[pn]]
    for (kk in seq_along(cn)) {
      i <- ((kk - 1L) %/% 2L) + 1L; j <- ((kk - 1L) %% 2L) + 1L
      truth <- if (!is.null(true_phi)) true_phi[i, j] else NA_real_
      lo <- P$lo[i, j]; hi <- P$hi[i, j]; pval <- P$p[i, j]; est <- P$est[i, j]
      ebc <- P$bc[i, j]
      rows[[length(rows) + 1L]] <- data.frame(
        procedure = pn, coef = cn[kk], est = est, est_bc = ebc,
        se = P$se[i, j], ci_lo = lo, ci_hi = hi, pval = pval,
        reject_05 = as.integer(!is.na(pval) && pval < 0.05),
        covered = if (!is.na(truth)) as.integer(lo <= truth & truth <= hi) else NA_integer_,
        bias = if (!is.na(truth)) est - truth else NA_real_,
        bc_bias = if (!is.na(truth) && !is.na(ebc)) ebc - truth else NA_real_,
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}


# ============================================================================
# Section 8. Smoke test (uncomment to run)
# ============================================================================
# phi <- matrix(c(0.5, 0.1, 0.1, 0.5), 2, byrow = TRUE)
#
# # all three error sources through the one simulator:
# Y_par  <- sim_var1(100, phi, sigma = diag(2))                      # parametric
# Y_fn   <- sim_var1(100, phi, innov_fn = function(n) matrix(rt(2*n, 4), n)) # custom
# pool   <- fit_var1(Y_par)$residuals
# Y_pool <- sim_var1(100, phi, residual_pool = pool, resample = "block")     # empirical
#
# # missingness + unified inference:
# Ym <- inject_mar(Y_par, rate = 0.30)
# print(boot_var1_all(Ym, B = 199, true_phi = phi, include_moving_block = TRUE))
