#' Moment-based measurement error correction for the concentration index
#' 
#' This function returns the moment-based concentration index estimates, correcting for measurement error in the ranking variable by incorporating partial validation data.
#'
#' @param outcome vector containing the outcome outcomes for all observations
#' @param unval_exposure vector (of the same length as `outcome`) containing the error-prone, unvalidated exposure on which all observations will be ranked
#' @param val_exposure vector (of the same length as `outcome`) containing the error-free, validated exposure on which all observations will be ranked. For observations that were not validated, this vector should contain `NA`.
#' @param return_naive logical for whether the naive estimate (based only on `unval_exposure`) should be returned (if `TRUE`). The default is `FALSE`.
#' @param return_cc logical for whether the complete-case estimate (based only on those with non-missing `val_exposure`) should be returned (if `TRUE`). The default is `FALSE`.
#' @param include_se logical for whether standard error estimate should be returned. The default is `FALSE`.
#' @param bootstraps scalar for how many bootstrap resamples should be used to estimate standard errors (if `include_se = TRUE`). For jackknife (leave-one-out) resampling, `bootstraps = 0` (the default).
#' @param conf_level scalar (between 0 and 1) for confidence level of confidence intervals (if `include_se = TRUE`). The default is `0.95` for 95% confidence intervals. 
#' @param rank_ascend logical for whether exposure is ranked in ascending or descending order. The default is `TRUE`.
#' @return a list with the concentration index estimates requested (and standard errors, where applicable)
#' @export
#' 
mochi <- function(outcome, unval_exposure, val_exposure, return_naive = FALSE, return_cc = FALSE, include_se = FALSE, bootstraps = 0, conf_level = 0.95, rank_ascend = TRUE) {
  # If requested, negate exposures to make larger values more disadvantaged (lower rank)
  if(!rank_ascend) {
    unval_exposure <- unval_exposure * -1
    val_exposure <- val_exposure * -1
  }
  
  # Estimate concentration index 
  c_hat <- mochi_estimate(outcome = outcome, 
                          unval_exposure = unval_exposure, 
                          val_exposure = val_exposure, 
                          return_naive = return_naive, 
                          return_cc = return_cc, 
                          conf_level = conf_level) 
  
  # If requested, calculate jacknife standard error
  if(include_se) {
    if (bootstraps == 0) {
      n <- length(outcome) ## total number of observations
      jack_ci_mb <- rep(NA, n) ## initialize empty vector
      for(i in 1:n) { ## loop over sample, leaving out each row one-by-one
        ## Define vectors excluding row i 
        Y <- outcome[-i] 
        Xstar <- unval_exposure[-i]
        X <- val_exposure[-i]
        
        ## Compute moment-based concentration index excluding row i
        jack_ci_mb[i] <- mochi_estimate(
          outcome = Y, 
          unval_exposure = Xstar,
          val_exposure = X,
          return_naive = FALSE)[[1]]
      }
      ## Calculate standard error from vector using jacknife formula
      se <- sqrt((n - 1) / n * sum((jack_ci_mb - mean(jack_ci_mb)) ^ 2))
      lb <- quantile(x = jack_ci_mb, probs = (1 - conf_level) / 2)
      ub <- quantile(x = jack_ci_mb, probs = (1 - (1 - conf_level) / 2))
    } else {
      n <- length(outcome) ## total number of observations
      nv <- sum(!is.na(val_exposure)) ## number of validated observations
      which_nv <- which(!is.na(val_exposure)) ## indices of validated observations
      boot_ci_mb <- rep(NA, bootstraps) ## initialize empty vector
      for(i in 1:bootstraps) { ## loop over sample, resampling n each time
        ## Resample indices (stratified by validation)
        id_nv <- sample(x = which_nv, size = nv, replace = TRUE)
        id <- sample(x = setdiff(1:n, which_nv), size = (n - nv), replace = TRUE)
        ## Define vectors excluding row i 
        Y <- outcome[c(id_nv, id)] 
        Xstar <- unval_exposure[c(id_nv, id)]
        X <- val_exposure[c(id_nv, id)]
        
        ## Compute moment-based concentration index excluding row i
        boot_ci_mb[i] <- mochi_estimate(
          outcome = Y, 
          unval_exposure = Xstar,
          val_exposure = X,
          return_naive = FALSE)[[1]]
      }
      ## Calculate standard error and central conf_level x 100% interval from vector
      se <- sd(boot_ci_mb)
      lb <- quantile(x = boot_ci_mb, probs = (1 - conf_level) / 2)
      ub <- quantile(x = boot_ci_mb, probs = (1 - (1 - conf_level) / 2))
    }
  } else {
    se <- lb <- ub <- NA
  }
  
  # Return estimates
  if (return_naive) {
    ## Calculate Wald-type confidence intervals for naive 
    lb_ci_naive <- c_hat[["ci_naive"]] - 1.96 * c_hat[["se_ci_naive"]]
    ub_ci_naive <- c_hat[["ci_naive"]] + 1.96 * c_hat[["se_ci_naive"]]
    if (return_cc) {
      return(
        list(
          ci_moment = c_hat[["ci_moment"]],
          se_ci_moment = se,
          confint_ci_moment = c(lb, ub), 
          ci_naive = c_hat[["ci_naive"]], 
          se_ci_naive = c_hat[["se_ci_naive"]],
          confint_ci_naive = c(lb_ci_naive, ub_ci_naive),
          ci_cc = c_hat[["ci_cc"]], 
          se_ci_cc = c_hat[["se_ci_cc"]],
          confint_ci_cc = c(c_hat[["lb_ci_cc"]], c_hat[["ub_ci_cc"]])
        )
      )
    } else if (return_cc) {
      return(
        list(
          ci_moment = c_hat[["ci_moment"]],
          se_ci_moment = se,
          confint_ci_moment = c(lb, ub), 
          ci_cc = c_hat[["ci_cc"]], 
          se_ci_cc = c_hat[["se_ci_cc"]],
          confint_ci_cc = c(c_hat[["lb_ci_cc"]], c_hat[["ub_ci_cc"]])
        )
      )
    } else {
      return(
        list(
          ci_moment = c_hat[["ci_moment"]],
          se_ci_moment = se,
          confint_ci_moment = c(lb, ub), 
          ci_naive = c_hat[["ci_naive"]], 
          se_ci_naive = c_hat[["se_ci_naive"]],
          confint_ci_naive = c(lb_ci_naive, ub_ci_naive)
        )
      )
    }
  } else {
    return(
      list(
        ci_moment = c_hat[["ci_moment"]],
        se_ci_moment = se, 
        confint_ci_moment = c(lb, ub)
        )
    )
  }
}
