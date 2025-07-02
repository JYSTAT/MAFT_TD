
required_packages <- c("nloptr", "survival", "dplyr", "nnls")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

#' Main MAFT Estimation Algorithm for Time-dependent Covariates
#'
#' @param dat Dataframe with columns: ID, start, stop, delta, covariates
#' @param X Character vector of covariate names
#' @param process Logical; if TRUE, iteration progress is printed
#' @param fig Logical; if TRUE, directional derivative plot is shown
#' @param D Threshold for convergence of directional derivative
#' @param CNM_iter Max number of CNM iterations
#' @param EM_itr Number of EM iterations per CNM step
#' @param w_criteria Threshold to remove components with small weights
#' @param algorithm Optimization algorithm (default: "NLOPT_LN_COBYLA")
#' @param tol Tolerance level for optimization
#' @param maxiter Max iterations for optimization
#' @return List with Final_Est, Support, W, History, D, iter_N
MAFT_TD <- function(dat, X, process = FALSE, fig = FALSE, D = 1e-5,
                    CNM_iter = 100, EM_itr = 3, w_criteria = 1e-10,
                    algorithm = "NLOPT_LN_COBYLA", tol = 1e-5, maxiter = 2000) {
  
  # ---- 1. Initial Values ---- #
  dat1 <- dat %>% group_by(ID) %>% arrange(stop) %>% slice_tail(n = 1) %>% ungroup()
  X.mat <- as.matrix(dat1[, X])
  PAFT <- survreg(Surv(dat1$stop, dat1$delta) ~ X.mat, dist = "lognormal")
  
  beta <- PAFT$coefficients
  gab <- PAFT$scale
  cons <- min(gab) / 5
  bb <- gab / 2.5
  weight <- 1
  beta_lik_value <- NA
  History <- NULL
  
  # ---- 2. Initial Weight Update ---- #
  final <- T0_cal(dat, beta, X) |> final_Inf()
  ini <- update_weight_univ2(dat, X, bb, weight, beta, cons, w_criteria)
  weight <- ini$weight
  bb <- ini$bb
  
  # ---- 3. CNM Algorithm Loop ---- #
  for (k in 1:CNM_iter) {
    
    # Step 1. Update support (sigma) given beta
    find_re <- find_sigma_univ2(dat, beta, X, fig, weight, bb, cons,
                                algorithm, tol, maxiter, k)
    new_sigma <- find_re$newx
    ttt <- find_re$ttt
    max_indicator <- find_re$max_indicator
    
    # Record history
    compo_N <- length(bb)
    Result_Est <- c(beta, compo_N, min(bb), max(bb), min(weight), max(weight),
                    max(ttt), -beta_lik_value, max_indicator)
    names(Result_Est) <- c(paste0("beta", 0:length(X)),
                           "compo_N", "minSigma", "maxSigma", "minW", "maxW",
                           "D(H;G)", "loglikelihood", "MaxSupport")
    History <- rbind(History, Result_Est)
    rownames(History) <- paste0("CNM_iter=", 1:nrow(History))
    
    if (process) cat("CNM_iter =", k, "\n")
    
    # Step 2. Check convergence
    if (max(ttt) < D && k > 1) break
    
    # Step 3. Add new support and update weight
    bb <- unique(c(bb, new_sigma))
    weight <- rep(1 / length(bb), length(bb))
    
    up <- update_weight_univ2(dat, X, bb, weight, beta, cons, w_criteria)
    bb <- up$bb
    weight <- up$weight
    
    # Step 4. EM Iterations: update beta and weight
    for (jj in 1:EM_itr) {
      up1 <- myfun_univ2(dat, X, bb, weight, beta, cons, algorithm, tol, maxiter)
      beta <- up1$upBeta
      up2 <- update_weight_univ2(dat, X, bb, weight, beta, cons, w_criteria)
      bb <- up2$bb
      weight <- up2$weight
    }
    
    beta_lik_value <- up1$likelihood
  }
  
  if (process) {
    cat("\n===== Finish MAFT algorithm =====\n",
        "* Number_itr =", nrow(History), "\n",
        "* Convergence =", ifelse(nrow(History) == CNM_iter, "Fail...", "Success!"), "\n")
  }
  
  Final_Est <- History[nrow(History), 1:(length(X) + 1)]
  iter_N <- data.frame(CNM_iter = nrow(History), CNM_tot = CNM_iter, EM_itr = EM_itr)
  
  return(list(Final_Est = Final_Est,
              Support = bb,
              W = weight,
              History = History,
              D = D,
              iter_N = iter_N))
}


#---------------------------- function part -------------------------------#

#------------------------------------#
### ----- Function1. T0_cal -----  ###
#------------------------------------#

T0_cal <- function(dat, beta, Z) {
  #' Compute baseline failure time quantities (psi, psi_d) for time-dependent AFT model
  #'
  #' @param dat  Long-format data with columns: ID, start, stop, delta, Z1, Z2, ...
  #' @param beta Coefficient vector (first element for X, remaining for time-dependent Z)
  #' @param Z    Character vector of column names for time-dependent covariates
  #'
  #' @return A data.frame with one row per subject containing:
  #'         - ID: subject identifier
  #'         - psi: cumulative risk integral (∫ exp(-Z(u)^T β) du)
  #'         - psi_d: risk at final time (exp(-Z(y)^T β))
  #'         - delta: event indicator
  
  ## -----------------------
  ## 1. Compute Capital Psi
  ## ∫ from start to stop of exp(-Z(u)^T β) du
  ## Approximated as (stop - start) * exp(-Z(u)^T β) per interval
  ## -----------------------
  
  # Extract time-dependent covariates matrix
  Z_mat <- as.matrix(dat[, Z, drop = FALSE])
  
  # Compute linear predictor: -Z(u)^T * beta (excluding intercept)
  psi_linear <- - Z_mat %*% beta[-1]
  
  # Compute ψ = (stop - start) * exp(-Z(u)^T β)
  dat$psi <- (dat$stop - dat$start) * exp(psi_linear)
  
  # Aggregate ψ by subject (ID)
  psi_df <- dat %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(psi = sum(psi), .groups = "drop")
  
  ## -----------------------
  ## 2. Compute Small Psi_d
  ## exp(-Z(y)^T β): risk at the last time point per subject
  ## -----------------------
  
  # Get last row per subject (y_i = last stop time)
  last_rows <- dat %>%
    dplyr::group_by(ID) %>%
    dplyr::slice_tail(n = 1) %>%
    dplyr::ungroup() %>%
    dplyr::select(ID, all_of(Z), delta) %>%
    as.data.frame()
  
  # Linear predictor at final time
  Z_last_mat <- as.matrix(last_rows[, Z, drop = FALSE])
  psi_d <- exp(- Z_last_mat %*% beta[-1])
  last_rows$psi_d <- as.vector(psi_d)
  
  ## -----------------------
  ## 3. Merge results: psi, psi_d, delta
  ## -----------------------
  
  final <- merge(psi_df, last_rows[, c("ID", "psi_d", "delta")], by = "ID", all.x = TRUE)
  final <- final[order(final$ID), ]
  
  return(final)
}

#' Adjust extreme values in psi and psi_d for numerical stability
#'
#' This function caps large values of \code{psi} and \code{psi_d} at \eqn{\pm 1e+300},
#' and replaces exact zeros with \eqn{1e-300}, to avoid numerical issues in likelihood computations.
#'
#' @param final A data.frame containing columns \code{psi} and \code{psi_d}.
#' @return A data.frame with the same structure, with extreme values adjusted.
#'
#' @examples
#' # final <- final_Inf(final)
#'
#' @export
final_Inf <- function(final) {
  # Cap very large or very small values in 'psi'
  final$psi[final$psi > 1e+300] <- 1e+300      # Upper bound
  final$psi[final$psi < -1e+300] <- -1e+300    # Lower bound
  final$psi[final$psi == 0] <- 1e-300          # Replace exact zero
  
  # Same adjustment for 'psi_d'
  final$psi_d[final$psi_d > 1e+300] <- 1e+300
  final$psi_d[final$psi_d < -1e+300] <- -1e+300
  final$psi_d[final$psi_d == 0] <- 1e-300
  
  return(final)
}


#---------------------------------------------------------------------------------#
### -- Function 2. Likelihood for log-Normal AFT model (with censoring) ### --  ###
#---------------------------------------------------------------------------------#

#' Compute individual likelihood contributions under a censored log-normal model
#'
#' This function calculates the likelihood for each subject under the assumption that
#' log(T0) follows a Gaussian distribution with censoring taken into account.
#'
#' @param final A data.frame from T0_cal() containing psi, psi_d, and delta.
#' @param beta0 Scalar mean parameter (intercept term) for log(T0).
#' @param se Scalar standard deviation for log(T0).
#'
#' @return A numeric vector of likelihood contributions per subject.
#'
#' @export
normal_cen <- function(final, beta0, se) {
  # Density for observed events
  den <- dnorm((log(final$psi) - beta0) / se) * (final$psi_d / (final$psi * se))
  
  # Likelihood function with right-censoring
  lik <- final$delta * den +
    (1 - final$delta) * (1 - pnorm((log(final$psi) - beta0) / se))
  
  # Prevent likelihood from becoming too small (for numerical stability)
  lik <- ifelse(lik <= 1.0e-20, 1.0e-20, lik)
  
  return(lik)
}


#-----------------------------------------------------------------#
### Function 3. Likelihood for 2-Component Gaussian Mixture AFT ###
#-----------------------------------------------------------------#

#' Compute likelihood under a 2-component Gaussian mixture AFT model
#'
#' This function computes the subject-wise likelihood assuming log(T0)
#' follows a 2-component Gaussian mixture distribution, incorporating censoring.
#'
#' @param final A data.frame from T0_cal() containing psi, psi_d, and delta
#' @param beta A coefficient vector (first element = intercept, shared across components)
#' @param se A numeric vector of standard deviations for each Gaussian component (length = 2)
#' @param weight A numeric vector of component weights (length = 2, should sum to 1)
#'
#' @return A numeric vector of mixture likelihood values for each subject
#'
#' @export
normal_mixture2 <- function(final, beta, se, weight) {
  # Initialize likelihood accumulator
  answer <- rep(0, nrow(final))
  
  # Add weighted likelihoods from each component
  for (i in 1:length(se)) {
    answer <- answer + weight[i] * normal_cen(final, beta[1], se[i])
  }
  
  # Apply lower bound to avoid numerical underflow
  answer1 <- ifelse(answer <= 1.0e-20, 1.0e-20, answer)
  
  return(answer1)
}


#-----------------------------------------------------------------------#
### Function 4. Directional Derivative for Candidate Sigma in Mixture ###
#-----------------------------------------------------------------------#

#' Compute the gradient function D(phi; G) for a candidate sigma in Gaussian mixture AFT
#'
#' This function evaluates the directional derivative (gradient function) when adding a candidate
#' Gaussian component with standard deviation \code{phi} to an existing mixture distribution G.
#'
#' @param final A data.frame from T0_cal() with psi, psi_d, and delta
#' @param phi Candidate standard deviation to test
#' @param beta Regression coefficient vector (intercept + others)
#' @param weight Numeric vector of weights for the current mixture components
#' @param se Numeric vector of standard deviations for the current mixture components
#'
#' @return A scalar value of the directional derivative D(phi; G)
#'
#' @export
gradient_sigma_univ2 <- function(final, phi, beta, weight, se) {
  n <- nrow(final)
  
  # Denominator: f(x; G) under current mixture (log-normal mixture)
  tmp1 <- normal_mixture2(final, beta, se, weight)
  
  # Numerator: f(x; H) where H is N(beta0, phi^2)
  tmp2 <- normal_cen(final, beta[1], phi)
  
  # Gradient function: D(phi; G) = sum(f_H / f_G) - n
  answer <- sum(tmp2 / tmp1) - n
  
  return(answer)
}


#--------------------------------------------------------------------#
### Function 5. Find Optimal Support Sigma in Gaussian Mixture AFT ###
#--------------------------------------------------------------------#

#' Find new support points (sigma) for Gaussian mixture AFT using gradient-based optimization
#'
#' This function searches for candidate standard deviations (support points) that improve the
#' mixture likelihood using the gradient function D(phi; G). It performs a 3-step process:
#' (1) evaluate directional derivatives, (2) find zero-crossings, and (3) optimize using `nloptr`.
#'
#' @param dat Long-format survival data (counting process style)
#' @param beta Coefficient vector (intercept + others)
#' @param X Character vector of time-dependent covariates
#' @param fig Logical; if TRUE, plot gradient function D(phi; G)
#' @param weight Current weights of mixture components
#' @param bb Current support points (standard deviations of existing components)
#' @param cons Minimum allowable support value
#' @param algorithm Optimization algorithm to use (default: "NLOPT_LN_COBYLA")
#' @param tol Relative tolerance for convergence
#' @param maxiter Maximum number of evaluations
#' @param k Iteration number (used in figure title)
#'
#' @return A list containing:
#' \describe{
#'   \item{newx}{Vector of selected new support points (sigma)}
#'   \item{ttt}{Directional derivative values at selected supports}
#'   \item{max_indicator}{Whether support search was successful ("Success"/"Fail")}
#'   \item{Optimum_algorithm}{Name of optimization algorithm used}
#'   \item{xtol}{Tolerance}
#'   \item{maxiter}{Maximum iterations}
#' }
#' @export
find_sigma_univ2 <- function(dat, beta, X, fig, weight, bb, cons,
                             algorithm = "NLOPT_LN_COBYLA", tol = 1e-5,
                             maxiter = 1000, k) {
  # Step 1: Preprocess
  final <- final_Inf(T0_cal(dat, beta, X))
  logPsi_center <- log(final$psi) - beta[1]
  
  # Construct candidate support set x
  max_range <- max(100 * IQR(logPsi_center), 2 * (diff(range(logPsi_center)))^2)
  min_range <- (min(abs(logPsi_center)) * 0.5)^2
  x <- sort(unique(c(cons,
                     seq(min_range, max_range, length.out = 100),
                     seq(min_range, max_range, length.out = 20),
                     bb)))
  
  # Step 2: Compute directional derivatives and their numerical gradients
  dd <- vapply(x, function(s) gradient_sigma_univ2(final, s, beta, weight, bb), numeric(1L))
  dd_deriv <- vapply(seq_along(x), function(i) {
    h <- 1e-6
    gradient_sigma_univ2(final, x[i] + h, beta, weight, bb) -
      gradient_sigma_univ2(final, x[i] - h, beta, weight, bb)
  }, numeric(1L)) / (2 * 1e-6)
  
  NF <- dd >= 0
  
  # Plot if needed
  if (fig) {
    plot(dd, pch = 19, xlab = "Index of support points", ylab = "D(support)",
         main = paste("Graph of D(S), CNM itr=", k), xaxt = "n")
    axis(1, at = seq(0, length(x), by = 20), labels = round(x[seq(1, length(x), by = 20)], 3))
    abline(h = 0, col = "red", lty = 2, lwd = 2)
  }
  
  # Indicator: flat tail & positive D somewhere
  tail_deriv <- dd_deriv[(length(x) - 6):length(x)]
  Indicator <- ifelse(all(tail_deriv >= 0) & any(NF), 1, 0)
  
  if (Indicator == 1) {
    xvec <- sort(unique(c(
      0.001 * 2^(-6:8),
      cons,
      seq(min_range, max_range, length.out = 100),
      bb
    )))
    x <- xvec
  }
  
  # Filter points with positive gradients
  D_vec <- vapply(x, function(s) gradient_sigma_univ2(final, s, beta, weight, bb), numeric(1L))
  valid_idx <- which(D_vec > 0)
  x <- x[valid_idx]
  D_vec <- D_vec[valid_idx]
  
  if (length(x) == 0) {
    return(list(newx = NULL, ttt = 0, max_indicator = "Success",
                Optimum_algorithm = algorithm, xtol = tol, maxiter = maxiter))
  }
  
  # Step 3: Refine candidate region
  if (length(x) == 1) {
    x <- sort(unique(c(seq(1e-10, x, length.out = 5), seq(x, x + 2, length.out = 50))))
  } else {
    x <- sort(unique(c(seq(1e-10, min(x), length.out = 5),
                       x, seq(min(x), max(x) + 2, length.out = 50))))
  }
  
  # Step 4: Find zero crossings of gradient (local maxima)
  br <- vapply(x, function(s) {
    h <- 1e-6
    (gradient_sigma_univ2(final, s + h, beta, weight, bb) -
        gradient_sigma_univ2(final, s - h, beta, weight, bb)) / (2 * h)
  }, numeric(1L))
  
  crossing_idx <- which(diff(sign(br)) < 0)
  newx <- (x[crossing_idx] + x[crossing_idx + 1]) / 2
  ran <- cbind(x[crossing_idx], x[crossing_idx + 1])
  
  # Step 5: Optimize to refine local maxima
  grad <- function(pointer) {
    tmp1 <- normal_mixture2(final, beta, bb, weight)
    tmp2 <- normal_cen(final, beta[1], pointer)
    -sum(tmp2 / tmp1) + nrow(final)
  }
  
  final_can <- c()
  for (i in seq_along(newx)) {
    res <- tryCatch({
      nloptr(
        x0 = newx[i],
        eval_f = grad,
        lb = ran[i, 1], ub = ran[i, 2],
        opts = list(algorithm = algorithm, xtol_rel = tol, maxeval = maxiter)
      )
    }, error = function(e) NULL)
    
    if (!is.null(res) && res$objective < -1e-7 &&
        min(abs(res$solution / bb - 1)) > 1e-6) {
      final_can <- c(final_can, res$solution)
    }
  }
  
  if (length(final_can) == 0) {
    return(list(newx = NULL, ttt = 0, max_indicator = "Success",
                Optimum_algorithm = algorithm, xtol = tol, maxiter = maxiter))
  }
  
  # Step 6: Validate final candidates
  final_can <- unique(round(final_can, 8))
  ttt <- vapply(final_can, function(s) gradient_sigma_univ2(final, s, beta, weight, bb), numeric(1L))
  newxx <- final_can[ttt > 1e-7]
  ttt <- ttt[ttt > 1e-7]
  
  # Add minimum if needed
  g_cons <- gradient_sigma_univ2(final, cons, beta, weight, bb)
  if (g_cons > 0) {
    newxx <- c(newxx, cons)
    ttt <- c(ttt, g_cons)
  }
  
  if (Indicator == 1) {
    g_tail <- gradient_sigma_univ2(final, max(xvec), beta, weight, bb)
    newxx <- c(newxx, max(xvec))
    ttt <- c(ttt, g_tail)
  }
  
  list(newx = unique(newxx), ttt = ttt,
       max_indicator = ifelse(Indicator == 0, "Success", "Fail"),
       Optimum_algorithm = algorithm, xtol = tol, maxiter = maxiter)
}


#--------------------------------------------------------------#
### Function 6. Mixture Weight Update via Armijo Line Search ###
#--------------------------------------------------------------#

#' Update weights in Gaussian mixture AFT via Armijo rule
#'
#' @param dat Long-format dataset
#' @param X Names of time-dependent covariates
#' @param bb Current support points (sigma)
#' @param weight Current weights
#' @param beta Coefficient vector
#' @param cons Minimum allowable sigma value
#' @param w_criteria Threshold to discard small weights
#' @param itr_n Maximum number of iterations
#'
#' @return List with updated weights and support points
#' @export
update_weight_univ2 <- function(dat, X, bb, weight, beta,
                                cons, w_criteria = 1e-20, itr_n = 200) {
  # Step 1: Preprocessing
  final <- final_Inf(T0_cal(dat, beta, X))
  n <- nrow(final)
  old_weight <- weight
  S <- matrix(0, nrow = n, ncol = length(bb))
  
  # Step 2: Iterative Weight Update
  for (itr in 1:itr_n) {
    logPsi <- log(final$psi)
    
    # Step 2.1: Compute score matrix
    denom <- normal_mixture2(final, beta, bb, weight)
    S <- vapply(bb, function(sigma) {
      normal_cen(final, beta[1], sigma) / denom
    }, numeric(n))
    
    # Step 2.2: Least squares update (with nonnegativity constraint)
    gamma <- sqrt(n * 1e-6)
    Sc <- rbind(gamma * S, rep(1, length(bb)))
    onec <- c(rep(2 * gamma, n), 1)
    weight_new <- nnls(as.matrix(Sc), onec)$x
    weight_new <- weight_new / sum(weight_new)
    
    # Step 2.3: Armijo rule
    eta <- weight_new - old_weight
    sig <- 0.5
    alpha <- 0.3
    
    old_lik <- sum(log(denom))
    direction <- S %*% eta
    for (k in 0:1000) {
      test_weight <- abs(old_weight + sig^k * eta)
      test_weight <- test_weight / sum(test_weight)
      new_lik <- sum(log(normal_mixture2(final, beta, bb, test_weight)))
      
      armijo_cond <- old_lik + sig^k * alpha * sum(direction)
      if (new_lik >= armijo_cond) {
        weight <- test_weight
        break
      }
    }
    
    # Convergence check
    if (abs(new_lik - old_lik) < 1e-8) {
      break
    } else {
      old_weight <- weight
    }
  }
  
  # Step 3: Filter small weights
  valid <- which(weight > w_criteria)
  weight <- weight[valid]
  bb <- bb[valid]
  
  return(list(
    weight = weight,
    bb = bb,
    iteration = ifelse(itr == itr_n, "Fail", "Success"),
    tot_n = itr_n,
    updatedWeight_iter = itr
  ))
}


#-----------------------------------------------------------------------------#
### Function 7. Update Regression Coefficients via Likelihood Maximization ###
#----------------------------------------------------------------------------#

#' Update regression coefficients (beta) in Gaussian mixture AFT model
#'
#' This function maximizes the log-likelihood to estimate regression coefficients
#' \code{beta}, given fixed mixture support points (sigmas) and weights.
#' It solves: \eqn{\hat{\beta} = \arg\max_\beta \ell(\beta; \text{G})}
#'
#' @param dat Long-format dataset
#' @param X Character vector of time-dependent covariate names
#' @param bb Numeric vector of support points (sigmas)
#' @param weight Mixture weights corresponding to \code{bb}
#' @param beta Initial values of regression coefficients
#' @param cons Minimum allowed sigma value (unused here but kept for compatibility)
#' @param algorithm Optimization algorithm used (default: "NLOPT_LN_COBYLA")
#' @param tol Convergence tolerance (default: 1e-5)
#' @param maxiter Maximum number of evaluations (default: 1000)
#'
#' @return A list containing:
#'   \item{upBeta}{Updated regression coefficient vector}
#'   \item{likelihood}{Final log-likelihood value}
#'   \item{itr}{Number of iterations}
#'   \item{status}{Integer code indicating solver status}
#'   \item{message}{Message string from optimizer}
#'
#' @export
myfun_univ2 <- function(dat, X, bb, weight, beta, cons,
                        algorithm = "NLOPT_LN_COBYLA", tol = 1e-5, maxiter = 1000) {
  
  # Step 1: Set search bounds for beta
  bound <- rep(10 * max(abs(beta)), length(beta))
  
  # Step 2: Define negative log-likelihood function
  objective_fn <- function(theta) {
    final <- final_Inf(T0_cal(dat, theta, X))
    -sum(log(normal_mixture2(final, theta, bb, weight)))
  }
  
  # Step 3: Optimize using nloptr
  res <- nloptr(
    x0 = beta,
    eval_f = objective_fn,
    lb = -bound, ub = bound,
    opts = list(algorithm = algorithm, xtol_rel = tol, maxeval = maxiter)
  )
  
  # Step 4: Return results
  list(
    upBeta = res$solution,
    likelihood = res$objective,
    itr = res$iterations,
    status = res$status,
    message = res$message
  )
}



































MAFT_TD <- function(dat,X,process=F,fig=F,D=0.00001,
                    CNM_iter=100,EM_itr=3,w_criteria=1.0e-10,algorithm="NLOPT_LN_COBYLA",tol=1.0e-5,maxiter=2000){
  
  #-------------------------------#
  # --- 1. Set initial values --- #
  #-------------------------------#
  dat1 <- as.data.frame(dat%>%group_by(ID)%>%arrange(stop)%>%slice(n()))
  X.mat <- as.matrix(dat1[,X])
  PAFT <- survreg(Surv(dat1[,"stop"],dat1[,"delta"]) ~ X.mat,dist="lognormal")
  
  # (1) initial betas
  beta <- PAFT$coefficients # initial betas
  # (2) initial support set(= sigma)
  gab <- PAFT$scale
  # (3) initial weight, bb, cons
  cons <- min(gab)/5; # minimum support set
  bb <- gab/2.5
  
  weight <- 1
  beta_lik_value <- NA # the likelihood of the AFT with lognormal
  
  #---------------------------#
  # --- 2. Updated weight ----#
  #---------------------------#
  # Before running the CNM algorithm, we update weight once 
  # and get the corresponding the support points
  final <- T0_cal(dat,beta,X);final <- final_Inf(final)
  logPsi <- log(final$psi)
  logPsi_center <- logPsi - beta[1]
  ini_sigma_weight <- update_weight_univ2(dat,X,bb,weight,beta,cons,w_criteria=w_criteria)
  weight <- ini_sigma_weight$weight # get updated weight
  bb <- ini_sigma_weight$bb # the corresponding support points
  
  
  #--------------------------#
  # --- 3. CNM algorithm ----#
  #--------------------------#
  
  History <- NULL
  
  for(k in 1:CNM_iter){
    
    ## Step1. Update G using the CNM algorithm for a fixed beta
    
    # update sigma
    find_re <- find_sigma_univ2(dat,beta,X,fig,weight,bb,cons,
                                algorithm=algorithm,tol=tol,maxiter=maxiter,k)
    new_sigma <- find_re$newx
    ttt <- find_re$ttt
    max_indicator <- find_re$max_indicator
    
    # the History of the Estimates at the i-th CNM itr
    Estimates <- beta
    compo_N <- length(bb)
    
    Result_Est <- c(Estimates,compo_N,min(bb),max(bb),min(weight),max(weight),max(ttt)
                    ,-beta_lik_value,max_indicator)
    names(Result_Est) <- c(paste("beta",c(0:length(X)),sep=""),
                           "compo_N","minSigma","maxSigma","minW","maxW",
                           "D(H;G)","loglikelihood","MaxSupport")
    
    History <- rbind(History,Result_Est)
    rownames(History) <- paste("CNM_iter=",c(1:nrow(History)),sep="")
    
    if(process==1){
      cat("CNM_itr=",k,'\n')
    }
    
    # If the directional derivative(=ttt) has so small(or negative),
    # then the algorithm is stopped, otherwise keep going.
    if(max(ttt) < D & k>1){
      break
    }else{
      bb <- unique(c(bb,new_sigma))
      weight <- rep(1,length(bb))/length(bb) 
      
      ## Step2. Update weight
      up <- update_weight_univ2(dat,X,bb,weight,beta,cons,w_criteria=w_criteria)
      bb <- up$bb
      weight <- up$weight
      
      ## Step3. Update beta using EM algorithms types for fixed G
      
      for(jj in 1:EM_itr){
        
        # updated beta
        up1 <- myfun_univ2(dat,X,bb,weight,beta,cons,
                           algorithm=algorithm,tol=tol,maxiter=maxiter)
        beta <- up1$upBeta
        
        # updated weight using updated beta
        up2 <- update_weight_univ2(dat,X,bb,weight,beta,cons,w_criteria=w_criteria)
        bb <- up2$bb
        weight <- up2$weight
      }
      beta_lik_value <- up1$likelihood
    }
  }
  
  if(process==1){
    cat("\n","===== Finish MAFT algorithm =====","\n","* Number_itr=",nrow(History),"\n",
        "* Convergence=",ifelse(nrow(History)==CNM_iter,"Fail...","Success!"),"\n")
  }
  
  Final_Est <- History[nrow(History),c(1:(length(X)+1))]
  iter_N <- data.frame(CNM_iter=nrow(History),CNM_tot=CNM_iter,EM_itr=EM_itr)
  
  # End <- Sys.time()
  
  # A <- list(Final_Est=Final_Est,History=History,D=D,iter_N=iter_N,Runnign_Time=End-Start)
  A <- list(Final_Est=Final_Est,Support=bb,W=weight,History=History,D=D,iter_N=iter_N)
}


#---------------------------- function part -------------------------------#

## generaion function ##
gen_fun <- function(dat,n_boot){
  boot_ind <- data.frame(ID=sample(1:n_boot, replace=TRUE),ID1=c(1:n_boot))
  data_boot <- merge(boot_ind,dat,by="ID",all.x=T)
  data_boot <- data_boot[order(data_boot$ID1,data_boot$start),]
  data_boot <- data_boot[,-which(names(data_boot)=="ID")]
  names(data_boot) <- names(dat)
  return(data_boot)
}

### Part1. MAFT algorithm

#------------------------------------#
### ----- Function1. T0_cal -----  ###
#------------------------------------#

# In order to calculate the likelihood, we should compute 
# the baseline failure time T0 = Psi
# This function calculates capital psi(T0) and small psi
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates

T0_cal <- function(dat,beta,Z){
  
  ## 1. Calculation of the baseline failulre time T0, psi
  
  ## Term1 : Calculation of the Capital psi 
  # : integrated from 0 to yi(observed time) of exp(-z(u)*beat^T)
  
  AA <- 0
  for(i in 1:length(Z)){
    AA <- -beta[i+1]*dat[,Z[i]]+AA
  }
  
  dat$psi <- (dat$stop-dat$start)*exp(AA)
  final <- aggregate(dat[,c("psi")],list(dat[,"ID"]),FUN=sum)
  names(final) <- c("ID","psi")
  final <- final[order(final$ID),] # ordering
  
  ## Term2 : Calculation of psi : exp(-coef*z(yi))
  last_dat <- as.data.frame(dat%>%group_by(ID)%>%slice(n()))[,c("ID",Z)]
  BB <- 0
  for(i in 1:length(Z)){
    BB <- -beta[i+1]*last_dat[,Z[i]]+BB
  }
  final$psi_d <- exp(BB)
  final <- final[order(final$ID),] # ordering
  
  ## Joint delta
  delta <- as.data.frame(dat%>%group_by(ID)%>%slice(n()))[,c("ID","delta")]
  final <- merge(final,delta,by="ID",all.x=T)
  final <- final[order(final$ID),] # ordering
  
  return(final)
  
}

## If Psi, Psi_d --> Inf then 1e+300

final_Inf <- function(final){
  
  # psi
  if(sum(final$psi>1.0e+300)!=0){
    final[which(final$psi>1.0e+300),"psi"] <- 1.0e+300
  }else{
    final <- final
  }
  if(sum(final$psi<(-1.0e+300))!=0){
    final[which(final$psi<(-1.0e+300)),"psi"] <- -1.0e+300
  }else{
    final <- final
  }
  
  # psi_d
  if(sum(final$psi_d>1.0e+300)!=0){
    final[which(final$psi_d>1.0e+300),"psi_d"] <- 1.0e+300
  }else{
    final <- final
  }
  if(sum(final$psi_d<(-1.0e+300))!=0){
    final[which(final$psi_d<(-1.0e+300)),"psi_d"] <- -1.0e+300
  }else{
    final <- final
  }
  
  ## 0 --> 1.0e-300
  
  # psi
  if(sum(final$psi==0)!=0){
    final[which(final$psi==0),"psi"] <- 1.0e-300
  }else{
    final <- final
  }
  # psi_d
  if(sum(final$psi_d==0)!=0){
    final[which(final$psi_d==0),"psi_d"] <- 1.0e-300
  }else{
    final <- final
  }
  
  return(final)
}



#----------------------------------------#
### ----- Function2. normal_cen -----  ###
#----------------------------------------#

# It is the likelihood function for one component Gaussian distribution
# final <- T0_cal(dat,beta,X)
# beta0 <- beta[1] # intercept : E(logT0)
# se <- 1

normal_cen <- function(final,beta0,se){
  
  # When the events occur, density function apply for the likelihood function
  den <- dnorm((log(final$psi)-beta0)/se)*(final$psi_d/(final$psi*se))
  
  # likelihood function with censoring
  lik <- final$delta*den + (1-final$delta)*(1-pnorm((log(final$psi)-beta0)/se))
  lik <- ifelse(lik<=1.0e-20,1.0e-20,lik)
  
  return(lik)
}

#---------------------------------------------#
### ----- Function3. normal_mixture2 -----  ###
#---------------------------------------------#

# It is the likelihood function for mixture Gaussian distributions
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates
# se <- c(1,3) # Gaussian mixture with two components
# weight <- c(0.8,0.2) # Gaussian mixture with two components

normal_mixture2 <- function(final,beta,se,weight){
  
  # mixture likelihood
  answer <- rep(0,nrow(final))
  for(i in 1:length(se)){
    answer <- answer + weight[i]*(normal_cen(final,beta[1],se[i])) 
  }
  answer1 <- ifelse(answer<=1.0e-20,1.0e-20,answer)
  
  return(answer1)
}



#--------------------------------------------#
### ----- Function4. gradient sigma -----  ###
#--------------------------------------------#

# It means the directional derivative from G to theta : d(theta;G)
# d(theta;G) is known as 'gradient function'
# phi <- 0.8 # candidate for the optimal Standard error
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates
# variance <- c(1,3) # Gaussian mixture with two components
# weight <- c(0.8,0.2) # Gaussian mixture with two components

gradient_sigma_univ2 <- function(final,phi,beta,weight,se){
  
  ## likelihood - directional derivative, D
  n <- nrow(final);tmp <- rep(0,n)
  
  # f(x;G) - denominator, where G : given sigma
  tmp1 <- normal_mixture2(final,beta,se,weight)
  
  # f(x;H) - numerator
  tmp2 <- normal_cen(final,beta[1],phi)
  
  # D : the directional derivative from G to H
  answer <- sum(tmp2/tmp1)-n
  
  return(answer)
}



#----------------------------------------------#
### ----- Function5. find_sigma_univ2 -----  ###
#----------------------------------------------#

# Given the beta, we find the optimization of sigma
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates
# weight <- c(0.8,0.2) # Gaussian mixture with two components
# bb <- c(1,2) # initial of the support points(sigma)
# cons <- 0.001 # minimum of the support points(sigma)

find_sigma_univ2 <- function(dat,beta,X,fig,weight,bb,cons,
                             algorithm="NLOPT_LN_COBYLA",tol=1.0e-5,maxiter=1000,k){
  
  
  #-- 1. find new support --#
  
  ## calculation of T0
  final <- T0_cal(dat,beta,X)
  final <- final_Inf(final)
  logPsi <- log(final$psi)
  logPsi_center <- logPsi - beta[1]
  
  # censoring of the biggest residual change no censoring!
  # --> one of the solving of the optimization problem
  # final[which.max(abs(logPsi_center)),"delta"] <- 1
  
  ## range : sigma set(support set)
  # the optimal newx candidate's support set was randomly assigned by researcher.
  # i.e. it is properly divided.
  max_range <- max(100*IQR(logPsi_center),2*(diff(range(logPsi_center)))^2)
  min_range <- (min(abs(logPsi_center))*0.5)^2
  x <- seq(min_range,max_range,max_range/100)
  x <- c(x,seq(x[1],x[2],(x[2]-x[1])/20))
  
  # x0 <- rep(0,15)
  # for(i in 1:15){
  #   x0[i] <- 0.001*2^(-6+i)
  # }
  # x <- sort(unique(c(x0,cons,x,bb)))
  x <- sort(unique(c(cons,x,bb)))
  
  #-- 2. find a set of all local maximizers of D(directional derivative from G to H) --#
  
  # direction derivative Graph
  dd <- c();dd_deriv <- c();NF <- c()
  for(ii in 1:length(x)){
    dd[ii] <- gradient_sigma_univ2(final,x[ii],beta,weight,bb)
    dd_deriv[ii] <- (gradient_sigma_univ2(final,x[ii]+1.0e-6,beta,weight,bb)-
                       gradient_sigma_univ2(final,x[ii]-1.0e-6,beta,weight,bb))/2.0e-6 
    NF[ii] <- ifelse(dd[ii]<0,F,T)
  }
  
  if(fig==T){
    main_title <- paste("Graph of D(S), CNM itr=",k)
    plot(dd,pch=19,xlab="Index of support points",ylab="D(support)",main=main_title,xaxt="n")
    x_axis <- c(0,20,40,60,80,100,120)
    axis(side=1,at=x_axis,labels=round(c(x[1],x[20],x[40],x[60],x[80],x[100],x[120]),3))
    abline(h=0,col="red",lty=2,lwd=2)
  }
  # main_title <- paste("Graph of D(S), CNM itr=",k)
  # plot(dd,pch=19,xlab="Index of suppor points",ylab="D(support)",main=main_title,xaxt="n")
  # x_axis <- c(0,20,40,60,80,100,120)
  # axis(side=1,at=x_axis,labels=round(c(x[1],x[20],x[40],x[60],x[80],x[100],x[120]),3))
  # abline(h=0,col="red",lty=2,lwd=2)
  # dev.copy(jpeg,filename=paste(main_title,".jpg",sep=""), res=300, height=20, width=28,  units='cm')
  # dev.off()
  
  
  ## The finding the optimal value equal means that
  ## The all the derivative have positive and the positive value have some 
  ## sum(diff(dd_deriv)<0)==0 & length(dd[NF])!=0
  ## So, we find the maximum residual and change the no censoring!
  
  ## First, we check the function of D shape using the residuals !!
  Indicator <- ifelse(sum(dd_deriv[(length(x)-6):length(x)]<0)==0 & length(dd[NF])!=0,1,0)
  
  if(Indicator==1){
    max_range <- max(100*IQR(logPsi_center),2*(diff(range(logPsi_center)))^2)
    min_range <- (min(abs(logPsi_center))*0.5)^2
    x <- seq(min_range,max_range,max_range/100)
    x <- c(x,seq(x[1],x[2],(x[2]-x[1])/20))
    
    x0 <- rep(0,15)
    for(i in 1:15){
      x0[i] <- 0.001*2^(-6+i)
    }
    x <- sort(unique(c(x0,cons,x,bb)))
    xvec <- x
  }else{
    x <- x
  }
  
  
  # select gradient values > 0
  gg <- c()
  for(ii in 1:length(x)){
    if(gradient_sigma_univ2(final,x[ii],beta,weight,bb)<0){
      gg[ii] <- NA
      x[ii] <- NA
    }else{
      gg[ii] <- gradient_sigma_univ2(final,x[ii],beta,weight,bb)
      x[ii] <- x[ii]
    }
  }
  x <- as.numeric(na.omit(x))
  
  ## According to the number of support points
  if(length(x)==0){
    newxx <- NULL
    ttt <- 0
  }else{
    
    # if x has only one component, then we have two side values
    if(length(x)==1){
      x1 <- x
      x <- sort(unique(c(seq(1.0e-10,x1,length=5),seq(x1,x1+2,length=50))))
    }else{
      x1 <- x
      x <- sort(unique(c(seq(1.0e-10,min(x1),length=5),x1,seq(min(x1),max(x1)+2,length=50))))
    }
    
    
    # br : the derivative value
    # newx <- local maximum value 
    # ran <- range of the newx(it's br values at newx)
    br <- c();newx <- c();ran <- c()
    
    ## the derivation of the gradient function by given the new support points
    for(ii in 1:length(x)){
      br[ii] <- (gradient_sigma_univ2(final,x[ii]+1.0e-6,beta,weight,bb)-
                   gradient_sigma_univ2(final,x[ii]-1.0e-6,beta,weight,bb))/2.0e-6 
      # find the value where the derivative is zero
      if(ii>1&&br[ii-1]>0&&br[ii]<0){
        newx <- c(newx,mean(c(x[ii-1],x[ii])))    
        ran <- rbind(ran,c(x[ii-1],x[ii]))
      }
    }
    
    
    # If the last value of the br is postive, we fail the finding the optimal value
    # So, we check using the indicator variable : max_indicator
    # If max_indicator = 1 then the last of the br is postive and the finding is fail...
    # If max_indicator = 1 then we replace the newx with maximum value of the x(candidate set)
    # and range is x[2] ~ max(x)
    # max_indicator <- ifelse(br[length(x)]>0,1,0)
    # if(max_indicator==1){
    #   newx <- c(newx,max_range+1)
    #   ran <- rbind(ran,c(x[2],max(max(x),max(newx))))
    # }
    # if(max_indicator==1){
    #   newx <- c(newx,10)
    #   ran <- rbind(ran,c(x[2],10))
    # }
    
    
    #-- 3. find support within a given interval - final candidates --#
    
    # We use nloptr(optimal algorithm) to find the optimal value of the newx once more.
    # i.e. find support within a given interval - final candidates
    # nloptr using grad function
    # grad function has the negative the directional derivative (-d(H;G))
    
    grad <- function(pointer){
      
      ## likelihood - directional derivative, D
      n <- dim(final)[1];tmp <- rep(0,n)
      
      # f(x;G) - denominator, where G : given sigma
      tmp1 <- normal_mixture2(final,beta,bb,weight)
      
      # f(x;H) - numerator
      tmp2 <- normal_cen(final,beta[1],pointer)
      
      answer <- sum(tmp2/tmp1)-n
      
      return(-answer)
    }
    
    
    final_can <- c() # final candidates
    
    for(t in 1:length(newx)){
      
      res1 <- nloptr(x0=newx[t],eval_f=grad,lb = ran[t,1],ub = ran[t,2],
                     opts = list("algorithm"=algorithm,"xtol_rel"=tol,maxeval=maxiter)) 
      
      fval <- res1$objective # the directional derivative(=nloptr valule)
      answer <- res1$solution # newx using by results of the nloptr
      
      
      # we can use the another optimum algorhtm  : optim 
      # res1 <- optim(newx[t],grad,lower=ran[t,1],upper=ran[t,2],method = "L-BFGS-B")
      # fval <- res1$value;answer <- res1$par
      
      
      # the optimal support(= final candidate) is satisfied
      # when the directional derivative(=fval) has negative.
      # fval means the netative d(H;G). So, the negative fval means that d(H;G) has postive. 
      # fval< -0.0001 && min(abs(answer/bb-1))!=0 --> we didn't satisfied the optimal value.
      # otherwise, the algorithm finised!.
      
      # fval< -0.0001
      if(fval< -0.0000001 && min(abs(answer/bb-1))!=0){
        final_can <- c(final_can,answer)
      }else{
        final_can <- final_can
      }
    }
    
    # If final_can didn't satistifed the optimal value(length(final_can)!=0),
    # we check again that the final_can has positive d(H;G).
    # Also, even if d(H;G) of the cons(minimum) has positive, 
    # it is added as a newx valule
    
    if(length(final_can)!=0){
      newx <- final_can
      newx <- unique(newx)
      
      newxx <- c();ttt <- c()
      for(jj in 1:length(newx)){
        gg <- gradient_sigma_univ2(final,newx[jj],beta,weight,bb)
        ttt <- c(ttt,gg)
        if(gg>0.0000001){newxx <- c(newxx,newx[jj])}
      }
      if(gradient_sigma_univ2(final,cons,beta,weight,bb)>0){
        ttt <- c(ttt,gradient_sigma_univ2(final,cons,beta,weight,bb))
        newxx <- c(newxx,cons)
      }
    }else{# else, we finished the algorithm with ttt=0
      newxx <- newx;
      ttt <- 0
    }
  }
  
  if(Indicator==1){
    newxx <- c(newxx,max(xvec))
    max_ttt <- gradient_sigma_univ2(final,max(xvec),beta,weight,bb)
    ttt <- c(ttt,max_ttt)
  }else{
    newxx <- newxx
  }
  
  result <- list(newx = unique(newxx),ttt=ttt,
                 max_indicator=ifelse(Indicator==0,"Success","Fail"),
                 Optimum_algorithm=algorithm,xtol=tol,maxiter=maxiter)
  return(result)
}


#-------------------------------------------------#
### ----- Function6. update_weight_univ2 -----  ###
#-------------------------------------------------#

# Given beta and updated sigma, we update the weights
# using 'Armijo rule'
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates
# weight <- c(0.8,0.2) # Gaussian mixture with two components
# bb <- c(1,2) # initial of the support points(sigma)
# cons <- 0.001 # minimum of the support points(sigma)

update_weight_univ2 <- function(dat,X,bb,weight,beta,cons,w_criteria=1.0e-20,itr_n=200){
  
  
  # ---- 1. Setting ----#
  final <- T0_cal(dat,beta,X);final <- final_Inf(final)
  n <- nrow(final)
  old_weight <- weight # before updated weight
  S <- matrix(rep(0,n*length(bb)),ncol=length(bb)) # score matrix
  W <- c()
  
  
  # ---- 2. Running algorithm ----#
  
  for(itr in 1:itr_n){   
    
    logPsi <- log(final$psi)
    
    # score matrix
    for(j in 1:length(bb)){
      S[,j] <- normal_cen(final,beta[1],bb[j])/normal_mixture2(final,beta,bb,weight)
    }
    
    gamma <- sqrt(n*10^(-6)) # gmma value
    Sc <- rbind(gamma*S,1)
    onec <- matrix(c(2*gamma*rep(1,n),1),ncol=1)
    
    ## update weight 
    # min|weight -1|^2 + gamma*|S*weight-2|^2, subject to weight>=0
    # where we set that the gamma valule is sqrt(n*10^(-6))
    # the solution to the least squares problem with non-negativity constraints
    weight <- nnls(as.matrix(Sc), as.vector(onec))$x
    weight <- weight/sum(weight)
    
    #--- Armijo rule ---#
    # To ensure the monotonic increase of log-likelihood at each iteration, 
    # we use the 'back-tracking line search strategy' guarded by 'Armijo ruel' : 
    
    eta <- weight - old_weight
    
    # In the study, we sets sigma=0.05
    # and the value of alpha is often chosen to be small
    # so, we use the alpha=0.3
    
    sig <- 1/2; alpha <- 0.3; 
    
    old_lik <- sum(log(normal_mixture2(final,beta,bb,old_weight)))
    
    for(k in 0:1000){   
      test_weight <- abs(old_weight+sig^k*eta)
      test_weight <- test_weight/sum(test_weight)
      new_lik <- sum(log(normal_mixture2(final,beta,bb,test_weight)))
      if(new_lik>=old_lik+matrix(rep(1,n),nrow=1)%*%S%*%matrix(eta,ncol=1)*sig^k*alpha){
        break
      }
    }
    
    # we check once more whether the new likelihood converges
    if(sqrt((new_lik-old_lik)^2)<1.0e-8){   
      break
    }else{
      old_weight <- weight
    }
  }
  
  # the number of iterations in the algorithm
  upWeight_iter <- itr 
  
  # If updated weight < w_criteria then we discard the updated weight and 
  # the corresponding sigma.
  ww <- which(weight<=w_criteria)
  bb[ww] <-  NA;bb <- as.vector(na.omit(bb))
  weight[ww] <- NA;weight <- as.vector(na.omit(weight))
  result <- list(weight=weight,bb=bb,iteration=ifelse(upWeight_iter==itr_n,"Fail","Success"),
                 tot_n=itr_n,updatedWeight_iter=upWeight_iter)
  
  return(result)
}

#-----------------------------------------#
### ----- Function7. myfun_univ2 -----  ###
#-----------------------------------------#
# Given updated G measures, we update the beta
# beta <- c(1,rep(0.5,length(X))) # the coefficients corresponding to the covariates
# weight <- c(0.8,0.2) # Gaussian mixture with two components
# cons <- 0.001 # minimum of the support points(sigma)

myfun_univ2 <- function(dat,X,bb,weight,beta,cons,
                        algorithm="NLOPT_LN_COBYLA",tol=1.0e-5,maxiter=1000){
  
  # bound
  bound <- rep(10*max(abs(beta)),length(beta))
  
  # loglikelihood for beta given by updated G components
  myfun2 <- function(theta){
    final <- T0_cal(dat,theta,X);final <- final_Inf(final)
    g = -sum(log(normal_mixture2(final,theta,bb,weight)))
    return(g)
  }
  
  res2 <- nloptr(x0=beta, eval_f=myfun2,lb = -bound, ub = bound,
                 opts=list("algorithm"=algorithm,"xtol_rel"=tol,"maxeval"=maxiter)) 
  
  # Results
  likelihood <- res2$objective # likelihood
  updatedBeta <- res2$solution # updated beta
  itr <- res2$iterations # number of iterations of the algorithm
  status <- res2$status # integer value with the status of the optimization
  message <- res2$message # more informative message with the status of the optimization
  
  result <- list(upBeta=updatedBeta,likelihood=likelihood,itr=itr,status=status,message=message)
  
  return(result)
}


#-------------------------------------------------------------#
### Function B1. Bootstrap Sample Generator for MAFT Estimation
#-------------------------------------------------------------#
#' Generate a bootstrap sample by sampling subject IDs with replacement.
#'
#' @param dat A data.frame in counting process format (with ID, start, stop, delta, covariates).
#' @param n_boot Number of subjects to resample (usually equal to number of unique IDs).
#' @return A bootstrap sample of the same format as `dat`.
#' @examples
#' boot_dat <- gen_fun(dat, n_boot = 200)
gen_fun <- function(dat, n_boot){
  # Sample subject IDs with replacement
  boot_ind <- data.frame(ID = sample(1:n_boot, replace = TRUE),
                         ID1 = seq_len(n_boot))
  
  # Merge to generate bootstrapped data by ID
  data_boot <- merge(boot_ind, dat, by = "ID", all.x = TRUE)
  
  # Reorder to maintain counting process format
  data_boot <- data_boot[order(data_boot$ID1, data_boot$start), ]
  
  # Remove temporary ID1 column and rename back
  data_boot$ID <- data_boot$ID1
  data_boot$ID1 <- NULL
  
  return(data_boot)
}

