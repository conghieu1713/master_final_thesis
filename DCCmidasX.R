###################################################################
###                                                             ###
#                     DCC-MIDAS-X FITTING                         #
###                                                             ###
###################################################################

# Author: Cong-Hieu, NGUYEN (524102110660)
# MSc. Candidate, University of Economics HCMC


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== RESET ALL =====###
  rm(list = ls(all.names = TRUE))
  graphics.off()
  invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)), 
                 detach, character.only = TRUE, unload = TRUE, force = TRUE))


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== LIBRARY DECLARATION =====###
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(dccmidas)
  library(rumidas)
  library(ggplot2)
  library(patchwork)
  library(zoo)
  library(xts)
  library(Rsolnp)

#######======================================================#######
#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== GET DATA FUNCTION =====###

file_path <- "D:\\Post-graduate programme\\Final Thesis\\Model\\official.xlsx" 

get_bivariate_data <- function(target_asset, file_path, K=12) {
  # ============================================
  # IMPORT DATA FROM EXCEL FILE, using guess_max
  # ============================================
    df_rate <- read_excel(file_path, sheet = "rate", guess_max = 100000) %>% 
      mutate(Date = as.Date(Date))
    df_rv_month <- read_excel(file_path, sheet = "rv_month", guess_max = 100000) %>% 
      mutate(Date = as.Date(Date))
  # ============================================
  # DYNAMIC BUFFER
  # Avoid Error related to Hessian (At least one parameter must not be fixed)
  # Omit NA value
  # ============================================
    valid_rv_data <- df_rv_month %>% 
      select(Date, VNIndex, all_of(target_asset)) %>% 
      drop_na()
  # Actual Begin Date
    start_mv <- min(valid_rv_data$Date)
  # Buffer for K
    start_rate_limit <- start_mv %m+% months(K+1)
  # =============================================
  # HIGH-FREQ DATA
  # =============================================
    rate_pair <- df_rate %>% 
      select(Date, VNIndex, all_of(target_asset)) %>% 
      filter(Date >= start_rate_limit) %>% 
      drop_na()
  
  # RETAIN COLUMN NAMES: use drop = FALSE to avoid error "replacement has length zero"
    xts_vn <- xts(rate_pair[, "VNIndex", drop = FALSE], order.by = rate_pair$Date)
    xts_target <- xts(rate_pair[, target_asset, drop = FALSE], order.by = rate_pair$Date)
  
  # r_t LIST: Avoid error "the condition has length > 1" of dcc_fit
  # Add white noise (sd = 1e-6) for SJC
    if (target_asset == "SJC") {
      set.seed(123) # Fix the seed so that the results don't change after each run
      noise <- rnorm(length(xts_target), mean = 0, sd = 1e-6)
      xts_target <- xts_target + noise
    }
    r_t_list <- list(xts_vn, xts_target)
  
  # ===========================================
  # LOW-FREQ DATA
  # ===========================================
    mv_xts <- xts(valid_rv_data[, -1], order.by = valid_rv_data$Date)
  # MV-INTO-MAT
    mv_vn <- mv_into_mat(r_t_list[[1]], mv_xts[, "VNIndex"], K = K, type = "monthly")
    mv_target <- mv_into_mat(r_t_list[[2]], mv_xts[, target_asset], K = K, type = "monthly")
  
  Z_list <- list(mv_vn, mv_target)
  
  # ===========================================
  # RETURN
  # ===========================================
  return(list(
    r_t = r_t_list,
    MV = Z_list,
    dates = index(xts_vn),
    target = target_asset,
    K_used = K
  ))
}


#######======================================================#######
#######======================================================#######
#######======================================================#######



###===== DCC-MIDAS FITTING AND RESULTS FUNCTION =====###


run_dcc_estimation <- function(target_asset, file_path, 
                               K = 36, N_c = 21, K_c = 504, 
                               univ_model = "GM_noskew") {
  # ====================================================
  # IMPORT DATA, using get_asset_data
  # ====================================================
    prepared_data <- get_bivariate_data(target_asset, file_path, K = K)
  
  if (is.null(prepared_data)) {
    cat(">> [ERROR] Data cannot be processed for", target_asset, "\n")
    return(NULL)
  }
  
  # ====================================================
  # MODEL FITTING
  # ====================================================
  cat(sprintf("\n[ESTIMATION] DCC-MIDAS for %s | K=%d, N_c=%d, Kc=%d, Model=%s\n", 
              target_asset, K, N_c, K_c, univ_model))
  
  fit <- tryCatch({
    dcc_fit(
      r_t = prepared_data$r_t, 
      MV = prepared_data$MV, 
      univ_model = univ_model, 
      corr_model = "DCCMIDAS", 
      K = prepared_data$K_used,
      N_c = N_c, 
      K_c = K_c, 
      lag_fun = "Beta",
      distribution = "norm"
    )
  }, error = function(e) {
    cat(">> [ESTIMATION ERROR]:", e$message, "\n")
    return(NULL)
  })
  
  if (!is.null(fit)) {cat(">> [ESTIMATION] SUCCESS!.\n")}
  
  # ====================================================
  # RETURN
  # ====================================================
  return(list(
    fit = fit,               
    data = prepared_data
  ))
}


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== EXTRACT NORMALIZED RESIDUALS FROM GARCH-MIDAS FUNCTION =====###
extract_std_residuals <- function(estimation_result) {
  # Extracting Components from Estimation Results
    fit <- estimation_result$fit
    data_info <- estimation_result$data
  # Asset name to name the column for the results
    target_name <- data_info$target
  # GET RETURNS - r_t
  # Convert to numeric form to ensure consistent calculations
    r_vn <- as.numeric(data_info$r_t[[1]])
    r_target <- as.numeric(data_info$r_t[[2]])
  # CONDITIONAL VOLATILITY - sqrt(h_t)
  # H_t in dccmidas contains variance
    vol_vn <- as.numeric(sqrt(fit$H_t[1, 1, ]))
    vol_target <- as.numeric(sqrt(fit$H_t[2, 2, ]))
  # CALCULATING STANDARDIZED RESIDUALS
  # eps = r_t / sqrt(h_t)
    std_resid_vn <- r_vn / vol_vn
    std_resid_target <- r_target / vol_target
  # RETURN
  # Return as a data.frame with the date for easier management.
    std_resid_df <- data.frame(
      Date = data_info$dates,
      VNIndex = std_resid_vn,
      Target = std_resid_target
    )
  # Rename the Target column to the actual asset name (BTC, XAUUSD, ...)
    colnames(std_resid_df)[3] <- target_name
  
    cat(sprintf(">> The %d observed normalized residuals for VNIndex and %s have been extracted.", 
              nrow(std_resid_df), target_name))
  
  return(std_resid_df)
}


#######======================================================#######
#######======================================================#######
#######======================================================#######



###===== PREPARE DATA FOR FITTING DCC-MIDAS-X =====###
prepare_dcc_x_data <- function(estimation_result, file_path, sheet_name, shock_name, K_x) {
  
  # READ MACROECONOMIC DATA DIRECTLY FROM EXCEL
    macro_df <- read_excel(file_path, sheet = sheet_name, guess_max = 100000) %>%
      mutate(Date = as.Date(Date))
  
  # Check the necessary columns.
    if (!"Date" %in% colnames(macro_df)) {
      stop("Error: Sheet '", sheet_name, "' must have columm 'Date'.")
    }
    if (!shock_name %in% colnames(macro_df)) {
      stop("Error: Varialbe '", shock_name, "' was not found.")
    }
  
  # EXTRACT STANDARDIZED RESIDUES (Call the extract_std_residuals function again)
    std_resid_df <- extract_std_residuals(estimation_result)
    res_matrix <- as.matrix(std_resid_df[, 2:3]) # Get the VNIndex and Target columns.
    rownames(res_matrix) <- as.character(std_resid_df$Date)
  
  # CREATE A MIDAS MATRIX FOR MACRO VARIABLES
    shock_xts <- xts(as.numeric(macro_df[[shock_name]]), order.by = macro_df$Date)
    X_mat <- mv_into_mat(
      x = estimation_result$data$r_t[[1]], 
      mv = shock_xts,  
      K = K_x, 
      type = "monthly"
    )
  
  # DATA ALIGNMENT
  # Ensure the time axis of the daily residuals and the monthly macrovariable align.
    X_mat_temp <- t(X_mat)
    daily_dates <- as.Date(index(estimation_result$data$r_t[[1]]))
    row.names(X_mat_temp) <- as.character(daily_dates)
    common_dates <- intersect(row.names(res_matrix), row.names(X_mat_temp))
    res_matrix_final <- res_matrix[common_dates, ]
    X_mat_final <- X_mat_temp[common_dates, ]
    cat(sprintf(">> [COMPLETED] %d observations have been prepared for variable %s.\n", 
                nrow(res_matrix_final), shock_name))
  
  # Returns a list containing two clean matrices to load into solnp.
  return(list(
    res_matrix = res_matrix_final,
    X_mat = X_mat_final
  ))
}


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== DECLARING THE MIDAS WEIGHT FUNCTION =====###
# Beta Weight Function (with a fixed parameter w1 = 1 for decreasing weights)
  beta_weights <- function(K, w) {
    k <- 1:K
    weights <- (1 - k/(K+1))^(w - 1)
    weights <- weights / sum(weights)
    return(weights)
  }


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== LOG-LIKELIHOOD FUNCTION FOR DCC-MIDAS-X =====###
  dcc_midas_x_loglik <- function(params, res_data, X_data) {
    a <- params[1]
    b <- params[2]
    m <- params[3]
    theta <- params[4]
    w <- params[5]
    # --- PROTECTION 1: Parameter Constraints ---
    # If numDeriv attempts to move the parameter outside the safe zone, return a huge penalty.
      if (a <= 0 || b <= 0 || (a + b) >= 0.999 || w <= 1.001) {
        return(1e6) 
      }
      T_obs <- nrow(res_data)
      K_lags <- ncol(X_data)
      weights <- beta_weights(K_lags, w)
      Z_t <- m + theta * (X_data %*% weights)
      rho_bar_t <- (exp(2 * Z_t) - 1) / (exp(2 * Z_t) + 1)
      Q_prev <- cov(res_data)
      LL <- 0
      
      for (t in 2:T_obs) {
        R_bar_t <- matrix(c(1, rho_bar_t[t], 
                            rho_bar_t[t], 1), nrow = 2, ncol = 2)
        
        eps_prev <- matrix(res_data[t-1, ], ncol = 1)
        Q_t <- (1 - a - b) * R_bar_t + a * (eps_prev %*% t(eps_prev)) + b * Q_prev
      
      # --- PROTECTION 2: Prevent matrix Q_t degeneration ---
        if (any(is.na(Q_t)) || Q_t[1,1] <= 0 || Q_t[2,2] <= 0) return(1e6)
        D_t <- diag(1 / sqrt(diag(Q_t)))
        R_t <- D_t %*% Q_t %*% D_t
        eps_t <- matrix(res_data[t, ], ncol = 1)
        det_R <- det(R_t)
      
      # --- PROTECTION 3: Matrix Determinant R_t ---
      # Instead of assigning 1e-8, we penalize directly to change the Hessian direction.
        if (is.na(det_R) || det_R <= 0) return(1e6) 
        inv_R <- solve(R_t)
        LL_t <- -0.5 * (log(det_R) + t(eps_t) %*% inv_R %*% eps_t - t(eps_t) %*% eps_t)
        LL <- LL + LL_t
        Q_prev <- Q_t
    }
    return(-as.numeric(LL)) 
  }


#######======================================================#######
#######======================================================#######
#######======================================================#######

 
 
###===== DCC-MIDAS-X FITTING AND RESULTS =====###

 
# ----------------------------------------------------------------- #
# --- 0. The lag values are used when fitting DCC-MIDAS-X model --- #
# ----------------------------------------------------------------- #
 
# --- Gold (XAUUSD, GC=F, SJC)
#       K   = 24 or 36 (or 12 if absolutely necessary)
#       N_c = 21
#       K_c = 126 or 252 or 504

# --- Crypto (BTC, ETH, BNB)
#       K   = 12 or 24
#       N_c = 21
#       K_c = 126 or 252

# --- Shocks (GPR, GPRT, GPRA, GEPU_current, GEPU_ppp, WUI)
#       K_x = 12 or 24 or 36

# -------------------------------------------------------------------- #
# --- 1. Extract Normalized Residuals from uni-variate GARCH-MIDAS --- #
# -------------------------------------------------------------------- #
  # Setting up input parameters
    target_asset = "BTC"; shock_name = "GPRA"
    K <- 12; N_c <- 21; K_c <- 252; K_x <- 24
    univ_model = "GM_noskew"
  # Uni-variate GARCH-MIDAS fitting
    res <- run_dcc_estimation(
       target_asset = target_asset,
       file_path = file_path,
       K = K, N_c = N_c, K_c = K_c,
       univ_model = univ_model
    )
  # Extract Residuals
    dccx_data <- prepare_dcc_x_data(
       estimation_result = res,
       file_path = file_path,
       sheet_name = "shocks",
       shock_name = shock_name,
       K_x = K_x
    )
  # Scale X_mat (if necessary)
    # --- Z-score scale
    # dccx_data$X_mat <- scale(dccx_data$X_mat)
    # --- [0,1] scale
    dccx_data$X_mat <- apply(dccx_data$X_mat, 2, scales::rescale, to = c(0, 1))
 
# ------------------------------------------------------------------- #
# ----------------------- 2. DCC-MIDAS-FITTING ---------------------- #
# ------------------------------------------------------------------- #
  # Set up initial values for the parameters
    par_init <- c(a = 0.05, b = 0.85, m = 0.15, theta = 0.2, w = 2)
  # Declaring conditional functions
    ineq_fun <- function(params, res_data, X_data) {
      return(params[1] + params[2]) 
    }
  # Optimization
   opt_results <- solnp(
     pars = par_init,                        # Initial values for pả
     fun = dcc_midas_x_loglik,               # Log-likelihood function
     ineqfun = ineq_fun, 
     ineqLB = c(0),                          # a + b > 0
     ineqUB = c(0.999),                      # a + b < 1
     LB = c(1e-3, 1e-3, -10, -10, 1.001),    # Lower bound for (a, b, m, theta, w)
     UB = c(0.2, 0.999, 10, 10, 50),         # Upper bound for (a, b, m, theta, w)
     res_data = dccx_data$res_matrix,        # Load the residual matrix
     X_data = dccx_data$X_mat                # Load macro matrix (X)
   )

# ------------------------------------------------------------------- #
# ------------ 3. CALCULATE STANDARD ERROR AND P VALUE -------------- #
# ------------------------------------------------------------------- #
  report_dcc_midas_results <- function(opt_res, res_matrix, x_matrix) {
    # 1. Extract basic parameters
      params <- opt_res$pars
      n_param <- length(params)
      n_obs <- nrow(res_matrix)
  
    # --- Corresponding parameter name (Automatically adapts if you lock 'a')
      if (n_param == 5) {
          param_names <- c("a (alpha)", "b (beta)", "m (mu)", "theta", "w (omega)")
      } else {
          param_names <- c("b (beta)", "m (mu)", "theta", "w (omega)")
      }
  
    # 2. Calculating the Hessian Matrix and Standard Error (Std. Error) 
    # --- Using ginv (Moore-Penrose) to process the nearly singular matrix
      h <- numDeriv::hessian(func = dcc_midas_x_loglik, x = params, 
                          res_data = res_matrix, X_data = x_matrix)
  
    # --- Add a small buffer layer (ridge) to ensure stability when inverted.
      cov_mat <- MASS::ginv(h + diag(1e-6, n_param))
      std_err <- sqrt(abs(diag(cov_mat)))
  
    # 3. Calculate t-stat và P-value
      t_stat <- params / std_err
      p_val <- 2 * (1 - pnorm(abs(t_stat)))
  
    # 4. Calculate information criteria index (Model Fit)
      ll <- -tail(opt_res$values, 1)
      aic <- 2 * n_param - 2 * ll
      bic <- n_param * log(n_obs) - 2 * ll
  
    # 5. Create Table of result
      results_df <- data.frame(
        Parameter = param_names,
        Estimate = round(params, 6),
        `Std. Error` = round(std_err, 6),
        `t-stat` = round(t_stat, 4),
        `P-value` = round(p_val, 4),
        Signif = ifelse(p_val < 0.01, "***", ifelse(p_val < 0.05, "**", ifelse(p_val < 0.1, "*", "")))
      )
  
    # 6. Print the table to the Console in academic style
      cat("\n==========================================================\n")
      cat("      ESTIMATION RESULTS: DCC-MIDAS-X MODEL\n")
      cat("==========================================================\n")
      print(results_df, row.names = FALSE)
      cat("----------------------------------------------------------\n")
      cat(sprintf("Log-Likelihood: %0.4f\n", ll))
      cat(sprintf("AIC           : %0.4f\n", aic))
      cat(sprintf("BIC           : %0.4f\n", bic))
      cat(sprintf("Observations  : %d\n", n_obs))
      cat("==========================================================\n")
      cat("Note: *** p<0.01, ** p<0.05, * p<0.1\n")
  
  # Returns an object for use if needed
  return(invisible(list(estimates = results_df, stats = c(LL=ll, AIC=aic, BIC=bic))))
  }

  cat(sprintf("\n[ESTIMATION] DCC-MIDAS-X for %s | K=%d, N_c=%d, Kc=%d, K_x=%d, Model=%s, Shock = %s\n", 
              target_asset, K, N_c, K_c, K_x, univ_model, shock_name))
  
  report_dcc_midas_results(
    opt_res = opt_results, 
    res_matrix = dccx_data$res_matrix, 
    x_matrix = dccx_data$X_mat
  )
  
  
  
  
  
  
  
  