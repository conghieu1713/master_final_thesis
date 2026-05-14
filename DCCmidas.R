###################################################################
###                                                             ###
#              DCC-MIDAS FITTING (using dccmidas)                 #
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


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== GET DATA FUNCTION =====###

file_path <- "D:\\Post-graduate programme\\Final Thesis\\Model\\official.xlsx" 
gold_assets <- c("XAUUSD", "GC=F", "SJC")
crypto_assets <- c("BTC", "ETH", "BNB")

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
  
  if (!is.null(fit)) {
    cat(">> [ESTIMATION] SUCCESS!.\n")
  }
  
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



###===== CHART FUNCTION =====###

plot_dcc_results <- function(estimation_result, filename = NULL) {
  # ====================================================
  # IMPORT DATA, using get_asset_data
  # ====================================================
    fit <- estimation_result$fit
    data_info <- estimation_result$data
    target_name <- data_info$target
    plot_dates <- data_info$dates
  
  # ====================================================
  # Calculate Volatility and Correlation
  # ====================================================
    vol_vn <- sqrt(fit$H_t[1, 1, ])
    vol_target <- sqrt(fit$H_t[2, 2, ])
    rho_sr <- fit$R_t[1, 2, ]
    rho_lr <- fit$R_t_bar[1, 2, ]
  # ====================================================
  # Create raw data.frame
  # ====================================================
    plot_df <- data.frame(
      Date = plot_dates,
      Vol_VN = as.numeric(vol_vn),
      Vol_Target = as.numeric(vol_target),
      SR_Corr = as.numeric(rho_sr),
      LR_Corr = as.numeric(rho_lr)
    )
  
  # ====================================================
  # CLEAR DATA - BURN IN
  # ====================================================
  # Automatically determine the starting point (remove burn-in lines)
  # I look for the first line where SR_Corr is not 0 or has a different initial value.
    burn_in_index <- which(plot_df$SR_Corr != 0 & !is.na(plot_df$SR_Corr))[1]
    if(is.na(burn_in_index)) burn_in_index <- 1
    plot_df_clean <- plot_df %>% slice(burn_in_index:nrow(plot_df))
  
  # ====================================================
  # THEME SETTING UP
  # ====================================================
  theme1 <- theme_bw() + 
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.major = element_line(linetype = "dotted", color = "darkgray"),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(), # Ẩn tên trục X
      axis.title.y = element_blank()  # Ẩn tên trục Y để giống hệt bài mẫu
    )
    
  # ====================================================
  # DRAW EACH PANEL (using ggplot2)
  # ====================================================
  # Panel A: VNIndex Volatility
    p1 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_VN)) +
      geom_line(color = "black", linewidth = 0.7) +
      labs(title = "VNIndex", y = "Conditional Volatility", x = "") +
      theme1
  # Panel B: Dynamic Correlation
    p2 <- ggplot(plot_df_clean, aes(x = Date)) +
      geom_line(aes(y = SR_Corr), color = "black", alpha = 0.7, linewidth = 0.6) +
      geom_line(aes(y = LR_Corr), color = "red", linewidth = 1) +
      labs(title = paste("VNIndex -", target_name), y = "Dynamic Correlation", x = "") +
      theme1
  # Panel C: Target Asset Volatility
    p3 <- ggplot(plot_df_clean, aes(x = Date, y = Vol_Target)) +
      geom_line(color = "black", linewidth = 0.7) +
      labs(title = target_name, y = "Conditional Volatility", x = "") +
      theme1
  
  # ====================================================
  # EXPORT
  # ====================================================
    des <- c(
      area(1,1),
      area(1,2),
      area(2,2)
    )
    final_plot <- p1 + p2 + p3 + 
      plot_layout(design = des) +
      plot_annotation(
        caption = "Fig. Conditional volatility, short- and long-run dynamic correlation.",
        theme = theme(plot.caption = element_text(hjust = 0.5))
      )
    
  # Print to RStudio screen
    print(final_plot)
    
  # If a file name is provided, proceed to save it as a high-quality PNG.
    if (!is.null(filename)) {
      ggsave(filename, plot = final_plot, width = 12, height = 8, dpi = 300)
      cat(sprintf(">> Direction: %s\n", filename))
    }
  
  return(final_plot)
}


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== RESULTS =====###

# ----------------------------------------------------------------- #
# ----- The lag values are used when fitting DCC-MIDAS-X model ---- #
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


  res <- run_dcc_estimation(
    target_asset = "SJC",
    file_path = file_path,
    K = 36, N_c = 21, K_c = 504,
    univ_model = "DAGM_noskew"
  )
  
  print(summary.dccmidas(res$fit))
  plot <- plot_dcc_results(estimation_result = res)
 
  
  
  
  
  
  
  