###################################################################
###                                                             ###
#                   GARCH-MIDAS SPECIFICATION                     #
###                                                             ###
###################################################################

# Author: Cong-Hieu, NGUYEN (524102110660)
# MSc. Candidate, University of Economics HCMC


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== RESET ALL =====###
  rm(list = ls(all.names = TRUE))


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== LIBRARY DECLARATION =====###
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(mfGARCH)
  library(tidyr)


#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== GARCH-MIDAS SPECIFICATION =====### 
# Folder Direction
  dir <- "D:\\Post-graduate programme\\Final Thesis\\Model"   
  setwd(dir)
  
# Asset List
  asset_list <- c("VNIndex", "XAUUSD", "GC=F", "SJC", "BTC", "ETH", "BNB", "BRENT", "WTI")

# Setting up
  lags_K <- c(6,12,24,36)
  lag_names <- c("Half year", "One year", "Two year", "Three year")
  
# Processing...
  all_results <- list()
  
  for (asset in asset_list) {
    cat("\n", asset, "processing...---\n")
    
    # Preliminary Analysis
      df_asset <- read_excel("official.xlsx", sheet = "rate") %>%
        select(Date, y = all_of(asset)) %>%
        inner_join(read_excel("official.xlsx", sheet = "rv") %>% select(Date, x = all_of(asset)), by = "Date") %>%
        rename(date = Date) %>%
        mutate(date = as.Date(date), year_month = floor_date(date, "month")) %>%
        drop_na() %>%
        arrange(date)
    
    # Fitting...
    for (asym_label in c("Without asymmetric", "With asymmetric")) {
      is_asym <- ifelse(asym_label == "With asymmetric", TRUE, FALSE)
      
      for (i in seq_along(lags_K)) {
        K_val <- lags_K[i]
        
        # Model Fitting...
          fit <- tryCatch({
            fit_mfgarch(data = df_asset, y = "y", x = "x",
                        low.freq = "year_month", var.ratio.freq = "year_month",
                        K = K_val, gamma = is_asym)
          }, error = function(e) {
            cat("     Error:", e$message, "\n") # Error
            return(NULL)
          })
        
        if (!is.null(fit)) {
          # Calculate LogLikelihood
            ll <- fit$llh
          # Calculate AIC and BIC
            k <- length(fit$par)
            aic <- 2 * k - 2 * ll
            bic <- fit$bic
          
          # Store results
          all_results[[length(all_results) + 1]] <- data.frame(
            Asset = asset,
            Specification = asym_label,
            Lag_Length = lag_names[i],
            LL = round(ll, 2),
            AIC = round(aic, 2),
            BIC = round(bic, 2)
          )
        }
      }
    }
  }

  # Result
    final_table <- bind_rows(all_results)
    print(final_table)
    write.csv(final_table, "DS Results\\GMspecification.csv", row.names = FALSE)


