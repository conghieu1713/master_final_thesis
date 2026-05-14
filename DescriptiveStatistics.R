###################################################################
###                                                             ###
#                    DESCRIPTIVE STATISTICS                       #
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
  library(moments)
  library(tseries)
  library(FinTS)
  library(dplyr)

#######======================================================#######
#######======================================================#######
#######======================================================#######


###===== ADD STARS FUNCTION =====###
  add_stars <- function(stat, p_val) {
    stars <- ifelse(p_val <= 0.01, "***",
                    ifelse(p_val <= 0.05, "**",
                           ifelse(p_val <= 0.1, "*", "")))
    return(paste0(sprintf("%.3f", stat), stars))
  }


#######======================================================#######
#######======================================================#######
#######======================================================#######
  

###===== DESCRIPTIVE STATISTICS FUNCTION =====###
  calc_stats <- function(x) {
    # Omit NA value
      x <- na.omit(x) 
    # Skip data with no more than 13 observations
      if(length(x) < 13) return(NULL)
    # Basic Descriptive Statistic
      mean_val <- mean(x)
      sd_val   <- sd(x)
      min_val  <- min(x)
      max_val  <- max(x)
      skew_val <- skewness(x)
      kurt_val <- kurtosis(x)
    # JB test
      jb <- jarque.bera.test(x)
      jb_str <- add_stars(jb$statistic, jb$p.value)
    # ADF test
      suppressWarnings({
        adf <- adf.test(x)
        adf_str <- add_stars(adf$statistic, adf$p.value)
    # PP test
      pp <- pp.test(x)
      pp_str <- add_stars(pp$statistic, pp$p.value)
    # KPSS test
      kpss <- kpss.test(x)
      kpss_str <- add_stars(kpss$statistic, kpss$p.value)
    })
    # ARCH(12) test
      arch <- ArchTest(x, lags = 12)
      arch_str <- add_stars(arch$statistic, arch$p.value)
    # LB(12) test (Ljung-Box)
      lb <- Box.test(x, lag = 12, type = "Ljung-Box")
      lb_str <- add_stars(lb$statistic, lb$p.value)
    # Results
    return(data.frame(
      Mean = sprintf("%.3f", mean_val),
      `Std. Dev` = sprintf("%.3f", sd_val),
      Min = sprintf("%.3f", min_val),
      Max = sprintf("%.3f", max_val),
      Skew = sprintf("%.3f", skew_val),
      Kurtosis = sprintf("%.3f", kurt_val),
      `JB test` = jb_str,
      `ADF test` = adf_str,
      `PP test` = pp_str,
      `KPSS test` = kpss_str,
      `ARCH (12)` = arch_str,
      `LB (12)` = lb_str,
      check.names = FALSE
    ))
  }
 
  
#######======================================================#######
#######======================================================#######
#######======================================================#######
   
 
###===== IMPORT DATA =====###
  # Folder Direction
    dir <- "D:\\Post-graduate programme\\Final Thesis\\Model"   
    setwd(dir)
  # High frequency data
    rate <- read_excel("official.xlsx", sheet = "rate", guess_max = 100000)
    rate$Date <- as.Date(rate$Date)
  # Low frequency data
    shocks <- read_excel("official.xlsx", sheet = "shocks", guess_max = 100000)
    shocks$Date <- as.Date(shocks$Date)


#######======================================================#######
#######======================================================#######
#######======================================================#######
  # Descriptive Statistics Result for High-Frequency Series
    # Setting up cols_to_test
      #cols_to_test <- c("VNIndex", "XAUUSD", "GC=F", "SJC", "BTC", "ETH", "BNB", "BRENT", "WTI")
      cols_to_test <- c("VNIndex", "XAUUSD", "GC=F", "SJC", "BTC", "ETH", "BNB")
    # Processing...
      results_list <- lapply(rate[cols_to_test], function(col) calc_stats(as.numeric(col)))
      final_table <- do.call(rbind, results_list)
    # Add Return Series Column
      final_table <- cbind(`Return Series` = rownames(final_table), final_table)
      rownames(final_table) <- NULL
    # Result
      print(final_table)
      write.csv(final_table, "DS Results\\Rates DS.csv", row.names = FALSE)
      
#######======================================================#######
#######======================================================#######
#######======================================================#######
  # Descriptive Statistics Result for Low-Frequency Series
    # Setting up cols_to_test
      #cols_to_test <- c("GPR", "GPRT", "GPRA", "GEPU_current", "GEPU_ppp", "WUI", "WPUI", "WSI", "WTUI")
      cols_to_test <- c("GPR", "GPRT", "GPRA", "GEPU_current", "GEPU_ppp", "WUI")
    # Processing...
      results_list <- lapply(shocks[cols_to_test], function(col) calc_stats(as.numeric(col)))
      final_table <- do.call(rbind, results_list)
    # Add Return Series Column
      final_table <- cbind(`Return Series` = rownames(final_table), final_table)
      rownames(final_table) <- NULL
    # Result
      print(final_table)
      write.csv(final_table, "DS Results\\Shocks DS.csv", row.names = FALSE)



      
      
      