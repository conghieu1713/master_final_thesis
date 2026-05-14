###################################################################
###                                                             ###
#                      CORRELATION MATRIX                         #
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


#######======================================================#######
#######======================================================#######
#######======================================================#######



###===== CORRELATION MATRIX FUNCTION =====###
  cor_with_stars <- function(data) {
    cols <- colnames(data)
    n <- length(cols)
    
    # Create a null matrix to store result
    res_matrix <- matrix("", n, n)
    rownames(res_matrix) <- cols
    colnames(res_matrix) <- cols
    
    for (i in 1:n) {
      for (j in 1:n) {
        if (i == j) {
          res_matrix[i, j] <- "1.000"
        } else {
          # Create a temporary data frame for col i,j and omit NA value
          pair_data <- na.omit(data.frame(x = data[[i]], y = data[[j]]))
          
          if (nrow(pair_data) > 2) {
            # Calculate Pearson and p.value
            test <- cor.test(pair_data$x, pair_data$y, method = "pearson")
            r <- test$estimate
            p <- test$p.value
            
            # Add stars
            stars <- ifelse(p < 0.01, "***", 
                            ifelse(p < 0.05, "**", 
                                   ifelse(p < 0.1, "*", "")))
            
            # Formatting
            res_matrix[i, j] <- paste0(sprintf("%.3f", r), stars)
          } else {
            # Return NA if data is not satisfactory
            res_matrix[i, j] <- "NA"
          }
        }
      }
    }
    return(as.data.frame(res_matrix))
  }
    
  
#######======================================================#######
#######======================================================#######
#######======================================================#######


  
###===== IMPORT DATA AND PRELIMINARY ANALYSIS =====###
  # Folder Direction
    dir <- "D:\\Post-graduate programme\\Final Thesis\\Model"   
    setwd(dir)
  # High frequency data
    cols_to_test <- c("VNIndex", "XAUUSD", "GC=F", "SJC", "BTC", "ETH", "BNB") 
    rate <- read_excel("official.xlsx", sheet = "rate", guess_max = 100000)
    rate <- rate[, cols_to_test]
    rate <- as.data.frame(sapply(rate, as.numeric))
  # Low frequency data
    cols_to_test <- c("GPR", "GPRT", "GPRA", "GEPU_current", "GEPU_ppp", "WUI")
    shocks <- read_excel("official.xlsx", sheet = "shocks", guess_max = 100000)
    shocks <- shocks[, cols_to_test]
    shocks <- as.data.frame(sapply(shocks, as.numeric))
  
  
#######======================================================#######
#######======================================================#######
#######======================================================#######
  

    
###===== CORRELATION MATRIX =====###
  # Correlation Matrix for High-Frequency Series
    matrix <- cor_with_stars(rate)
    print(matrix)
    write.csv(matrix, "DS Results\\CorMat HF.csv")
  # Correlation Matrix for Low-Frequency Series
    matrix <- cor_with_stars(shocks)
    print(matrix)
    write.csv(matrix, "DS Results\\CorMat LF.csv")

    

    